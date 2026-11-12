import argparse
import importlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from kolach import package_version


def positive_float(value: str) -> float:
    try:
        val = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value: '{value}'")
    if not math.isfinite(val) or val <= 0:
        raise argparse.ArgumentTypeError(f"Value must be positive (> 0), got {value}")
    return val


def nonnegative_float(value: str) -> float:
    try:
        val = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value: '{value}'")
    if not math.isfinite(val) or val < 0:
        raise argparse.ArgumentTypeError(f"Value must be finite and nonnegative, got {value}")
    return val


def fraction(value: str) -> float:
    val = nonnegative_float(value)
    if val > 1:
        raise argparse.ArgumentTypeError(f"Value must be within [0, 1], got {value}")
    return val


def positive_int(value: str) -> int:
    try:
        val = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: '{value}'")
    if val <= 0:
        raise argparse.ArgumentTypeError(f"Value must be a positive integer (> 0), got {value}")
    return val


CONFLICT_STRATEGY_HELP = (
    "Consensus conflict strategy for disjoint calls: multiple (default: sets accepted_ko and ko to '-', "
    "records conflicting alternatives in alternative_kos; recommended for downstream pathway tools to avoid "
    "false multifunctional enzyme inference) or priority (selects top method in hierarchy: kofam > eggnog > deepkoala)."
)


def _add_integration_options(parser, heuristic_note=""):
    """Options shared by `annotate` and `integrate` that control consensus integration."""
    parser.add_argument(
        "--heuristic-bitscore-fraction",
        type=fraction,
        default=0.75,
        help=f"Fraction of the KOfam profile threshold used to rescue below-threshold hits{heuristic_note} "
        "(default: 0.75).",
    )
    parser.add_argument(
        "--heuristic-e-value",
        type=nonnegative_float,
        default=1e-5,
        help=f"Maximum E-value for KOfam hit rescue consideration{heuristic_note} (default: 1e-5).",
    )
    parser.add_argument(
        "--eggnog-min-bitscore",
        type=nonnegative_float,
        default=60.0,
        help="Minimum eggNOG bit score to retain KO assignment during integration (default: 60.0).",
    )
    parser.add_argument(
        "--eggnog-max-evalue",
        type=nonnegative_float,
        default=1e-5,
        help="Maximum eggNOG e-value to retain KO assignment during integration (default: 1e-5).",
    )
    parser.add_argument(
        "--eggnog-filter-multi",
        choices=["disambiguate", "none"],
        default="disambiguate",
        help="eggNOG multi-KO filtering strategy: disambiguate against other tools (default) or none.",
    )
    parser.add_argument(
        "--conflict-strategy",
        choices=["multiple", "priority"],
        default="multiple",
        help=CONFLICT_STRATEGY_HELP,
    )
    parser.add_argument(
        "--add-pathways",
        action="store_true",
        default=False,
        help="Annotate accepted KOs with KEGG categories, subcategories, and pathways in the output table.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="kolach",
        description="Assign KEGG Orthology identifiers to proteins.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"kolach {package_version()}",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # ------------------------------------------------------
    # Download subcommand
    # ------------------------------------------------------
    download_parser = subparsers.add_parser(
        "download",
        help="Download databases.",
    )

    download_parser.add_argument(
        "--database-dir",
        required=True,
        help="Directory where databases will be stored.",
    )

    download_parser.add_argument(
        "--databases",
        nargs="+",
        choices=["eggnog", "kofam", "deepkoala"],
        required=True,
        help="Databases to download. The KEGG Orthology hierarchy (ko00001.json) used by "
        "--add-pathways is always fetched alongside them.",
    )

    download_parser.add_argument(
        "--deepkoala-release",
        type=str,
        default="latest",
        help="DeepKOALA model release date (e.g. 202608 or latest, default: latest). "
        "'latest' is resolved against GenomeNet on each run, so newer releases are downloaded.",
    )

    # ------------------------------------------------------
    # Annotate subcommand
    # ------------------------------------------------------
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate input proteins with KEGG Orthology identifiers.",
    )

    annotate_parser.add_argument(
        "--protein-fasta",
        type=str,
        required=True,
        help="Path to the input protein FASTA file.",
    )

    annotate_parser.add_argument(
        "--database-dir",
        type=str,
        required=True,
        help="Path to the directory containing databases.",
    )

    annotate_parser.add_argument(
        "--databases",
        nargs="+",
        choices=["kofam", "eggnog", "deepkoala"],
        default=["kofam"],
        help="Databases to use for annotation (default: kofam).",
    )

    annotate_parser.add_argument(
        "--output-dir",
        type=str,
        default="kolach_output",
        help="Directory where output files will be written (default: kolach_output).",
    )

    annotate_parser.add_argument(
        "--threads",
        type=positive_int,
        default=1,
        help="Number of CPU threads to use (default: 1).",
    )

    annotate_parser.add_argument(
        "--skip-bitscore-heuristic",
        action="store_true",
        help="Skip the 75 percent bitscore relaxation heuristic for KOfam.",
    )

    annotate_parser.add_argument(
        "--no-hmmer-prefiltering",
        action="store_true",
        help="Disable HMMER e-value prefiltering for large datasets.",
    )

    annotate_parser.add_argument(
        "--deepkoala-model",
        choices=["full", "frag"],
        default="full",
        help="DeepKOALA model variant (full or frag, default: full).",
    )

    annotate_parser.add_argument(
        "--deepkoala-release",
        type=str,
        default="latest",
        help="DeepKOALA model release date (e.g. 202608 or latest, default: latest).",
    )

    annotate_parser.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda"],
        default="auto",
        help="Inference device for deep learning methods (default: auto).",
    )

    annotate_parser.add_argument(
        "--batch-size",
        type=positive_int,
        default=64,
        help="Batch size for DeepKOALA inference (default: 64).",
    )

    annotate_parser.add_argument(
        "--detail",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include detailed probabilities, thresholds, and marks in DeepKOALA output (default: True).",
    )

    annotate_parser.add_argument(
        "--eggnog-mode",
        choices=["diamond", "mmseqs"],
        default="diamond",
        help="Search mode for eggNOG-mapper (diamond or mmseqs, default: diamond).",
    )

    annotate_parser.add_argument(
        "--eggnog-sensmode",
        choices=["default", "fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"],
        default=None,
        help="Diamond sensitivity mode for eggNOG-mapper. If unspecified, eggNOG-mapper's default (sensitive, iterative) is used. Selecting 'default' explicitly requests DIAMOND's default sensitivity mode.",
    )

    annotate_parser.add_argument(
        "--eggnog-temp-dir",
        type=str,
        default=None,
        help="Base directory for eggNOG temporary files and DIAMOND scratch space (default: system temporary directory).",
    )

    annotate_parser.add_argument(
        "--eggnog-dmnd-block-size",
        type=positive_float,
        default=None,
        help="DIAMOND block size in billions of sequence letters (--dmnd_block_size). If unspecified, upstream automatic tuning based on host RAM is used. Note: host RAM tuning may exceed allocated memory in Slurm jobs.",
    )

    annotate_parser.add_argument(
        "--eggnog-dmnd-index-chunks",
        type=positive_int,
        default=None,
        help="Number of chunks for processing DIAMOND seed index (--dmnd_index_chunks). If unspecified, upstream automatic tuning is used.",
    )

    _add_integration_options(annotate_parser)

    # ------------------------------------------------------
    # Integrate subcommand
    # ------------------------------------------------------
    integrate_parser = subparsers.add_parser(
        "integrate",
        help="Integrate annotation tables into a unified consensus table.",
    )

    integrate_parser.add_argument(
        "--protein-fasta",
        type=str,
        default=None,
        help="Path to the input protein FASTA file (preserves gene order and full gene inventory).",
    )

    integrate_parser.add_argument(
        "--kofam-table",
        type=str,
        default=None,
        help="Path to KOfam annotations TSV.",
    )

    integrate_parser.add_argument(
        "--deepkoala-table",
        type=str,
        default=None,
        help="Path to DeepKOALA annotations TSV.",
    )

    integrate_parser.add_argument(
        "--eggnog-table",
        type=str,
        default=None,
        help="Path to eggNOG annotations TSV.",
    )

    integrate_parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path for integrated output TSV.",
    )

    integrate_parser.add_argument(
        "--evidence-file",
        type=str,
        default=None,
        help="Path for long-form evidence output TSV (default: kolach_evidence.tsv alongside output-file).",
    )

    integrate_parser.add_argument(
        "--database-dir",
        type=str,
        default=None,
        help="Path to database directory (to load KO definitions from ko_list).",
    )

    _add_integration_options(integrate_parser, heuristic_note=", reported in the evidence table")

    # ------------------------------------------------------
    # Pathway subcommand
    # ------------------------------------------------------
    pathway_parser = subparsers.add_parser(
        "pathway",
        help="Annotate KO identifiers in an annotation table with KEGG categories, subcategories, and pathways.",
    )

    pathway_parser.add_argument(
        "--annotation-table",
        type=str,
        required=True,
        help="Path to the input annotations TSV table (e.g. kolach_annotations.tsv).",
    )

    pathway_parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path for output TSV table with KEGG pathway annotations.",
    )

    pathway_parser.add_argument(
        "--database-dir",
        type=str,
        required=True,
        help="Directory where databases are stored (contains or will download ko00001.json).",
    )

    pathway_parser.add_argument(
        "--ko-column",
        type=str,
        default="accepted_ko",
        help="Column name containing KO identifiers to annotate (default: 'accepted_ko').",
    )

    return parser


def download_config(args) -> dict:
    """Workflow config for `kolach download`, with native JSON types."""
    return {
        "database_dir": args.database_dir,
        "databases": list(args.databases),
        "deepkoala_release": args.deepkoala_release,
    }


def annotate_config(args) -> dict:
    """Workflow config for `kolach annotate`, with native JSON types.

    Optional eggNOG settings are omitted when unset so the rules fall back to
    upstream defaults.
    """
    config = {
        "database_dir": args.database_dir,
        "protein_fasta": args.protein_fasta,
        "output_dir": args.output_dir,
        "databases": list(args.databases),
        "threads": args.threads,
        "skip_bitscore_heuristic": args.skip_bitscore_heuristic,
        "no_hmmer_prefiltering": args.no_hmmer_prefiltering,
        "heuristic_bitscore_fraction": args.heuristic_bitscore_fraction,
        "heuristic_e_value": args.heuristic_e_value,
        "deepkoala_model": args.deepkoala_model,
        "deepkoala_release": args.deepkoala_release,
        "device": args.device,
        "batch_size": args.batch_size,
        "detail": args.detail,
        "eggnog_mode": args.eggnog_mode,
        "eggnog_min_bitscore": args.eggnog_min_bitscore,
        "eggnog_max_evalue": args.eggnog_max_evalue,
        "eggnog_filter_multi": args.eggnog_filter_multi,
        "conflict_strategy": args.conflict_strategy,
        "add_pathways": args.add_pathways,
    }
    for key in ("eggnog_sensmode", "eggnog_temp_dir", "eggnog_dmnd_block_size", "eggnog_dmnd_index_chunks"):
        value = getattr(args, key)
        if value is not None:
            config[key] = value
    return config


def check_annotation_dependencies(databases) -> None:
    """Fail before launching the workflow if a selected method cannot run."""
    for db in databases:
        importlib.import_module(f"kolach.methods.{db}").check_dependencies()


def run_workflow(snakefile: str, config: dict, cores: int) -> int:
    """Run a bundled Snakefile with `config` passed as a JSON config file.

    A config file (rather than `--config key=value`) keeps values typed and stops
    Snakemake from YAML-parsing paths and strings into other types.
    """
    if not shutil.which("snakemake"):
        raise FileNotFoundError(
            "The 'snakemake' executable was not found in PATH. "
            "Please activate the kolach conda environment."
        )
    snakefile_path = Path(__file__).parent / "workflows" / snakefile
    fd, config_path = tempfile.mkstemp(prefix="kolach_config_", suffix=".json")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(config, handle, indent=2)
        cmd = [
            "snakemake",
            "-s", str(snakefile_path),
            "--cores", str(cores),
            "--configfile", config_path,
        ]
        return subprocess.run(cmd).returncode
    finally:
        Path(config_path).unlink(missing_ok=True)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "integrate" and args.add_pathways and not args.database_dir:
        parser.error("integrate --add-pathways requires --database-dir")

    if args.command == "download":
        returncode = run_workflow("download.smk", download_config(args), cores=1)
        if returncode != 0:
            sys.exit(returncode)

    elif args.command == "annotate":
        try:
            check_annotation_dependencies(args.databases)
        except (FileNotFoundError, ImportError) as exc:
            parser.exit(1, f"kolach annotate: error: {exc}\n")
        returncode = run_workflow("annotate.smk", annotate_config(args), cores=args.threads)
        if returncode != 0:
            sys.exit(returncode)

    elif args.command == "integrate":
        from kolach.integrate import integrate_annotations

        integrate_annotations(
            protein_fasta=args.protein_fasta,
            kofam_tsv=args.kofam_table,
            deepkoala_tsv=args.deepkoala_table,
            eggnog_tsv=args.eggnog_table,
            output_tsv=args.output_file,
            evidence_tsv=args.evidence_file,
            database_dir=args.database_dir,
            eggnog_min_bitscore=args.eggnog_min_bitscore,
            eggnog_max_evalue=args.eggnog_max_evalue,
            eggnog_filter_multi=args.eggnog_filter_multi,
            conflict_strategy=args.conflict_strategy,
            heuristic_bitscore_fraction=args.heuristic_bitscore_fraction,
            heuristic_e_value=args.heuristic_e_value,
        )

        if args.add_pathways:
            from kolach.pathways import annotate_pathways

            annotate_pathways(
                annotation_table=args.output_file,
                output_tsv=args.output_file,
                database_dir=args.database_dir,
            )

    elif args.command == "pathway":
        from kolach.pathways import annotate_pathways

        annotate_pathways(
            annotation_table=args.annotation_table,
            output_tsv=args.output_file,
            database_dir=args.database_dir,
            ko_col=args.ko_column,
        )


if __name__ == "__main__":
    main()
