import argparse
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import subprocess
import sys


def _package_version() -> str:
    try:
        return version("kolach")
    except PackageNotFoundError:
        return "0.0.0"


def positive_float(value: str) -> float:
    try:
        val = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value: '{value}'")
    if val <= 0:
        raise argparse.ArgumentTypeError(f"Value must be positive (> 0), got {value}")
    return val


def positive_int(value: str) -> int:
    try:
        val = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: '{value}'")
    if val <= 0:
        raise argparse.ArgumentTypeError(f"Value must be a positive integer (> 0), got {value}")
    return val


def main():
    parser = argparse.ArgumentParser(
        prog="kolach",
        description="Assign KEGG Orthology identifiers to proteins.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"kolach {_package_version()}",
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
        help="Databases to download.",
    )

    download_parser.add_argument(
        "--deepkoala-release",
        type=str,
        default="latest",
        help="DeepKOALA model release date (e.g. 202608 or latest, default: latest).",
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
        type=int,
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
        type=int,
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

    annotate_parser.add_argument(
        "--eggnog-min-bitscore",
        type=float,
        default=60.0,
        help="Minimum eggNOG bit score to retain KO assignment during integration (default: 60.0).",
    )

    annotate_parser.add_argument(
        "--eggnog-max-evalue",
        type=float,
        default=1e-5,
        help="Maximum eggNOG e-value to retain KO assignment during integration (default: 1e-5).",
    )

    annotate_parser.add_argument(
        "--eggnog-filter-multi",
        choices=["disambiguate", "none"],
        default="disambiguate",
        help="eggNOG multi-KO filtering strategy: disambiguate against other tools (default) or none.",
    )

    annotate_parser.add_argument(
        "--conflict-strategy",
        choices=["multiple", "priority"],
        default="multiple",
        help=(
            "Consensus conflict strategy for disjoint calls: multiple (default: sets accepted_ko and ko to '-', "
            "records conflicting alternatives in alternative_kos; recommended for downstream pathway tools to avoid "
            "false multifunctional enzyme inference) or priority (selects top method in hierarchy: kofam > eggnog > deepkoala)."
        ),
    )

    annotate_parser.add_argument(
        "--add-pathways",
        action="store_true",
        default=False,
        help="Annotate accepted KOs with KEGG categories, subcategories, and pathways in kolach_annotations.tsv.",
    )

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

    integrate_parser.add_argument(
        "--eggnog-min-bitscore",
        type=float,
        default=60.0,
        help="Minimum eggNOG bit score for orthology retention (default: 60.0).",
    )

    integrate_parser.add_argument(
        "--eggnog-max-evalue",
        type=float,
        default=1e-5,
        help="Maximum eggNOG e-value for orthology retention (default: 1e-5).",
    )

    integrate_parser.add_argument(
        "--eggnog-filter-multi",
        choices=["disambiguate", "none"],
        default="disambiguate",
        help="eggNOG multi-KO disambiguation strategy: disambiguate against other tools (default) or none.",
    )

    integrate_parser.add_argument(
        "--conflict-strategy",
        choices=["multiple", "priority"],
        default="multiple",
        help=(
            "Consensus conflict strategy for disjoint calls: multiple (default: sets accepted_ko and ko to '-', "
            "records conflicting alternatives in alternative_kos; recommended for downstream pathway tools to avoid "
            "false multifunctional enzyme inference) or priority (selects top method in hierarchy: kofam > eggnog > deepkoala)."
        ),
    )

    integrate_parser.add_argument(
        "--add-pathways",
        action="store_true",
        default=False,
        help="Annotate accepted KOs with KEGG categories, subcategories, and pathways in the output table.",
    )

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

    args = parser.parse_args()

    if args.command == "download":
        snakefile_path = Path(__file__).parent / "workflows" / "download.smk"
        db_list_repr = f"[{','.join(repr(db) for db in args.databases)}]"
        config_args = [
            f"database_dir={args.database_dir}",
            f"databases={db_list_repr}",
            f"deepkoala_release={args.deepkoala_release}",
        ]
        cmd = [
            "snakemake",
            "-s", str(snakefile_path),
            "--cores", "1",
            "--config", *config_args,
        ]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)

    elif args.command == "annotate":
        snakefile_path = Path(__file__).parent / "workflows" / "annotate.smk"
        db_list_repr = f"[{','.join(repr(db) for db in args.databases)}]"
        config_args = [
            f"database_dir={args.database_dir}",
            f"protein_fasta={args.protein_fasta}",
            f"output_dir={args.output_dir}",
            f"databases={db_list_repr}",
            f"threads={args.threads}",
            f"skip_bitscore_heuristic={args.skip_bitscore_heuristic}",
            f"no_hmmer_prefiltering={args.no_hmmer_prefiltering}",
            f"deepkoala_model={args.deepkoala_model}",
            f"deepkoala_release={args.deepkoala_release}",
            f"device={args.device}",
            f"batch_size={args.batch_size}",
            f"detail={args.detail}",
            f"eggnog_mode={args.eggnog_mode}",
            f"eggnog_min_bitscore={args.eggnog_min_bitscore}",
            f"eggnog_max_evalue={args.eggnog_max_evalue}",
            f"eggnog_filter_multi={args.eggnog_filter_multi}",
            f"conflict_strategy={args.conflict_strategy}",
            f"add_pathways={args.add_pathways}",
        ]
        if args.eggnog_sensmode is not None:
            config_args.append(f"eggnog_sensmode={args.eggnog_sensmode}")
        if args.eggnog_temp_dir is not None:
            config_args.append(f"eggnog_temp_dir={args.eggnog_temp_dir}")
        if args.eggnog_dmnd_block_size is not None:
            config_args.append(f"eggnog_dmnd_block_size={args.eggnog_dmnd_block_size}")
        if args.eggnog_dmnd_index_chunks is not None:
            config_args.append(f"eggnog_dmnd_index_chunks={args.eggnog_dmnd_index_chunks}")
        cmd = [
            "snakemake",
            "-s", str(snakefile_path),
            "--cores", str(args.threads),
            "--config", *config_args,
        ]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)

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

