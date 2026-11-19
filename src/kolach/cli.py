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


SLURM_PROFILE = Path(__file__).parent / "workflows" / "profiles" / "slurm"
DEFAULT_SLURM_JOBS = 10


def _add_executor_options(parser, gpu=False):
    """Options for running the workflow's heavy steps as separate Slurm jobs."""
    group = parser.add_argument_group(
        "Slurm execution",
        "Run kolach in one small controller job; each heavy step (each tool or download) is submitted "
        "as its own Slurm job so they run in parallel.")
    group.add_argument(
        "--executor", choices=["local", "slurm"], default="local",
        help="Where workflow steps run: local (default) or slurm (one Slurm job per heavy step).")
    group.add_argument("--slurm-account", help="Slurm account charged for submitted jobs (required with --executor slurm).")
    group.add_argument("--slurm-partition", help="Slurm partition for submitted jobs (default: the cluster default).")
    group.add_argument(
        "--jobs", type=positive_int, default=None,
        help=f"Maximum number of Slurm jobs queued or running at once (default: {DEFAULT_SLURM_JOBS}).")
    group.add_argument(
        "--slurm-profile",
        help="Snakemake profile directory with per-step CPUs, memory and walltime, replacing kolach's bundled "
        f"one (copy {SLURM_PROFILE / 'config.yaml'} as a template).")
    if gpu:
        group.add_argument(
            "--slurm-gpu-partition",
            help="Partition for the DeepKOALA job when --device cuda (default: --slurm-partition).")


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
        choices=["eggnog", "kofam", "deepkoala", "dbcan", "graftm"],
        required=True,
        help="Databases to download. The KEGG Orthology hierarchy (ko00001.json) used by "
        "--add-pathways is always fetched alongside them. Specialisation needs dbcan and graftm.",
    )

    download_parser.add_argument(
        "--deepkoala-release",
        type=str,
        default="latest",
        help="DeepKOALA model release date (e.g. 202608 or latest, default: latest). "
        "'latest' is resolved against GenomeNet on each run, so newer releases are downloaded.",
    )

    _add_executor_options(download_parser)

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
        choices=["kofam", "eggnog", "deepkoala", "dbcan"],
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
        "--dbcan-keep-raw",
        action="store_true",
        help="Keep dbCAN native search files and run metadata (default: discard after integration).",
    )

    annotate_parser.add_argument(
        "--threads",
        type=positive_int,
        default=1,
        help="Number of CPU threads to use (default: 1). With --executor slurm, the CPUs the controller job "
        "uses for light steps (integration, pathways, specialisation); per-tool CPUs come from the Slurm "
        "profile's set-threads (a tool missing from a custom profile falls back to this value).",
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

    annotate_parser.add_argument(
        "--add-specialisations", action="store_true",
        help="Predict genome-level pathways and AnnoGuild guilds. Also runs dbCAN (CAZy substrates) "
        "and GraftM (marker refinement); needs `kolach download --databases dbcan graftm`.")
    annotate_group = annotate_parser.add_mutually_exclusive_group()
    annotate_group.add_argument(
        "--genome-id",
        help="Genome name for specialisation (default: the protein FASTA file name without extension).")
    annotate_group.add_argument("--gene-genome-map", help="TSV mapping gene_id to genome when the FASTA holds several genomes.")

    _add_executor_options(annotate_parser, gpu=True)

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
        "--dbcan-table", default=None,
        help="Normalized dbcan_overview.tsv (supplementary evidence, not a KO voter).",
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

    # Re-run specialisation on a finished `annotate --add-specialisations` output
    # directory, e.g. with other thresholds. It needs no databases.
    specialise_parser = subparsers.add_parser(
        "specialise",
        help="Re-evaluate pathways and AnnoGuild specialisations from an `annotate --add-specialisations` output directory.",
    )
    specialise_parser.add_argument(
        "--annotation-dir", required=True,
        help="Output directory of `kolach annotate --add-specialisations` (reads kolach_annotations.tsv "
        "and kolach_graftm_markers.tsv).")
    specialise_parser.add_argument(
        "--output-dir", required=True,
        help="Where to write results; must differ from --annotation-dir so the workflow's own reports stay intact.")
    specialise_group = specialise_parser.add_mutually_exclusive_group()
    specialise_group.add_argument(
        "--genome-id", help="All annotation rows belong to this genome (default: the genome recorded by the annotate run).")
    specialise_group.add_argument("--gene-genome-map", help="TSV with gene_id and genome columns.")
    specialise_parser.add_argument("--reaction-threshold", type=fraction, default=0.7,
                                   help="Minimum fraction of a route's reactions present (default: 0.7).")
    specialise_parser.add_argument("--ko-threshold", type=fraction, default=0.6,
                                   help="Minimum KO coverage of the best route (default: 0.6).")
    specialise_parser.add_argument("--reactions-file", help="Custom reaction TSV (pathway, reaction_id, definition).")
    specialise_parser.add_argument("--pathways-file", help="Custom EMERGE refined pathway TSV.")

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
        "dbcan_keep_raw": args.dbcan_keep_raw,
        "eggnog_mode": args.eggnog_mode,
        "eggnog_min_bitscore": args.eggnog_min_bitscore,
        "eggnog_max_evalue": args.eggnog_max_evalue,
        "eggnog_filter_multi": args.eggnog_filter_multi,
        "conflict_strategy": args.conflict_strategy,
        "add_pathways": args.add_pathways,
        "add_specialisations": args.add_specialisations,
    }
    if args.add_specialisations:
        # Specialisation needs dbCAN's CAZy evidence; add it rather than fail.
        if "dbcan" not in config["databases"]:
            config["databases"].append("dbcan")
        if not args.gene_genome_map:
            config["genome_id"] = args.genome_id or genome_name(args.protein_fasta)
    for key in ("eggnog_sensmode", "eggnog_temp_dir", "eggnog_dmnd_block_size", "eggnog_dmnd_index_chunks",
                "gene_genome_map"):
        value = getattr(args, key)
        if value is not None:
            config[key] = value
    if args.executor == "slurm" and args.device == "cuda" and "deepkoala" in config["databases"]:
        resources = {"gpu": 1}
        partition = args.slurm_gpu_partition or args.slurm_partition
        if partition:
            resources["slurm_partition"] = partition
        config["deepkoala_resources"] = resources
    return config


def load_profile(profile: Path) -> dict:
    """Settings of a Snakemake profile directory's config.yaml."""
    import yaml

    with open(profile / "config.yaml", encoding="utf-8") as handle:
        try:
            settings = yaml.safe_load(handle) or {}
        except yaml.YAMLError as exc:
            raise ValueError(f"invalid YAML in {profile / 'config.yaml'}: {exc}") from exc
    if not isinstance(settings, dict):
        raise ValueError(f"{profile / 'config.yaml'} must be a mapping of Snakemake options")
    return settings


def _resource_pairs(value) -> dict:
    """A profile resource setting (mapping, list or string of name=value) as a dict."""
    if not value:
        return {}
    if isinstance(value, dict):
        return {str(name): str(item) for name, item in value.items()}
    if isinstance(value, str):
        value = value.split()
    pairs = {}
    for item in value:
        name, sep, expression = str(item).partition("=")
        if not sep:
            raise ValueError(f"resource setting {item!r} is not name=value")
        pairs[name] = expression
    return pairs


def profile_default_resources(profile: Path) -> dict:
    """`default-resources` of a Snakemake profile as {name: expression}."""
    return _resource_pairs(load_profile(profile).get("default-resources"))


def check_slurm_extra(profile: Path) -> None:
    """Reject `slurm_extra` options the Slurm executor manages itself.

    The plugin raises for these (e.g. --comment, --qos, --mem) only once it tries
    to submit, and the controller then waits forever instead of exiting.
    """
    import re

    try:
        from snakemake_executor_plugin_slurm.validation import get_forbidden_slurm_options
    except ImportError:
        return
    settings = load_profile(profile)
    extras = [("default-resources", _resource_pairs(settings.get("default-resources")).get("slurm_extra"))]
    rule_resources = settings.get("set-resources") or {}
    if isinstance(rule_resources, dict):
        for rule, values in rule_resources.items():
            extras.append((f"set-resources {rule}", _resource_pairs(values).get("slurm_extra")))
    else:
        for item in ([rule_resources] if isinstance(rule_resources, str) else rule_resources):
            target, _, expression = str(item).partition("=")
            if target.endswith(":slurm_extra"):
                extras.append((f"set-resources {target.split(':')[0]}", expression))
    for where, extra in extras:
        if not extra:
            continue
        for pattern, option in get_forbidden_slurm_options().items():
            if re.search(pattern, extra):
                raise ValueError(
                    f"slurm_extra in {where} sets the {option} ({extra}), which the Slurm executor manages; "
                    "use the matching kolach option or profile resource instead")


def slurm_arguments(args) -> list:
    """Snakemake arguments that submit each non-local rule as its own Slurm job."""
    profile = Path(args.slurm_profile).expanduser().resolve() if args.slurm_profile else SLURM_PROFILE
    # A command-line --default-resources replaces the profile's whole map, so the
    # profile's own defaults are carried over. Snakemake evaluates resource values
    # as Python, so account and partition are quoted to stay literal strings
    # (e.g. a partition named 'all' or an account '0123').
    check_slurm_extra(profile)
    resources = profile_default_resources(profile)
    resources["slurm_account"] = repr(args.slurm_account)
    if args.slurm_partition:
        resources["slurm_partition"] = repr(args.slurm_partition)
    return [
        "--executor", "slurm",
        "--jobs", str(args.jobs or DEFAULT_SLURM_JOBS),
        "--workflow-profile", str(profile),
        "--default-resources", *(f"{name}={value}" for name, value in resources.items()),
    ]


FASTA_SUFFIXES = (".faa", ".fa", ".fasta", ".fas", ".pep", ".aa", ".prot", ".protein")


def genome_name(protein_fasta) -> str:
    """Genome name from the FASTA file name: drop .gz, then one FASTA suffix."""
    name = Path(protein_fasta).name
    if name.lower().endswith((".gz", ".gzip")):
        name = name.rsplit(".", 1)[0]
    stem, dot, suffix = name.rpartition(".")
    if dot and stem and f".{suffix.lower()}" in FASTA_SUFFIXES:
        name = stem
    return name


def check_gene_genome_map(path, gene_ids) -> None:
    """Fail before any search if the map does not cover every FASTA gene."""
    from kolach.metabolism.evaluate import read_tsv

    columns, rows = read_tsv(path)
    if not {"gene_id", "genome"} <= set(columns):
        raise ValueError(f"{path} needs gene_id and genome columns")
    mapped = set()
    for row in rows:
        # Same per-row rules as specialisation, so a bad map fails up front.
        if not row["gene_id"] or not row["genome"] or row["gene_id"] in mapped:
            raise ValueError(f"Empty or duplicate gene ID in {path}: {row['gene_id']!r}")
        mapped.add(row["gene_id"])
    missing = [gene for gene in gene_ids if gene not in mapped]
    if missing:
        raise ValueError(f"{len(missing)} FASTA gene(s) have no genome in {path}, e.g. {missing[:3]}")


def check_annotation_dependencies(databases, specialisations=False) -> None:
    """Fail before launching the workflow if a selected method cannot run."""
    for db in databases:
        importlib.import_module(f"kolach.methods.{db}").check_dependencies()
    if specialisations:
        importlib.import_module("kolach.methods.graftm").check_dependencies()


def run_workflow(snakefile: str, config: dict, cores: int, slurm_args=None, config_dir=None) -> int:
    """Run a bundled Snakefile with `config` passed as a JSON config file.

    A config file (rather than `--config key=value`) keeps values typed and stops
    Snakemake from YAML-parsing paths and strings into other types.

    With `slurm_args`, heavy rules are submitted as Slurm jobs and `cores` only
    covers the local rules. Those jobs re-read the config file on other nodes, so
    it is written to `config_dir` (a shared directory) rather than /tmp, under a
    per-run name so another invocation cannot overwrite it mid-run.
    """
    if not shutil.which("snakemake"):
        raise FileNotFoundError(
            "The 'snakemake' executable was not found in PATH. "
            "Please activate the kolach conda environment."
        )
    snakefile_path = Path(__file__).parent / "workflows" / snakefile
    if slurm_args:
        config_dir = Path(config_dir).expanduser().resolve()
        config_dir.mkdir(parents=True, exist_ok=True)
        fd, config_path = tempfile.mkstemp(prefix=f"{Path(snakefile).stem}_config_", suffix=".json",
                                           dir=config_dir)
        os.chmod(config_path, 0o644)
    else:
        fd, config_path = tempfile.mkstemp(prefix="kolach_config_", suffix=".json")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(config, handle, indent=2)
        cmd = ["snakemake", "-s", str(snakefile_path)]
        if not slurm_args:
            cmd += ["--cores", str(cores), "--configfile", config_path]
            return subprocess.run(cmd).returncode
        cmd += ["--local-cores", str(cores), *slurm_args, "--configfile", config_path]
        return run_controller(cmd)
    finally:
        Path(config_path).unlink(missing_ok=True)


def run_controller(cmd) -> int:
    """Run a Slurm-executor Snakemake, cancelling its jobs if kolach is terminated.

    On SIGTERM (scancel, or the controller job reaching its walltime) Snakemake
    would keep waiting for its submitted jobs. SIGINT makes it cancel them and
    release its lock, so a later run of the same command can resume.
    """
    import signal

    process = subprocess.Popen(cmd)

    def interrupt(signum, frame):
        process.send_signal(signal.SIGINT)

    previous = signal.signal(signal.SIGTERM, interrupt)
    try:
        while True:
            try:
                return process.wait()
            except KeyboardInterrupt:
                # Ctrl-C reaches Snakemake too; keep waiting for its cleanup.
                continue
    finally:
        signal.signal(signal.SIGTERM, previous)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "integrate" and args.add_pathways and not args.database_dir:
        parser.error("integrate --add-pathways requires --database-dir")

    if args.command == "annotate":
        if not args.add_specialisations and (args.genome_id or args.gene_genome_map):
            parser.error("--genome-id and --gene-genome-map require --add-specialisations")
        if args.genome_id is not None and not args.genome_id.strip():
            parser.error("--genome-id cannot be empty")

    slurm_args = None
    if args.command in ("download", "annotate"):
        gpu_partition = getattr(args, "slurm_gpu_partition", None)
        slurm_only = [flag for flag, value in (
            ("--slurm-account", args.slurm_account), ("--slurm-partition", args.slurm_partition),
            ("--slurm-profile", args.slurm_profile), ("--jobs", args.jobs),
            ("--slurm-gpu-partition", gpu_partition)) if value]
        if args.executor == "slurm":
            if not args.slurm_account:
                parser.error("--executor slurm requires --slurm-account")
            if args.slurm_profile and not (Path(args.slurm_profile).expanduser() / "config.yaml").is_file():
                parser.error(f"--slurm-profile must be a directory containing config.yaml: {args.slurm_profile}")
            if gpu_partition and not (args.device == "cuda" and "deepkoala" in args.databases):
                parser.error("--slurm-gpu-partition requires --device cuda and deepkoala in --databases")
            try:
                slurm_args = slurm_arguments(args)
            except (OSError, ValueError, AttributeError) as exc:
                parser.error(f"invalid Slurm profile: {exc}")
        elif slurm_only:
            parser.error(f"{', '.join(slurm_only)} require(s) --executor slurm")

    if args.command == "download":
        config = download_config(args)
        if slurm_args and "deepkoala" in config["databases"]:
            from kolach.methods.deepkoala import resolve_release

            try:
                config["deepkoala_release"] = resolve_release(args.deepkoala_release)
            except (OSError, ValueError) as exc:
                parser.exit(1, f"kolach download: error: {exc}\n")
            config["deepkoala_release_resolved"] = True
        returncode = run_workflow("download.smk", config, cores=1, slurm_args=slurm_args,
                                  config_dir=Path(args.database_dir) / ".kolach")
        if returncode != 0:
            sys.exit(returncode)

    elif args.command == "annotate":
        from kolach.integrate import read_fasta_ids

        config = annotate_config(args)
        try:
            check_annotation_dependencies(config["databases"], args.add_specialisations)
            # Duplicate IDs would only be rejected at the integrate step, after
            # every annotation method has already run.
            gene_ids = read_fasta_ids(args.protein_fasta)
            if args.gene_genome_map:
                check_gene_genome_map(args.gene_genome_map, gene_ids)
        except (OSError, ImportError, ValueError) as exc:
            parser.exit(1, f"kolach annotate: error: {exc}\n")
        returncode = run_workflow("annotate.smk", config, cores=args.threads, slurm_args=slurm_args,
                                  config_dir=Path(args.output_dir) / ".kolach")
        if returncode != 0:
            sys.exit(returncode)

    elif args.command == "integrate":
        from kolach.integrate import integrate_annotations

        integrate_annotations(
            protein_fasta=args.protein_fasta,
            kofam_tsv=args.kofam_table,
            deepkoala_tsv=args.deepkoala_table,
            eggnog_tsv=args.eggnog_table,
            dbcan_tsv=args.dbcan_table,
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

    elif args.command == "specialise":
        from kolach.metabolism.specialise import specialise

        directory = Path(args.annotation_dir)
        annotations = directory / "kolach_annotations.tsv"
        markers = directory / "kolach_graftm_markers.tsv"
        if Path(args.output_dir).resolve() == directory.resolve():
            parser.error("--output-dir must differ from --annotation-dir")
        try:
            for path in (annotations, markers):
                if not path.is_file():
                    raise FileNotFoundError(f"{path} not found; run `kolach annotate --add-specialisations` first")
            genome_id, genome_map = args.genome_id, args.gene_genome_map
            metadata = directory / "kolach_specialisation_metadata.json"
            recorded = json.loads(metadata.read_text(encoding="utf-8")) if metadata.is_file() else {}
            if genome_id is None and genome_map is None:
                # Reuse the grouping the annotate run recorded.
                genome_id = recorded.get("genome_id")
                genome_map = (recorded.get("inputs", {}).get("genome_map") or {}).get("path")
                if genome_id is None and genome_map is None:
                    raise ValueError("No recorded genome; give --genome-id or --gene-genome-map")
            paths = specialise(
                annotations, args.output_dir,
                genome_id=genome_id, gene_genome_map=genome_map,
                markers=markers, cazy_from_dbcan=True,
                reaction_threshold=args.reaction_threshold, ko_threshold=args.ko_threshold,
                reactions_file=args.reactions_file, pathways_file=args.pathways_file,
                graftm_bundle=recorded.get("graftm_bundle"),
            )
        except (ValueError, OSError, KeyError) as error:
            parser.exit(1, f"kolach specialise: error: {error}\n")
        print(f"[kolach] Wrote specialisations to {paths['specialisations']}")

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
