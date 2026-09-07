import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="kolach",
        description="Assign KEGG Orthology identifiers to proteins.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version="kolach 0.0.0",
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
        choices=["eggnog", "kofam"],
        required=True,
        help="Databases to download.",
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
        choices=["kofam", "eggnog"],
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
        help="Skip the 75% bitscore relaxation heuristic for KOfam.",
    )

    annotate_parser.add_argument(
        "--no-hmmer-prefiltering",
        action="store_true",
        help="Disable HMMER e-value prefiltering for large datasets.",
    )

    args = parser.parse_args()

    if args.command == "download":
        snakefile_path = Path(__file__).parent / "workflows" / "download.smk"
        db_list_repr = f"[{','.join(repr(db) for db in args.databases)}]"
        config_args = [
            f"database_dir={args.database_dir}",
            f"databases={db_list_repr}",
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
        ]
        cmd = [
            "snakemake",
            "-s", str(snakefile_path),
            "--cores", str(args.threads),
            "--config", *config_args,
        ]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)


if __name__ == "__main__":
    main()
