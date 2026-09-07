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
        choices=["eggnog"],
        required=True,
        help="Databases to download.",
    )
    # ------------------------------------------------------
    # end of download subcommand
    # ------------------------------------------------------
    
    # ------------------------------------------------------
    # annotate subcommand
    # ------------------------------------------------------
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate input proteins with KEGG Orthology identifiers.",
    )

    annotate_parser.add_argument(
        "--protein-fasta",
        type=str,
        required=True,
        help="Path to the input protein FASTA file (or stdout).",
    )

    annotate_parser.add_argument(
        "--database-dir",
        type=str,
        required=True,
        help="Path to the directory containing databases.",
    )

    # annotate_parser.add_argument(
    #     "--hmm-evalue",
    #     type=float,
    #     default=1e-3,
    #     help="E-value threshold for HMMER (default: 1e-3).",
    # )

    annotate_parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (or stdout if not provided).",
    )
    # ------------------------------------------------------
    # end of annotate subcommand
    # ------------------------------------------------------


    args = parser.parse_args()

    if args.command == "download":
        snakefile_path = Path(__file__).parent / "workflows" / "download.smk"

        # Format config arguments for Snakemake CLI
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


if __name__ == "__main__":
    main()
