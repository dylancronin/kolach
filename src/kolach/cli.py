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
        help="Skip the 75% bitscore relaxation heuristic for KOfam.",
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
        action="store_true",
        help="Include detailed probabilities, thresholds, and marks in DeepKOALA output.",
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
        default="default",
        help="Diamond sensitivity mode for eggNOG-mapper (default: default).",
    )

    annotate_parser.add_argument(
        "--eggnog-dbmem",
        action="store_true",
        help="Load eggNOG diamond database into memory for faster execution.",
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
            f"eggnog_sensmode={args.eggnog_sensmode}",
            f"eggnog_dbmem={args.eggnog_dbmem}",
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
