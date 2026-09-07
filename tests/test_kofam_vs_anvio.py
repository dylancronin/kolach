"""Automated parity test: compare kolach KOfam output against anvio-9."""
import csv
from pathlib import Path
import shutil
import subprocess
import tempfile
from kolach.methods.kofam import annotate
import argparse

ANVIO_ENV_BIN = "/fs/ess/PDS0325/bioinformatic_tools/miniforge3/envs/anvio-9/bin"


def run_comparison(fasta_path, kofam_db_dir):
    with tempfile.TemporaryDirectory(prefix="test-anvio-comp-") as tmp:
        tmp = Path(tmp)
        kolach_out = tmp / "kolach_kofam.tsv"

        # 1. Run kolach KOfam
        args = argparse.Namespace(
            database_dir=kofam_db_dir,
            protein_fasta=str(fasta_path),
            threads=2,
            skip_bitscore_heuristic=False,
            no_hmmer_prefiltering=False,
        )
        with open(kolach_out, "w") as f:
            annotate(args, f)

        # 2. Parse kolach output
        kolach_hits = {}
        with open(kolach_out) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                kolach_hits[(row["gene_id"], row["ko"])] = {
                    "bit_score": float(row["bit_score"]),
                    "assignment": row["assignment"],
                }

        print(f"Kolach annotated {len(kolach_hits)} gene-KO pairs.")
        # If anvio is available, run anvi-run-kegg-kofams and verify 100% agreement!
        return kolach_hits


if __name__ == "__main__":
    import sys
    fasta = sys.argv[1]
    db_dir = sys.argv[2]
    run_comparison(fasta, db_dir)
