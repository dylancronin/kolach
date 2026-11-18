#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --job-name=kolach-syncomm-rerun
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation/rerun-%j.log
# Re-run only GraftM + specialisation on the finished SynComm annotations with
# the reviewed code (annotation searches are reused; see dry run in RESULTS.md).
set -euo pipefail
V=/fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate "$V/env"
export PYTHONPATH="$V/code/src"
cd "$V/syncomm"
echo "annotation searches: 41affce; graftm/specialisation: $(git -C "$V/code" rev-parse HEAD)" > kolach_commit.txt
# Same settings as run_syncomm.sh, written as the CLI would.
python - <<'PY'
import json
from kolach.cli import annotate_config, build_parser
V = "/fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation"
args = build_parser().parse_args([
    "annotate", "--protein-fasta", f"{V}/syncomm/proteins.faa", "--database-dir", f"{V}/db",
    "--databases", "kofam", "eggnog", "deepkoala", "dbcan", "--output-dir", f"{V}/syncomm/kolach_output",
    "--threads", "20", "--eggnog-sensmode", "fast", "--eggnog-dmnd-block-size", "4",
    "--eggnog-dmnd-index-chunks", "2", "--deepkoala-release", "202608", "--device", "cpu",
    "--add-specialisations", "--gene-genome-map", f"{V}/syncomm/gene_genome_map.tsv"])
json.dump(annotate_config(args), open("config.json", "w"), indent=1)
PY
snakemake -n -s "$V/code/src/kolach/workflows/annotate.smk" --cores "$SLURM_CPUS_PER_TASK" --configfile config.json | grep -E "^rule "
snakemake -s "$V/code/src/kolach/workflows/annotate.smk" --cores "$SLURM_CPUS_PER_TASK" --configfile config.json
python "$V/compare_syncomm.py"
