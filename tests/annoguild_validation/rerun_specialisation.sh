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
git -C "$V/code" rev-parse HEAD > kolach_commit.txt
snakemake -s "$V/code/src/kolach/workflows/annotate.smk" --cores "$SLURM_CPUS_PER_TASK" --configfile config.json
python "$V/compare_syncomm.py"
