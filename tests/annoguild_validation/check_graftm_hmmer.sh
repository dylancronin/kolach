#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --job-name=kolach-graftm-hmmer
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
set -euo pipefail
cd /fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate "$PWD/env"
python check_graftm_hmmer.py "$SLURM_CPUS_PER_TASK"
