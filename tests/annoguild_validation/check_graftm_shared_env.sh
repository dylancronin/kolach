#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --job-name=kolach-graftm-shared-env
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
# Same placement check as check_graftm_hmmer.sh, using the shared kolach
# environment (pplacer 1.1.alpha23) instead of the validation clone (alpha19).
set -euo pipefail
cd /fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate kolach
python check_graftm_hmmer.py "$SLURM_CPUS_PER_TASK" shared_env
