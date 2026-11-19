#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --job-name=kolach-shared-env-e2e
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
# End-to-end check of the shared kolach environment and shared databases:
# real KOfam + dbCAN + GraftM + specialisation on a 60-protein Afacilis subset.
set -euo pipefail
cd /fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation/shared_env_e2e
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate kolach
rm -rf out
kolach annotate --protein-fasta Afacilis_subset.faa --database-dir /fs/scratch/PDS0325/Cronin/kolach-db \
  --databases kofam --output-dir out --threads "$SLURM_CPUS_PER_TASK" --add-specialisations
