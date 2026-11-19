#!/bin/bash
#SBATCH --job-name=kolach-slurm-cancel
#SBATCH --account=PDS0325
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/fs/scratch/PDS0325/Cronin/kolach/tests/slurm_validation/cancel_check/controller-%j.log
# Cancel check: scancel this controller while its kofam job runs.
set -euo pipefail
cd /fs/scratch/PDS0325/Cronin/kolach/tests/slurm_validation/cancel_check
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate kolach
export PYTHONPATH=/fs/scratch/PDS0325/Cronin/kolach/src
exec kolach annotate --protein-fasta ../../anvio_validation/proteins.faa \
    --database-dir /fs/scratch/PDS0325/Cronin/kolach-db --databases kofam \
    --output-dir "$PWD/final_output" --threads 1 \
    --executor slurm --slurm-account PDS0325 --slurm-profile "$PWD/../final_profile" --jobs 2
