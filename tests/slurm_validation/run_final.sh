#!/bin/bash
#SBATCH --job-name=kolach-slurm-final
#SBATCH --account=PDS0325
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/fs/scratch/PDS0325/Cronin/kolach/tests/slurm_validation/final-%j.log
# Final-code check after review fixes: one DeepKOALA job, custom profile with default-resources.
set -euo pipefail
cd /fs/scratch/PDS0325/Cronin/kolach/tests/slurm_validation
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate kolach
export PYTHONPATH=/fs/scratch/PDS0325/Cronin/kolach/src
kolach annotate --protein-fasta ../anvio_validation/proteins.faa \
    --database-dir /fs/scratch/PDS0325/Cronin/kolach-db --databases deepkoala \
    --output-dir "$PWD/final_output" --threads 1 --deepkoala-release 202608 --device cpu \
    --executor slurm --slurm-account PDS0325 --slurm-profile "$PWD/final_profile" --jobs 2
