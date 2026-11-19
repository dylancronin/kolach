#!/bin/bash
#SBATCH --job-name=kolach-slurm-controller
#SBATCH --account=PDS0325
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH --output=/fs/scratch/PDS0325/Cronin/kolach/tests/slurm_validation/controller-%j.log
# Controller job: kolach submits kofam, eggnog and deepkoala as separate Slurm jobs.
# Tool settings match tests/anvio_validation/run.sh so per-tool tables are comparable.
set -euo pipefail
cd /fs/scratch/PDS0325/Cronin/kolach/tests/slurm_validation
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate kolach
export PYTHONPATH=/fs/scratch/PDS0325/Cronin/kolach/src
python -c 'import kolach; print("Validated package:", kolach.__file__)'
kolach annotate --protein-fasta ../anvio_validation/proteins.faa \
    --database-dir /fs/scratch/PDS0325/Cronin/kolach-db \
    --databases kofam eggnog deepkoala --output-dir "$PWD/output" --threads 2 \
    --eggnog-sensmode fast --eggnog-dmnd-block-size 4 --eggnog-dmnd-index-chunks 2 \
    --deepkoala-release 202608 --device cpu \
    --executor slurm --slurm-account PDS0325
