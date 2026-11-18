#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --job-name=kolach-syncomm-specialisation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --output=/fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation/syncomm-%j.log
# Full Kolach run on the six SynComm genomes AnnoGuild was run on, using the
# pinned checkout in code/ and the validation environment in env/.
set -euo pipefail
V=/fs/scratch/PDS0325/Cronin/kolach/tests/annoguild_validation
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate "$V/env"
export PYTHONPATH="$V/code/src"
cd "$V/syncomm"
git -C "$V/code" rev-parse HEAD > kolach_commit.txt
python -m kolach.cli annotate \
  --protein-fasta "$V/syncomm/proteins.faa" --database-dir "$V/db" \
  --databases kofam eggnog deepkoala dbcan --output-dir "$V/syncomm/kolach_output" \
  --threads "$SLURM_CPUS_PER_TASK" --eggnog-sensmode fast --eggnog-dmnd-block-size 4 \
  --eggnog-dmnd-index-chunks 2 --deepkoala-release 202608 --device cpu \
  --add-specialisations --gene-genome-map "$V/syncomm/gene_genome_map.tsv"
python "$V/compare_syncomm.py"
