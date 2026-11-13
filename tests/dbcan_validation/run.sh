#!/bin/bash
#SBATCH --job-name=kolach_dbcan_validation
#SBATCH --account=PDS0325
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=tests/dbcan_validation/slurm-%j.out
#SBATCH --error=tests/dbcan_validation/slurm-%j.err

set -euo pipefail
TASK_ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}
VALIDATION_DIR="$TASK_ROOT/tests/dbcan_validation"
export PATH="$VALIDATION_DIR/env/bin:/fs/ess/PDS0325/bioinformatic_tools/miniforge3/envs/kolach/bin:$PATH"
export PYTHONPATH="$TASK_ROOT/src"
export XDG_CACHE_HOME="$VALIDATION_DIR/cache"
export NXF_HOME="$VALIDATION_DIR/.nextflow_home"
export NXF_OFFLINE=true
export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS='-Xms256m -Xmx2g'
export OMP_NUM_THREADS=4
cd "$VALIDATION_DIR"
test -f db/dbcan/.download_complete
run_dbcan version
diamond version
./nextflow -version

# Real Kolach workflow; searches only dbCAN and retains the original gene IDs.
/usr/bin/time -v python -m kolach.cli annotate \
    --databases dbcan --protein-fasta ../anvio_validation/proteins.faa \
    --database-dir db --output-dir kolach_output --threads 4 --dbcan-keep-raw \
    > kolach.log 2>&1

# Actual, unmodified DRAM2 Nextflow process, with the same host tools and DBs.
/usr/bin/time -v ./nextflow run reference.nf -c reference.config \
    --fasta "$VALIDATION_DIR/kolach_output/dbcan/input.faa" \
    --db "$VALIDATION_DIR/db/dbcan" -work-dir work -resume \
    -with-trace trace.tsv > reference.log 2>&1

python compare.py
