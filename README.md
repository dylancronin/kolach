# kolach
Genome annotation targeting KEGG Orthologies (KOs).

`kolach` provides unified workflows to download databases and assign KO identifiers to protein sequences using multiple methods:
- **KOfam**: Profile HMM searches with bitscore threshold and rescue heuristics.
- **DeepKOALA**: Ultra-fast deep learning GRU models for full-length or metagenomic fragment sequences.
- **eggNOG**: Orthology assignment based on eggNOG orthologous groups.

---

## Installation

```bash
conda env create -f environment.yaml
conda activate kolach
```

---

## Usage

### 1. Download Databases

Download required databases (e.g. `kofam`, `deepkoala`, or `eggnog`) into a specified directory:

```bash
kolach download --database-dir /path/to/databases --databases deepkoala kofam
```

Optional DeepKOALA download flags:
- `--deepkoala-release`: Specify release date tag (e.g. `202608` or `latest`, default: `latest`).

### 2. Annotate Proteins

Annotate a protein FASTA file with chosen methods:

```bash
kolach annotate \
    --protein-fasta proteins.faa \
    --database-dir /path/to/databases \
    --databases deepkoala kofam eggnog \
    --output-dir kolach_output \
    --threads 8
```

DeepKOALA-specific options:
- `--deepkoala-model`: Model variant to use (`full` for complete sequences, `frag` for metagenomic fragments; default: `full`).
- `--deepkoala-release`: Model release tag to use (default: `latest`).
- `--device`: Compute device (`auto`, `cpu`, or `cuda`; default: `auto`).
- `--batch-size`: Batch size for inference (default: `64`).
- `--detail`: Emit detailed probabilities, thresholds, and boundary marks.

eggNOG-specific options:
- `--eggnog-mode`: Search mode (`diamond` or `mmseqs`, default: `diamond`).
- `--eggnog-sensmode`: Diamond sensitivity mode (`default`, `fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive`). When left unspecified, eggNOG-mapper's default sensitivity (sensitive, iterative search) is used. Selecting `default` explicitly requests DIAMOND's default (faster, non-sensitive) sensitivity mode. Note: faster sensitivity settings reduce runtime but can lower distant-homolog and KO recovery rates.
- `--eggnog-temp-dir`: Path to base directory for eggNOG temporary files and DIAMOND scratch files (`emappertmp_dmdn_*`). Defaults to the system temporary directory.
- `--eggnog-dmnd-block-size`: DIAMOND block size in billions of sequence letters (`--dmnd_block_size`). When unspecified, DIAMOND automatically tunes this based on host RAM. In HPC/Slurm jobs with constrained cgroup memory limits, automatic host-RAM tuning may exceed job memory allocations, leading to OOM termination; setting this parameter controls memory usage.
- `--eggnog-dmnd-index-chunks`: Number of chunks for processing the DIAMOND seed index (`--dmnd_index_chunks`). When unspecified, upstream automatic tuning is used.

---

### High-Performance eggNOG Execution on HPC & Slurm

eggNOG searches using DIAMOND are compute- and I/O-intensive. For optimal throughput on compute clusters, follow these recommendations:

#### 1. Allocate Multi-Core CPUs and Local Scratch Space

Use `--eggnog-temp-dir` to direct DIAMOND temporary scratch files to fast node-local storage (such as `$TMPDIR` or NVMe scratch) rather than a shared network filesystem (NFS/Lustre):

```bash
kolach annotate \
    --protein-fasta proteins.faa \
    --database-dir /path/to/databases \
    --databases eggnog \
    --output-dir kolach_output \
    --threads 16 \
    --eggnog-temp-dir "$TMPDIR"
```

#### 2. Staging Databases on Node-Local Storage

If node-local disk is available, staging databases locally eliminates network filesystem latency. eggNOG requires databases to reside inside an `eggnog/` subdirectory within the database directory. `kolach` does not automatically duplicate multi-gigabyte database files on every run. You can stage databases manually before executing `kolach`:

```bash
# 1. Create staged directory structure preserving the eggnog/ subdirectory layout
mkdir -p "$TMPDIR/databases/eggnog"

# 2. Stage database files to node-local scratch
cp /shared/databases/eggnog/* "$TMPDIR/databases/eggnog/"

# 3. Point --database-dir to the staging root
kolach annotate \
    --protein-fasta proteins.faa \
    --database-dir "$TMPDIR/databases" \
    --databases eggnog \
    --output-dir kolach_output \
    --threads 16 \
    --eggnog-temp-dir "$TMPDIR"
```

#### 3. Tuning Memory within Slurm Job Allocations

DIAMOND automatically tunes memory according to total physical host RAM. On shared compute nodes, a Slurm job with a restricted memory allocation (e.g. `--mem=32G`) may be killed by Slurm's OOM killer if DIAMOND scales up to total node RAM. Explicitly set `--eggnog-dmnd-block-size` and `--eggnog-dmnd-index-chunks` to constrain memory:

```bash
kolach annotate \
    --protein-fasta proteins.faa \
    --database-dir "$TMPDIR/databases" \
    --databases eggnog \
    --output-dir kolach_output \
    --threads 16 \
    --eggnog-temp-dir "$TMPDIR" \
    --eggnog-dmnd-block-size 2.0 \
    --eggnog-dmnd-index-chunks 4
```
