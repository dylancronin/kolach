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
