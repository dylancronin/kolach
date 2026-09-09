"""eggNOG search and annotation wrapper for kolach.

This module provides annotation using eggNOG-mapper (emapper.py), mapping
proteins to eggNOG orthologous groups and KEGG Orthologies (KOs).
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

REFERENCE = "https://github.com/eggnogdb/eggnog-mapper"


def annotate(args, output_file):
    """Run upstream eggNOG-mapper (emapper.py) on protein sequences.

    Args:
        args: Namespace containing:
            - database_dir: Base directory containing databases
            - protein_fasta: Input protein FASTA file
            - threads: Number of CPUs / threads
            - eggnog_mode: Search mode (diamond, mmseqs)
            - eggnog_sensmode: Sensitivity mode for diamond (or None for eggNOG default)
            - eggnog_temp_dir: Base directory for temporary files (or None for system default)
            - eggnog_dmnd_block_size: DIAMOND block size in billions of sequence letters (or None)
            - eggnog_dmnd_index_chunks: Number of chunks for processing DIAMOND seed index (or None)
        output_file: Path to destination TSV output file.
    """
    if not shutil.which("emapper.py"):
        raise FileNotFoundError(
            "The 'emapper.py' executable was not found in PATH. "
            "Please install eggnog-mapper (e.g. conda install -c bioconda eggnog-mapper) "
            "or ensure your conda environment is activated."
        )

    db_dir = Path(args.database_dir).expanduser().resolve() / "eggnog"
    if not db_dir.is_dir():
        raise FileNotFoundError(
            f"No eggNOG database directory found at {db_dir}. "
            "Please run 'kolach download --databases eggnog' first."
        )

    if not any(db_dir.iterdir()):
        raise FileNotFoundError(
            f"The eggNOG database directory at {db_dir} is empty. "
            "Please run 'kolach download --databases eggnog' first."
        )

    output_path = Path(output_file).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    base_tmp_dir = getattr(args, "eggnog_temp_dir", None)
    if str(base_tmp_dir).lower() in ("none", ""):
        base_tmp_dir = None
    if base_tmp_dir is not None:
        base_tmp_dir = str(Path(base_tmp_dir).expanduser().resolve())
        Path(base_tmp_dir).mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="kolach_eggnog_", dir=base_tmp_dir) as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        output_prefix = "eggnog_res"

        threads = getattr(args, "threads", 1)
        mode = getattr(args, "eggnog_mode", "diamond")
        sensmode = getattr(args, "eggnog_sensmode", None)
        if str(sensmode).lower() in ("none", ""):
            sensmode = None

        block_size = getattr(args, "eggnog_dmnd_block_size", None)
        if str(block_size).lower() in ("none", ""):
            block_size = None
        if block_size is not None:
            block_size = float(block_size)
            if block_size <= 0:
                raise ValueError(f"eggnog_dmnd_block_size must be positive, got {block_size}")

        index_chunks = getattr(args, "eggnog_dmnd_index_chunks", None)
        if str(index_chunks).lower() in ("none", ""):
            index_chunks = None
        if index_chunks is not None:
            index_chunks = int(index_chunks)
            if index_chunks <= 0:
                raise ValueError(f"eggnog_dmnd_index_chunks must be a positive integer, got {index_chunks}")

        cmd = [
            "emapper.py",
            "-i", str(Path(args.protein_fasta).expanduser().resolve()),
            "--itype", "proteins",
            "-m", mode,
            "--data_dir", str(db_dir),
            "--cpu", str(threads),
            "-o", output_prefix,
            "--output_dir", str(tmp_dir_path),
            "--temp_dir", str(tmp_dir_path),
            "--override",
        ]

        if mode == "diamond":
            if sensmode is not None:
                cmd.extend(["--dmnd_sensmode", str(sensmode)])
            if block_size is not None:
                cmd.extend(["--dmnd_block_size", str(block_size)])
            if index_chunks is not None:
                cmd.extend(["--dmnd_index_chunks", str(index_chunks)])

        env = os.environ.copy()
        env["EGGNOG_DATA_DIR"] = str(db_dir)

        result = subprocess.run(cmd, env=env)
        if result.returncode != 0:
            raise RuntimeError(f"emapper.py failed with return code {result.returncode}")

        annotations_file = tmp_dir_path / f"{output_prefix}.emapper.annotations"
        if not annotations_file.exists():
            raise FileNotFoundError(
                f"Expected emapper output file {annotations_file} was not generated."
            )

        with open(annotations_file, "r", encoding="utf-8") as in_f, open(
            output_path, "w", encoding="utf-8", newline=""
        ) as out_f:
            for line in in_f:
                if line.startswith("##"):
                    continue
                if line.startswith("#query"):
                    line = line[1:]
                out_f.write(line)
