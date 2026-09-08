from pathlib import Path

DB_DIR = Path(config["database_dir"]).expanduser().resolve()
PROTEIN_FASTA = str(Path(config["protein_fasta"]).expanduser().resolve())
OUTDIR = str(Path(config["output_dir"]).expanduser().resolve())
SELECTED_DBS = set(config.get("databases", ["kofam"]))

targets = []

if "kofam" in SELECTED_DBS:
    KOFAM_DIR = str(DB_DIR / "kofam")
    targets.append(f"{OUTDIR}/kofam_annotations.tsv")
    include: "rules/annotate_kofam.smk"

if "deepkoala" in SELECTED_DBS:
    DEEPKOALA_DIR = str(DB_DIR / "deepkoala")
    targets.append(f"{OUTDIR}/deepkoala_annotations.tsv")
    include: "rules/annotate_deepkoala.smk"

rule all:
    input:
        targets
