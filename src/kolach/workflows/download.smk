from pathlib import Path

DB_DIR = Path(config["database_dir"]).expanduser().resolve()
SELECTED_DBS = set(config.get("databases", []))

targets = []

if "eggnog" in SELECTED_DBS:
    EGGNOG_DIR = str(DB_DIR / "eggnog")
    targets.append(f"{EGGNOG_DIR}/.download_complete")
    include: "rules/download_eggnog.smk"

if "kofam" in SELECTED_DBS:
    KOFAM_DIR = str(DB_DIR / "kofam")
    targets.append(f"{KOFAM_DIR}/.download_complete")
    include: "rules/download_kofam.smk"

rule all:
    input:
        targets
