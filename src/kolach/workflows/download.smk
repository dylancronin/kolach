from pathlib import Path

DB_DIR = Path(config["database_dir"]).expanduser().resolve()
SELECTED_DBS = set(config.get("databases", []))

targets = []

if "dbcan" in SELECTED_DBS:
    DBCAN_DIR = str(DB_DIR / "dbcan")
    targets.append(f"{DBCAN_DIR}/kolach_dbcan_database.json")
    include: "rules/download_dbcan.smk"

if "eggnog" in SELECTED_DBS:
    EGGNOG_DIR = str(DB_DIR / "eggnog")
    targets.append(f"{EGGNOG_DIR}/.download_complete")
    include: "rules/download_eggnog.smk"

if "kofam" in SELECTED_DBS:
    KOFAM_DIR = str(DB_DIR / "kofam")
    targets.append(f"{KOFAM_DIR}/.download_complete")
    include: "rules/download_kofam.smk"

if "deepkoala" in SELECTED_DBS:
    DEEPKOALA_DIR = str(DB_DIR / "deepkoala")
    targets.append(f"{DEEPKOALA_DIR}/.download_complete")
    include: "rules/download_deepkoala.smk"

rule all:
    default_target: True
    input:
        targets
