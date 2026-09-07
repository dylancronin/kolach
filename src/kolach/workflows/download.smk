from pathlib import Path

DB_DIR = Path(config["database_dir"]).expanduser().resolve()
SELECTED_DBS = set(config.get("databases", []))

targets = []

if "eggnog" in SELECTED_DBS:
    EGGNOG_DIR = str(DB_DIR / "eggnog")
    targets.append(f"{EGGNOG_DIR}/.download_complete")
    # Import rules from rules/download_eggnog.smk relative to this file.
    # Inherits global workflow scope (EGGNOG_DIR, config).
    include: "rules/download_eggnog.smk"

rule all:
    input:
        targets
