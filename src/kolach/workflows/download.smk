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

if "deepkoala" in SELECTED_DBS:
    from kolach.methods.deepkoala import resolve_release

    DEEPKOALA_DIR = str(DB_DIR / "deepkoala")
    # Resolve 'latest' now so a newly published release becomes a new target
    # instead of being masked by an existing download marker.
    DEEPKOALA_RELEASE = resolve_release(config.get("deepkoala_release", "latest"))
    targets.append(f"{DEEPKOALA_DIR}/{DEEPKOALA_RELEASE}/manifest.json")
    include: "rules/download_deepkoala.smk"

if targets:
    # KEGG Orthology hierarchy used by --add-pathways, fetched up front so the
    # annotation workflow does not need network access.
    KEGG_HIERARCHY = str(DB_DIR / "ko00001.json")
    targets.append(KEGG_HIERARCHY)
    include: "rules/download_kegg.smk"

rule all:
    input:
        targets
