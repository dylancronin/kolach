from pathlib import Path

DB_DIR = Path(config["database_dir"]).expanduser().resolve()
PROTEIN_FASTA = str(Path(config["protein_fasta"]).expanduser().resolve())
OUTDIR = str(Path(config["output_dir"]).expanduser().resolve())
SELECTED_DBS = set(config.get("databases", ["kofam"]))
ADD_PATHWAYS = str(config.get("add_pathways", "false")).lower() in ("true", "1", "yes")

targets = []
supplementary_targets = []

if "dbcan" in SELECTED_DBS:
    DBCAN_DIR = str(DB_DIR / "dbcan")
    targets.append(f"{OUTDIR}/dbcan_annotations.tsv")
    include: "rules/annotate_dbcan.smk"
    supplementary_targets.extend(path for path in rules.annotate_dbcan.output if path not in targets)

if "kofam" in SELECTED_DBS:
    KOFAM_DIR = str(DB_DIR / "kofam")
    targets.append(f"{OUTDIR}/kofam_annotations.tsv")
    include: "rules/annotate_kofam.smk"

if "deepkoala" in SELECTED_DBS:
    DEEPKOALA_DIR = str(DB_DIR / "deepkoala")
    targets.append(f"{OUTDIR}/deepkoala_annotations.tsv")
    include: "rules/annotate_deepkoala.smk"

if "eggnog" in SELECTED_DBS:
    EGGNOG_DIR = str(DB_DIR / "eggnog")
    targets.append(f"{OUTDIR}/eggnog_annotations.tsv")
    include: "rules/annotate_eggnog.smk"

if targets:
    include: "rules/integrate.smk"
    if ADD_PATHWAYS:
        include: "rules/pathway.smk"
    all_targets = targets + supplementary_targets + [f"{OUTDIR}/kolach_annotations.tsv", f"{OUTDIR}/kolach_evidence.tsv"]
else:
    all_targets = []

rule all:
    default_target: True
    input:
        all_targets
