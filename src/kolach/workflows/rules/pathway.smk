from pathlib import Path
from kolach.pathways import annotate_pathways

rule annotate_pathways:
    input:
        tsv=f"{OUTDIR}/.kolach_annotations_raw.tsv"
    output:
        tsv=f"{OUTDIR}/kolach_annotations.tsv"
    run:
        annotate_pathways(
            annotation_table=input.tsv,
            output_tsv=output.tsv,
            database_dir=config["database_dir"],
            ko_col=config.get("ko_column", "accepted_ko"),
        )
