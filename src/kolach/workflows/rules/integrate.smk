from pathlib import Path
from kolach.integrate import integrate_annotations

rule integrate_annotations:
    input:
        fasta=PROTEIN_FASTA,
        tables=targets
    output:
        tsv=f"{OUTDIR}/kolach_annotations.tsv",
        evidence=f"{OUTDIR}/kolach_evidence.tsv"
    run:
        kofam_file = f"{OUTDIR}/kofam_annotations.tsv" if "kofam" in SELECTED_DBS else None
        deepkoala_file = f"{OUTDIR}/deepkoala_annotations.tsv" if "deepkoala" in SELECTED_DBS else None
        eggnog_file = f"{OUTDIR}/eggnog_annotations.tsv" if "eggnog" in SELECTED_DBS else None

        integrate_annotations(
            protein_fasta=input.fasta,
            kofam_tsv=kofam_file,
            deepkoala_tsv=deepkoala_file,
            eggnog_tsv=eggnog_file,
            output_tsv=output.tsv,
            evidence_tsv=output.evidence,
            database_dir=config.get("database_dir"),
            eggnog_min_bitscore=float(config.get("eggnog_min_bitscore", 60.0)),
            eggnog_max_evalue=float(config.get("eggnog_max_evalue", 1e-5)),
            eggnog_filter_multi=config.get("eggnog_filter_multi", "disambiguate"),
            conflict_strategy=config.get("conflict_strategy", "multiple"),
        )
