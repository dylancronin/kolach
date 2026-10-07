from kolach.metabolism.dbcan import run_dbcan, DATABASE_FILES

rule annotate_dbcan:
    input:
        fasta=PROTEIN_FASTA,
        references=[f"{DBCAN_DIR}/{name}" for name in DATABASE_FILES]
    output:
        overview=f"{OUTDIR}/dbcan/overview.tsv",
        metadata=f"{OUTDIR}/dbcan/kolach_dbcan_run.json"
    threads:
        int(config.get("threads", 1))
    run:
        run_dbcan(input.fasta, DBCAN_DIR, f"{OUTDIR}/dbcan", threads)
