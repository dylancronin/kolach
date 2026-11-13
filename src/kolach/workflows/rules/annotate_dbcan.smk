import argparse
from kolach.methods.dbcan import ASSETS, annotate as annotate_dbcan

DBCAN_KEEP_RAW = str(config.get("dbcan_keep_raw", "false")).lower() in ("true", "1", "yes")
DBCAN_RAW_OUTPUTS = ([f"{OUTDIR}/dbcan_evidence.tsv", f"{OUTDIR}/dbcan_metadata.json"] + [
    f"{OUTDIR}/dbcan/{name}" for name in (
        "overview.tsv", "dbCAN_hmm_results.tsv", "dbCANsub_hmm_results.tsv",
        "diamond.out", "input.faa", "gene_ids.tsv", "run.log",
    )
]) if DBCAN_KEEP_RAW else []

rule annotate_dbcan:
    input:
        fasta=PROTEIN_FASTA,
        assets=[f"{DBCAN_DIR}/{name}" for name in ASSETS]
    output:
        tsv=f"{OUTDIR}/dbcan_overview.tsv",
        evidence=temp(f"{OUTDIR}/.dbcan_evidence.tsv"),
        raw=DBCAN_RAW_OUTPUTS
    threads:
        int(config.get("threads", 1))
    run:
        args = argparse.Namespace(database_dir=str(DB_DIR), protein_fasta=input.fasta, threads=threads,
                                  dbcan_keep_raw=DBCAN_KEEP_RAW)
        annotate_dbcan(args, output.tsv)
