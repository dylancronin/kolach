import argparse
from kolach.methods.dbcan import ASSETS, annotate as annotate_dbcan

rule annotate_dbcan:
    input:
        fasta=PROTEIN_FASTA,
        assets=[f"{DBCAN_DIR}/{name}" for name in ASSETS]
    output:
        tsv=f"{OUTDIR}/dbcan_annotations.tsv",
        evidence=f"{OUTDIR}/dbcan_evidence.tsv",
        metadata=f"{OUTDIR}/dbcan_metadata.json",
        overview=f"{OUTDIR}/dbcan/overview.tsv",
        hmm=f"{OUTDIR}/dbcan/dbCAN_hmm_results.tsv",
        sub=f"{OUTDIR}/dbcan/dbCANsub_hmm_results.tsv",
        diamond=f"{OUTDIR}/dbcan/diamond.out",
        search_input=f"{OUTDIR}/dbcan/input.faa",
        ids=f"{OUTDIR}/dbcan/gene_ids.tsv",
        log=f"{OUTDIR}/dbcan/run.log"
    threads:
        int(config.get("threads", 1))
    run:
        args = argparse.Namespace(database_dir=str(DB_DIR), protein_fasta=input.fasta, threads=threads)
        annotate_dbcan(args, output.tsv)
