import argparse
from kolach.methods.deepkoala import annotate

rule annotate_deepkoala:
    input:
        fasta=PROTEIN_FASTA
    output:
        tsv=f"{OUTDIR}/deepkoala_annotations.tsv"
    threads:
        int(config.get("threads", 1))
    run:
        args = argparse.Namespace(
            database_dir=config["database_dir"],
            protein_fasta=input.fasta,
            threads=threads,
            deepkoala_model=config.get("deepkoala_model", "full"),
            deepkoala_release=config.get("deepkoala_release", "latest"),
            device=config.get("device", "auto"),
            batch_size=int(config.get("batch_size", 64)),
            detail=config.get("detail", False),
        )
        annotate(args, output.tsv)
