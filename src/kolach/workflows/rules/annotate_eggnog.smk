import argparse
from kolach.methods.eggnog import annotate

rule annotate_eggnog:
    input:
        fasta=PROTEIN_FASTA
    output:
        tsv=f"{OUTDIR}/eggnog_annotations.tsv"
    threads:
        int(config.get("threads", 1))
    run:
        args = argparse.Namespace(
            database_dir=config["database_dir"],
            protein_fasta=input.fasta,
            threads=threads,
            eggnog_mode=config.get("eggnog_mode", "diamond"),
            eggnog_sensmode=config.get("eggnog_sensmode", "default"),
            eggnog_dbmem=config.get("eggnog_dbmem", False),
        )
        annotate(args, output.tsv)
