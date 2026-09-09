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
        sensmode = config.get("eggnog_sensmode", None)
        if sensmode in ("None", "", None):
            sensmode = None

        temp_dir = config.get("eggnog_temp_dir", None)
        if temp_dir in ("None", "", None):
            temp_dir = None

        block_size = config.get("eggnog_dmnd_block_size", None)
        if block_size in ("None", "", None):
            block_size = None
        else:
            block_size = float(block_size)

        index_chunks = config.get("eggnog_dmnd_index_chunks", None)
        if index_chunks in ("None", "", None):
            index_chunks = None
        else:
            index_chunks = int(index_chunks)

        args = argparse.Namespace(
            database_dir=config["database_dir"],
            protein_fasta=input.fasta,
            threads=threads,
            eggnog_mode=config.get("eggnog_mode", "diamond"),
            eggnog_sensmode=sensmode,
            eggnog_temp_dir=temp_dir,
            eggnog_dmnd_block_size=block_size,
            eggnog_dmnd_index_chunks=index_chunks,
        )
        annotate(args, output.tsv)
