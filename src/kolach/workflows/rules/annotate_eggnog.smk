import argparse
from pathlib import Path
from kolach.methods.eggnog import annotate as annotate_eggnog

def get_eggnog_db_input(wildcards):
    dl_marker = Path(EGGNOG_DIR) / ".download_complete"
    if dl_marker.exists():
        return str(dl_marker)
    dmnd = Path(EGGNOG_DIR) / "eggnog_proteins.dmnd"
    if dmnd.exists():
        return str(dmnd)
    db = Path(EGGNOG_DIR) / "eggnog.db"
    if db.exists():
        return str(db)
    return str(dl_marker)

rule annotate_eggnog:
    input:
        fasta=PROTEIN_FASTA,
        db=get_eggnog_db_input
    output:
        tsv=f"{OUTDIR}/eggnog_annotations.tsv"
    threads:
        int(config.get("threads", 1))
    params:
        eggnog_mode=config.get("eggnog_mode", "diamond"),
        eggnog_sensmode=config.get("eggnog_sensmode", None),
        eggnog_temp_dir=config.get("eggnog_temp_dir", None),
        eggnog_dmnd_block_size=config.get("eggnog_dmnd_block_size", None),
        eggnog_dmnd_index_chunks=config.get("eggnog_dmnd_index_chunks", None)
    run:
        # Option normalization (None/'' handling and type coercion) happens once in
        # kolach.methods.eggnog.annotate; raw config values are forwarded as-is.
        args = argparse.Namespace(
            database_dir=config["database_dir"],
            protein_fasta=input.fasta,
            threads=threads,
            eggnog_mode=params.eggnog_mode,
            eggnog_sensmode=params.eggnog_sensmode,
            eggnog_temp_dir=params.eggnog_temp_dir,
            eggnog_dmnd_block_size=params.eggnog_dmnd_block_size,
            eggnog_dmnd_index_chunks=params.eggnog_dmnd_index_chunks,
        )
        annotate_eggnog(args, output.tsv)