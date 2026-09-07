import argparse
from kolach.methods.kofam import annotate

rule annotate_kofam:
    input:
        fasta=PROTEIN_FASTA,
        profiles=f"{KOFAM_DIR}/profiles.hmm",
        ko_list=f"{KOFAM_DIR}/ko_list"
    output:
        tsv=f"{OUTDIR}/kofam_annotations.tsv"
    threads:
        int(config.get("threads", 1))
    run:
        args = argparse.Namespace(
            database_dir=config["database_dir"],
            protein_fasta=input.fasta,
            threads=threads,
            skip_bitscore_heuristic=config.get("skip_bitscore_heuristic", False),
            no_hmmer_prefiltering=config.get("no_hmmer_prefiltering", False),
            heuristic_bitscore_fraction=config.get("heuristic_bitscore_fraction", 0.75),
            heuristic_e_value=config.get("heuristic_e_value", 1e-5),
        )
        with open(output.tsv, "w") as out_f:
            annotate(args, out_f)
