import argparse
from kolach.methods.kofam import annotate as annotate_kofam

rule annotate_kofam:
    input:
        fasta=PROTEIN_FASTA,
        profiles=f"{KOFAM_DIR}/profiles.hmm",
        ko_list=f"{KOFAM_DIR}/ko_list"
    output:
        tsv=f"{OUTDIR}/kofam_annotations.tsv"
    threads:
        int(config.get("threads", 1))
    params:
        skip_bitscore_heuristic=str(config.get("skip_bitscore_heuristic", False)).lower() in ("true", "1"),
        no_hmmer_prefiltering=str(config.get("no_hmmer_prefiltering", False)).lower() in ("true", "1"),
        heuristic_bitscore_fraction=float(config.get("heuristic_bitscore_fraction", 0.75)),
        heuristic_e_value=float(config.get("heuristic_e_value", 1e-5))
    run:
        args = argparse.Namespace(
            database_dir=config["database_dir"],
            protein_fasta=input.fasta,
            threads=threads,
            skip_bitscore_heuristic=params.skip_bitscore_heuristic,
            no_hmmer_prefiltering=params.no_hmmer_prefiltering,
            heuristic_bitscore_fraction=params.heuristic_bitscore_fraction,
            heuristic_e_value=params.heuristic_e_value,
        )
        with open(output.tsv, "w") as out_f:
            annotate_kofam(args, out_f)