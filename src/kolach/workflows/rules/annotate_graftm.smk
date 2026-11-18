from kolach.methods.graftm import PLACEMENTS_FILE, annotate as annotate_graftm

# Runs after integration: proteins are chosen by their accepted or alternative KO.
rule annotate_graftm:
    input:
        fasta=PROTEIN_FASTA,
        annotations=f"{OUTDIR}/kolach_annotations.tsv",
        packages=f"{GRAFTM_DIR}/.download_complete",
        placements=str(PLACEMENTS_FILE)
    output:
        markers=f"{OUTDIR}/kolach_graftm_markers.tsv",
        log=f"{OUTDIR}/kolach_graftm_markers.tsv.log"
    params:
        package_dir=GRAFTM_DIR
    threads:
        int(config.get("threads", 1))
    run:
        annotate_graftm(input.fasta, input.annotations, params.package_dir, output.markers,
                        threads=threads, workdir=f"{OUTDIR}/.graftm_work")
