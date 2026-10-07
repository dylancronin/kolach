from kolach.metabolism.dbcan import DEFINITIONS, make_cazy_product

CAZY_FILES = ["kolach_cazy_product.tsv", "kolach_cazy_pathways.tsv", "kolach_cazy_metadata.json"]

rule distill_cazy:
    input:
        annotations=f"{OUTDIR}/kolach_annotations.tsv",
        definitions=str(DEFINITIONS),
        genome_map=[SPECIALISATION_INPUTS["gene_genome_map"]] if "gene_genome_map" in SPECIALISATION_INPUTS else []
    output:
        [f"{OUTDIR}/{name}" for name in CAZY_FILES]
    params:
        genome_id=GENOME_ID,
        genome_map=SPECIALISATION_INPUTS.get("gene_genome_map")
    run:
        make_cazy_product(input.annotations, OUTDIR, genome_id=params.genome_id,
                          gene_genome_map=params.genome_map)
