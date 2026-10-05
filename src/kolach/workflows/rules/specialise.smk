from kolach.metabolism.specialise import specialise
from kolach.metabolism.evaluate import DATA_DIR

# Register supplemental files as inputs so Snakemake reruns when they change.
SPECIALISATION_FILES = [
    "kolach_product_refined.tsv", "kolach_specialisations.tsv",
    "kolach_metabolic_pathways.tsv", "kolach_specialisation_assessment.tsv",
    "kolach_specialisation_metadata.json",
]
SPECIALISATION_INPUTS = [str(Path(config[name]).expanduser().resolve())
                         for name in ("gene_genome_map", "markers", "cazy_product")
                         if config.get(name)]
SPECIALISATION_RULES = list(DATA_DIR.glob("*.tsv")) + [DATA_DIR / "manifest.json"]

rule specialise_genomes:
    input:
        annotations=f"{OUTDIR}/kolach_annotations.tsv",
        supplemental=SPECIALISATION_INPUTS,
        rules=[str(path) for path in SPECIALISATION_RULES]
    output:
        [f"{OUTDIR}/{filename}" for filename in SPECIALISATION_FILES]
    run:
        specialise(
            input.annotations, OUTDIR,
            genome_id=config.get("genome_id"),
            gene_genome_map=config.get("gene_genome_map"),
            markers=config.get("markers"),
            cazy_product=config.get("cazy_product"),
        )
