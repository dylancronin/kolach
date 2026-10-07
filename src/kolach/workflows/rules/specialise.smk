from kolach.metabolism.specialise import specialise
from kolach.metabolism.evaluate import DATA_DIR

# Register supplemental files as inputs so Snakemake reruns when they change.
SPECIALISATION_FILES = [
    "kolach_product_refined.tsv", "kolach_specialisations.tsv",
    "kolach_metabolic_pathways.tsv", "kolach_specialisation_assessment.tsv",
    "kolach_specialisation_metadata.json",
]
if "dbcan" in SELECTED_DBS:
    if config.get("cazy_product"):
        raise ValueError("dbCAN annotation and an external cazy_product cannot both supply CAZy calls")
    SPECIALISATION_FILES += ["kolach_cazy_product.tsv", "kolach_cazy_pathways.tsv", "kolach_cazy_metadata.json"]
SPECIALISATION_INPUTS = {name: str(Path(config[name]).expanduser().resolve())
                         for name in ("gene_genome_map", "markers", "cazy_product")
                         if config.get(name)}
SPECIALISATION_RULES = list(DATA_DIR.glob("*.tsv")) + [DATA_DIR / "manifest.json"]
GENOME_ID = config.get("genome_id")
# The CLI wraps string IDs to survive Snakemake's child-job config coercion.
# Direct workflow users can still provide an ordinary scalar genome_id.
if isinstance(GENOME_ID, list):
    GENOME_ID = GENOME_ID[0]
if GENOME_ID is not None:
    GENOME_ID = str(GENOME_ID)

rule specialise_genomes:
    input:
        annotations=f"{OUTDIR}/kolach_annotations.tsv",
        supplemental=list(SPECIALISATION_INPUTS.values()),
        rules=[str(path) for path in SPECIALISATION_RULES]
    output:
        [f"{OUTDIR}/{filename}" for filename in SPECIALISATION_FILES]
    params:
        # Config used only inside run is invisible to Snakemake's change tracking.
        # Track grouping and input roles as well as the input files themselves.
        genome_id=GENOME_ID,
        genome_map=SPECIALISATION_INPUTS.get("gene_genome_map"),
        markers=SPECIALISATION_INPUTS.get("markers"),
        cazy_product=SPECIALISATION_INPUTS.get("cazy_product")
    run:
        cazy_product = params.cazy_product
        if "dbcan" in SELECTED_DBS:
            from kolach.metabolism.dbcan import make_cazy_product
            cazy_product = make_cazy_product(input.annotations, OUTDIR, genome_id=params.genome_id,
                                             gene_genome_map=params.genome_map)
        specialise(
            input.annotations, OUTDIR,
            genome_id=params.genome_id,
            gene_genome_map=params.genome_map,
            markers=params.markers,
            cazy_product=cazy_product,
        )
