from kolach.metabolism.specialise import specialise
from kolach.metabolism.evaluate import DATA_DIR

SPECIALISATION_FILES = [
    "kolach_product_refined.tsv", "kolach_specialisations.tsv",
    "kolach_metabolic_pathways.tsv", "kolach_specialisation_assessment.tsv",
    "kolach_specialisation_metadata.json", "kolach_cazy_substrates.tsv",
]
GENE_GENOME_MAP = (str(Path(config["gene_genome_map"]).expanduser().resolve())
                   if config.get("gene_genome_map") else None)
# Rule data are inputs so edits to bundled definitions trigger a rerun.
SPECIALISATION_RULES = sorted(str(path) for path in DATA_DIR.rglob("*") if path.is_file())
GENOME_ID = config.get("genome_id")
# The CLI passes a JSON config file, so string IDs such as '001' stay strings.
if GENOME_ID is not None:
    GENOME_ID = str(GENOME_ID)

# Light step: runs in the controller job under --executor slurm.
localrules: specialise_genomes

rule specialise_genomes:
    input:
        annotations=f"{OUTDIR}/kolach_annotations.tsv",
        markers=f"{OUTDIR}/kolach_graftm_markers.tsv",
        genome_map=[GENE_GENOME_MAP] if GENE_GENOME_MAP else [],
        graftm_manifest=f"{GRAFTM_DIR}/manifest.json",
        rules=SPECIALISATION_RULES
    output:
        [f"{OUTDIR}/{filename}" for filename in SPECIALISATION_FILES]
    params:
        # Config used only inside run is invisible to Snakemake's change tracking.
        genome_id=GENOME_ID,
        genome_map=GENE_GENOME_MAP
    run:
        specialise(
            input.annotations, OUTDIR,
            genome_id=params.genome_id,
            gene_genome_map=params.genome_map,
            markers=input.markers,
            cazy_from_dbcan=True,
            graftm_manifest=input.graftm_manifest,
        )
