from kolach.pathways import download_kegg_hierarchy

rule download_kegg_hierarchy:
    output:
        KEGG_HIERARCHY
    params:
        data_dir=str(DB_DIR)
    run:
        download_kegg_hierarchy(params.data_dir, force_download=True)
