from kolach.methods.graftm import PACKAGES, download_graftm

rule download_graftm:
    output:
        marker=f"{GRAFTM_DIR}/.download_complete",
        manifest=f"{GRAFTM_DIR}/manifest.json"
    params:
        data_dir=GRAFTM_DIR
    run:
        download_graftm(params.data_dir)
