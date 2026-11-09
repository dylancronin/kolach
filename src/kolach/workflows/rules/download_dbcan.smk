from kolach.methods.dbcan import ASSETS, download_dbcan

rule download_dbcan:
    output:
        marker=f"{DBCAN_DIR}/.download_complete",
        manifest=f"{DBCAN_DIR}/manifest.json",
        assets=[f"{DBCAN_DIR}/{name}" for name in ASSETS]
    params:
        data_dir=DBCAN_DIR
    run:
        download_dbcan(params.data_dir)
