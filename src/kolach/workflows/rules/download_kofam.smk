from kolach.methods.kofam import download_kofam

rule download_kofam:
    output:
        f"{KOFAM_DIR}/.download_complete"
    params:
        data_dir=KOFAM_DIR
    run:
        download_kofam(params.data_dir)
        Path(output[0]).touch()
