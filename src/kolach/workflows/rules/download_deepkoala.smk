# DEEPKOALA_RELEASE is already resolved against GenomeNet in download.smk.
from kolach.methods.deepkoala import download_release

rule download_deepkoala:
    output:
        f"{DEEPKOALA_DIR}/{DEEPKOALA_RELEASE}/manifest.json"
    params:
        data_dir=DEEPKOALA_DIR,
        release=DEEPKOALA_RELEASE
    run:
        download_release(params.data_dir, params.release)
