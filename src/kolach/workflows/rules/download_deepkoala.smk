from pathlib import Path
from kolach.methods.deepkoala import download_deepkoala

rule download_deepkoala:
    output:
        f"{DEEPKOALA_DIR}/{DEEPKOALA_RELEASE}/manifest.json"
    params:
        data_dir=DEEPKOALA_DIR,
        release=DEEPKOALA_RELEASE
    run:
        download_deepkoala(params.data_dir, release=params.release)
        Path(params.data_dir, ".download_complete").touch()
