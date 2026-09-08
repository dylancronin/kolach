from pathlib import Path
from kolach.methods.deepkoala import download_deepkoala

rule download_deepkoala:
    output:
        f"{DEEPKOALA_DIR}/.download_complete"
    params:
        data_dir=DEEPKOALA_DIR,
        release=config.get("deepkoala_release", "latest")
    run:
        download_deepkoala(params.data_dir, release=params.release)
        Path(output[0]).touch()
