import os
from pathlib import Path
import subprocess
from kolach.methods.eggnog import ensure_package_data_dir

rule download_eggnog:
    output:
        f"{EGGNOG_DIR}/.download_complete"
    params:
        data_dir=EGGNOG_DIR
    run:
        # Create the site-packages/data dir at job runtime (not workflow parse time)
        ensure_package_data_dir()
        Path(params.data_dir).mkdir(parents=True, exist_ok=True)
        os.environ["EGGNOG_DATA_DIR"] = params.data_dir
        subprocess.run(
            [
                "download_eggnog_data.py",
                "-y",
                "--release", "3.0",
                "--data_dir", params.data_dir,
            ],
            check=True,
        )
        Path(output[0]).touch()
