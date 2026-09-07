import sys
from pathlib import Path

# Create default site-packages data dir if missing (fixes eggnog-mapper v3 package bug)
(Path(sys.prefix) / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages" / "data").mkdir(parents=True, exist_ok=True)

rule download_eggnog:
    output:
        f"{EGGNOG_DIR}/.download_complete"
    params:
        data_dir=EGGNOG_DIR
    shell:
        """
        mkdir -p {params.data_dir:q}
        export EGGNOG_DATA_DIR={params.data_dir:q}
        download_eggnog_data.py -y --release 3.0 --data_dir {params.data_dir:q} 
        touch {output:q}
        """
