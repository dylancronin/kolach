"""Fetch the unmodified DRAM2 dbCAN module and its table combiner."""
from pathlib import Path
from urllib.request import urlopen
import hashlib
import json

ROOT = Path(__file__).resolve().parent
COMMIT = "4536cd93c97c07d783fab9e9ef99e56b271f3782"
FILES = {
    "module/main.nf": "modules/nf-core/rundbcan/cazymeannotation/main.nf",
    "module/environment.yml": "modules/nf-core/rundbcan/cazymeannotation/environment.yml",
    "combine_annotations.py": "bin/combine_annotations.py",
    "utils/logger.py": "bin/utils/logger.py",
}


def main():
    manifest = {"commit": COMMIT, "files": {}}
    for target, source in FILES.items():
        url = f"https://raw.githubusercontent.com/BortonWrightonLabs/DRAM/{COMMIT}/{source}"
        path = ROOT / "upstream" / target
        path.parent.mkdir(parents=True, exist_ok=True)
        with urlopen(url, timeout=120) as response:
            data = response.read()
        path.write_bytes(data)
        manifest["files"][source] = {"url": url, "sha256": hashlib.sha256(data).hexdigest()}
    (ROOT / "reference_sources.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
