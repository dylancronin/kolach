"""Copy upstream rule assets verbatim and record their commits and SHA-256s.

Usage: python scripts/import_dram2_rules.py DRAM_CHECKOUT DRAM_VIZ_CHECKOUT
The DRAM checkout must have its pinned bin/rule_parser submodule initialized.
Run verify_dram2_parity.py after every import, before committing changed assets.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def commit(path):
    return subprocess.check_output(['git', '-C', str(path), 'rev-parse', 'HEAD'], text=True).strip()


def main():
    dram, viz = map(Path, sys.argv[1:])
    parser = dram / 'bin/rule_parser'
    destination = Path(__file__).resolve().parents[1] / 'src/kolach/metabolism/data/dram2'
    destination.mkdir(parents=True, exist_ok=True)
    sources = {
        'traits.tsv': (dram, 'bin/assets/traits_rules.tsv'),
        'common.tsv': (parser, 'src/rules_common.tsv'),
        'LICENSE.GPL-3.0': (parser, 'LICENSE'),
        'LICENSE.dram-viz': (viz, 'LICENSE'),
        **{p.name: (viz, 'dram_viz/data/' + p.name) for p in (viz / 'dram_viz/data').glob('*.tsv')},
    }
    manifest = {}
    for name, (checkout, path) in sources.items():
        data = subprocess.check_output(['git', '-C', str(checkout), 'show', 'HEAD:' + path])
        (destination / name).write_bytes(data)
        manifest[name] = {'commit': commit(checkout), 'path': path,
                          'sha256': hashlib.sha256(data).hexdigest()}
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
