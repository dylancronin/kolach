"""Run released dbCAN's real overview generator and AnnoGuild CAZy gate.

Usage: PYTHONPATH=src python scripts/verify_dbcan_parity.py DBCAN_CHECKOUT APPEND_CAZY_PY
Requires pandas, numpy and biopython. No HMM/DIAMOND search is simulated here:
these checks start at synthetic, post-filter search hits and test their import
and metabolic interpretation. Full database searches require a Linux runtime.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

import pandas as pd
from kolach.metabolism.dbcan import DEFINITIONS, cazy_calls, read_overview


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('dbcan', type=Path)
    parser.add_argument('append_cazy', type=Path)
    parser.add_argument('--write-fixture', action='store_true')
    args = parser.parse_args()
    sha = subprocess.check_output(['git', '-C', str(args.dbcan), 'rev-parse', 'HEAD'], text=True).strip()
    assert sha == '614c93f896939042ae5bd574b9c6b971e80803f6', 'Use released dbCAN v5.2.9'
    for path, expected in ((args.append_cazy, 'ce715b6fe25de02bf4c7ea7aa1539953a386104f'),
                           (DEFINITIONS, '508a9bb2e6502697739bef0e6225e36144953d49')):
        data = path.read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        assert actual == expected, (path, 'Changed upstream CAZy source')
    sys.path.insert(0, str(args.dbcan.resolve()))
    from dbcan.IO.OverviewGenerator import OverviewGenerator
    from dbcan.configs.base_config import OverviewGeneratorConfig
    from dbcan.constants.databases_constants import S3_DATABASE_RELEASE
    assert S3_DATABASE_RELEASE == 'db_v5-2-9_5-5-2026'
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        generator = OverviewGenerator(OverviewGeneratorConfig(output_dir=temporary))
        indexed = {'dbcan_hmm': {}, 'dbcan_sub': {}, 'diamond': {}}
        for i in range(8):
            gene = f'g{i}'
            if i & 1:
                indexed['dbcan_hmm'][gene] = pd.DataFrame([{'HMM Name': 'GH13', 'Target From': 1, 'Target To': 100}])
            if i & 2:
                indexed['dbcan_sub'][gene] = pd.DataFrame([{'Subfam Name': 'GH13_e1', 'Target From': 1, 'Target To': 100,
                                                           'Subfam EC': '3.2.1.1', 'Substrate': 'starch'}])
            if i & 4:
                indexed['diamond'][gene] = pd.DataFrame([{'CAZy ID': 'GH13'}])
        indexed['dbcan_hmm']['multi'] = pd.DataFrame([{'HMM Name': 'GH13', 'Target From': 1, 'Target To': 100},
                                                    {'HMM Name': 'CBM2', 'Target From': 200, 'Target To': 250}])
        indexed['dbcan_sub']['multi'] = pd.DataFrame([{'Subfam Name': 'GH13_e1', 'Target From': 1, 'Target To': 100,
                                                      'Subfam EC': '3.2.1.1', 'Substrate': 'starch'},
                                                     {'Subfam Name': 'GH5_e2', 'Target From': 300, 'Target To': 400,
                                                      'Subfam EC': '3.2.1.4', 'Substrate': 'cellulose'}])
        genes = [f'g{i}' for i in range(8)] + ['multi']
        overview = generator.aggregate_data(genes, indexed)
        path = directory / 'overview.tsv'
        overview.to_csv(path, sep='\t', index=False)
        normalized = read_overview(path)
        for gene in genes:
            original = overview.loc[overview['Gene ID'] == gene].iloc[0]
            call = normalized[gene]
            assert call['dbcan_tools'] == int(original['#ofTools'])
            assert call['dbcan_recommended'] == (int(original['#ofTools']) >= 2)
            if int(original['#ofTools']) < 2:
                assert call['dbcan_id'] == '-'
        if args.write_fixture:
            fixture = Path(__file__).resolve().parents[1] / 'tests/fixtures'
            fixture.mkdir(exist_ok=True)
            (fixture / 'dbcan_overview.tsv').write_bytes(path.read_bytes())
            (fixture / 'dbcan_normalized.json').write_text(json.dumps(normalized, indent=2) + '\n')

        spec = importlib.util.spec_from_file_location('upstream_append_cazy', args.append_cazy)
        upstream = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(upstream)
        definitions = upstream.read_cazy_definitions(DEFINITIONS)
        raw = pd.read_csv(DEFINITIONS, sep='\t')
        features = sorted(set(feature.strip() for text in raw['function_ids'] for feature in text.split(',')))
        rng = random.Random(723)
        profiles = {'empty': set(), 'all': set(features)}
        profiles.update({f'single_{i}': {feature} for i, feature in enumerate(features)})
        profiles.update({f'random_{i}': {feature for feature in features if rng.random() < density}
                         for i, density in enumerate([0.1, 0.5, 0.9] * 10)})
        # AnnoGuild consumes DRAM1's all-cleavage Boolean product (>=0.6).
        # Build that input explicitly, then execute AnnoGuild's original gate.
        product = []
        distillate = []
        for genome, observed in profiles.items():
            row = {'genome': genome}
            for name, group in raw.groupby('function_name', sort=False):
                row['CAZy: ' + name] = all(observed & {x.strip() for x in text.split(',')} for text in group['function_ids'])
            product.append(row)
            distillate.extend({'genome': genome, 'gene_id': feature, 'call': int(feature in observed)} for feature in features)
        product_path = directory / 'product.tsv'
        pd.DataFrame(product).to_csv(product_path, sep='\t', index=False)
        expected = upstream.pivot_cazy_to_wide(upstream.read_dram_product(definitions, pd.DataFrame(distillate), product_path))
        actual, _ = cazy_calls(profiles)
        expected = expected.set_index('genome')
        comparisons = 0
        for row in actual:
            for name, value in row.items():
                if name != 'genome':
                    assert bool(expected.loc[row['genome'], name]) == value, (row['genome'], name)
                    comparisons += 1
        if args.write_fixture:
            snapshot = {'dbcan_commit': sha, 'annoguild_commit': '291b232907e6568597178e21384ff0f964993a82',
                        'profiles': {genome: sorted(values) for genome, values in profiles.items()},
                        'expected': {genome: {name: bool(value) for name, value in row.items()}
                                     for genome, row in expected.iterrows()}}
            (fixture / 'cazy_parity.json').write_text(json.dumps(snapshot, indent=2) + '\n')
        print(f'PASS: 9 real dbCAN overview profiles and {comparisons:,} AnnoGuild CAZy comparisons')


if __name__ == '__main__':
    main()
