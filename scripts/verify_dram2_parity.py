"""Compare Kolach with the actual pinned DRAM2 parser (no parser stubs).

Usage: PYTHONPATH=src python scripts/verify_dram2_parity.py DRAM_CHECKOUT
Requires upstream's lark, networkx, numpy and polars. Fails on changed source
commits/assets or any unlisted upstream failure or numeric/Boolean mismatch.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys

import numpy as np
import polars as pl
from kolach.metabolism.dram2 import DATA, evaluate, load_rules


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('dram', type=Path)
    parser.add_argument('--write-fixture', action='store_true')
    parser.add_argument('--dram-viz', type=Path, help='Also verify product assets against their source checkout')
    args = parser.parse_args()
    root = args.dram / 'bin/rule_parser'
    manifest = json.loads((DATA / 'manifest.json').read_text())
    sha = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    assert sha == manifest['common.tsv']['commit'], 'Wrong upstream parser checkout'
    for name, source in manifest.items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == source['sha256'], name
        checkout = root if source['path'].startswith('src/') or name == 'LICENSE.GPL-3.0' else args.dram if name == 'traits.tsv' else args.dram_viz
        if checkout is not None:
            actual_sha = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
            assert actual_sha == source['commit'], (name, 'Wrong source commit')
            data = subprocess.check_output(['git', '-C', str(checkout), 'show', 'HEAD:' + source['path']])
            assert hashlib.sha256(data).hexdigest() == source['sha256'], (name, 'Changed upstream asset')
    spec = importlib.util.spec_from_file_location('upstream_rules', root / 'src/rules.py')
    upstream = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = upstream
    spec.loader.exec_module(upstream)
    rng = random.Random(723)
    comparisons = 0
    fixture = {'parser_commit': sha, 'tables': {}, 'known_upstream_errors': {}}
    for path in sorted(DATA.glob('*.tsv')):
        if path.name == 'common.tsv':
            continue
        compiled = upstream.CompiledRules.from_rules(rules_path=path, common_rules_path=DATA / 'common.tsv', allow_visualize_functions=True)
        features = sorted(compiled.needed_features)
        # Empty, complete, sparse/dense sets, plus every single feature. Counts
        # include repeated genes; at_least COUNTS currently coerces to bool.
        profiles = [{}, dict.fromkeys(features, 1)]
        profiles += [{f: rng.randrange(1, 5) for f in features if rng.random() < density}
                     for density in (0.05, 0.2, 0.5, 0.8) for _ in range(8)]
        profiles += [{feature: 1} for feature in features]
        samples = [str(i) for i in range(len(profiles))]
        annotations = pl.DataFrame({'genome': samples, 'taxonomy': ['g__Methanosarcina'] * len(samples),
                                    'kegg_description': ['nitrate reductase'] * len(samples),
                                    'heme_regulatory_motif_count': [4] * len(samples)})
        present_map = {f: np.array([p.get(f, 0) for p in profiles]) for f in features}
        ev = upstream.Evaluator(samples, present_map, 'genome', annotations)
        local = {row['name']: node for row, node in load_rules(path)}
        expected = {}
        errors = {}
        for name, expression in compiled.rules.items():
            if name == '':
                # Upstream includes the last common alias as an unnamed output
                # because it does not normalize common-sheet empty names.
                errors[name] = 'Unnamed common-sheet output omitted'
                continue
            try:
                result = ev.eval_bool(expression)
            except TypeError as error:
                # Pinned upstream passes masks= to not_(), which has no masks
                # parameter. Only this precise defect is an allowed exclusion.
                if name != 'probable-iron-reduction' or "unexpected keyword argument 'masks'" not in str(error):
                    raise
                errors[name] = str(error)
                continue
            wanted = result.to_dicts() if isinstance(result, pl.DataFrame) else [bool(x) for x in result]
            captured = []
            for i, profile in enumerate(profiles):
                rows = [{'taxonomy': 'g__Methanosarcina', 'kegg_description': 'nitrate reductase', 'heme_regulatory_motif_count': '4'}]
                actual = evaluate(local[name], profile, rows)
                target = wanted[i]
                if isinstance(actual, dict):
                    target = {key: value for key, value in target.items() if key != 'genome'}
                else:
                    actual = bool(actual)
                assert actual == target, (path.name, name, i, actual, target)
                captured.append(target)
                comparisons += 1
            expected[name] = captured[:6]
        fixture['tables'][path.name] = {'profiles': profiles[:6], 'expected': expected}
        fixture['known_upstream_errors'][path.name] = errors
        print(f'{path.name}: {len(compiled.rules)} rules, {len(profiles)} profiles, {len(errors)} known upstream errors')
    if args.write_fixture:
        target = Path(__file__).resolve().parents[1] / 'tests/fixtures/dram2_parity.json'
        target.write_text(json.dumps(fixture, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(f'PASS: {comparisons:,} independent comparisons')


if __name__ == '__main__':
    main()
