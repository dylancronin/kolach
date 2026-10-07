import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from kolach.metabolism.dram2 import DATA, annotation_features, evaluate, load_rules, parse_rule, run_dram2

ROOT = Path(__file__).resolve().parents[1]


class TestDram2(unittest.TestCase):
    def test_upstream_snapshots_and_asset_hashes(self):
        for name, metadata in json.loads((DATA / 'manifest.json').read_text()).items():
            self.assertEqual(hashlib.sha256((DATA / name).read_bytes()).hexdigest(), metadata['sha256'])
        fixture = json.loads((ROOT / 'tests/fixtures/dram2_parity.json').read_text())
        rows = [{'taxonomy': 'g__Methanosarcina', 'kegg_description': 'nitrate reductase', 'heme_regulatory_motif_count': '4'}]
        for table, tests in fixture['tables'].items():
            rules = {row['name']: node for row, node in load_rules(DATA / table)}
            for name, expected in tests['expected'].items():
                for counts, result in zip(tests['profiles'], expected):
                    actual = evaluate(rules[name], counts, rows)
                    self.assertEqual(actual if isinstance(actual, dict) else bool(actual), result, (table, name))

    def test_boolean_steps_subunits_and_thresholds(self):
        counts = {'A': 1, 'B': 1}
        self.assertFalse(evaluate(parse_rule('A,B,C'), counts))
        self.assertTrue(evaluate(parse_rule('A | [B & C]'), counts))
        self.assertTrue(evaluate(parse_rule('percent(66, [A,B,C])'), counts))
        self.assertFalse(evaluate(parse_rule('percent(67, [A,B,C])'), counts))
        self.assertEqual(evaluate(parse_rule('path_steps(A & B,C)'), counts),
                         {'steps': 2, 'steps_present': 1, 'coverage_percentage': 0.5})
        self.assertEqual(evaluate(parse_rule('path_subunits(A & B & C)'), counts),
                         {'subunits': 3, 'subunits_present': 2, 'coverage_percentage': 2 / 3})
        # Explicitly preserve the upstream COUNTS coercion.
        self.assertFalse(evaluate(parse_rule('at_least(2,COUNTS,[A,B])'), {'A': 10}))
        for text in ('A | B & C', '[A', 'A+', 'percent(50,[A,B]) trailing'):
            with self.assertRaises(ValueError):
                parse_rule(text)

    def test_taxonomy_and_gene_filters(self):
        rows = [{'taxonomy': 'g__Methanosarcina', 'description': 'nitrate reductase', 'hemes': '9'},
                {'taxonomy': '', 'description': 'cytochrome', 'hemes': '4'}]
        self.assertTrue(evaluate(parse_rule('tax(g__Methanosarcina)'), {}, rows))
        self.assertFalse(evaluate(parse_rule('tax(g__Methanosarcina)'), {}, []))
        rule = 'not(filter_contains(description,"nitrate reductase")) -> column_count_values(hemes,ge,4,ge,2)'
        self.assertFalse(evaluate(parse_rule(rule), {}, rows))
        self.assertTrue(evaluate(parse_rule(rule.replace('ge,2)', 'ge,1)')), {}, rows))
        with self.assertRaises(ValueError):
            evaluate(parse_rule('filter_compare(missing,ge,2)'), {}, rows)

    def test_candidates_are_not_evidence(self):
        features = annotation_features({'accepted_ko': 'K00001', 'alternative_kos': 'K00002',
                                        'kofam_id': 'K00003', 'dbcan_id': 'GH13;GH5'})
        self.assertEqual(features, {'K00001', 'GH13', 'GH5'})

    def test_cli_output_grouping_and_taxonomy_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            annotation = path / 'annotations.tsv'
            annotation.write_text('genome\tgene_id\taccepted_ko\n001\tg1\tK00399\n002\tg1\t-\n')
            result = run_dram2(annotation, path / 'out')
            text = result.read_text()
            self.assertIn('001\ttraits\tmethanogen or ANME', text)
            self.assertNotIn('\ttraits\thydrogenotrophic-methanogen', text)
            self.assertIn('002\tproduct\t', text)
            annotation.write_text('genome\tgene_id\taccepted_ko\ttaxonomy\n001\tg1\tK00399\tg__Methanosarcina\n')
            text = run_dram2(annotation, path / 'out').read_text()
            self.assertNotIn('\ttraits\tmethanogen or ANME', text)
            self.assertIn('\ttraits\thydrogenotrophic-methanogen', text)


if __name__ == '__main__':
    unittest.main()
