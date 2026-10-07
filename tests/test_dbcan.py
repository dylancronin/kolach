import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from kolach.metabolism.dbcan import DATABASE_FILES, DBCAN_COLUMNS, cazy_calls, download_dbcan, families, make_cazy_product, merge_overview, read_overview, run_dbcan
from kolach.metabolism.specialise import specialise
from kolach.metabolism.evaluate import read_tsv

FIXTURES = Path(__file__).parent / 'fixtures'


class TestDbcan(unittest.TestCase):
    def test_annoguild_signature_gate_snapshots(self):
        fixture = json.loads((FIXTURES / 'cazy_parity.json').read_text())
        actual, _ = cazy_calls({genome: set(values) for genome, values in fixture['profiles'].items()})
        for row in actual:
            genome = row.pop('genome')
            self.assertEqual(row, fixture['expected'][genome], genome)

    def test_real_upstream_overview_snapshot(self):
        self.assertEqual(read_overview(FIXTURES / 'dbcan_overview.tsv'), json.loads((FIXTURES / 'dbcan_normalized.json').read_text()))
        self.assertEqual(read_overview(FIXTURES / 'dbcan_overview.tsv')['multi']['dbcan_family_ec'],
                         'GH13;EC:3.2.1.1|GH5;EC:3.2.1.4')

    def test_family_boundaries_and_signatures(self):
        self.assertEqual(families('GH13_e1(1-100)+CBM2(200-250)|GH5'), {'GH13_e1', 'CBM2', 'GH5'})
        with self.assertRaises(ValueError):
            families('GH13unexpected')
        calls, _ = cazy_calls({'missing_oligo': {'GH13'}, 'starch': {'GH13', 'GH15'},
                               'low_signature': {'GH23', 'CE4'}, 'chitin': {'GH18', 'GH20'}, 'empty': set()})
        rows = {row['genome']: row for row in calls}
        self.assertFalse(rows['missing_oligo']['CAZy-Starch'])
        self.assertTrue(rows['starch']['CAZy-Starch'])
        self.assertFalse(rows['low_signature']['CAZy-Chitin'])
        self.assertTrue(rows['chitin']['CAZy-Chitin'])
        self.assertFalse(any(value for key, value in rows['empty'].items() if key != 'genome'))

    def test_join_inventory_candidates_and_metabolism(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            annotation = path / 'input.tsv'
            annotation.write_text('gene_id\taccepted_ko\talternative_kos\ng0\t-\tK05909\ng1\t-\t-\ng3\t-\t-\ng7\t-\t-\nmulti\t-\t-\nmissing\t-\t-\n')
            # A dbCAN file with extra genes is rejected, not silently trimmed.
            with self.assertRaises(ValueError):
                merge_overview(annotation, FIXTURES / 'dbcan_overview.tsv', path / 'merged.tsv')
            annotation.write_text('gene_id\taccepted_ko\talternative_kos\n' + ''.join(f'{gene}\t-\tK05909\n' for gene in [*[f'g{i}' for i in range(8)], 'multi', 'missing']))
            merged = merge_overview(annotation, FIXTURES / 'dbcan_overview.tsv', path / 'merged.tsv')
            text = merged.read_text()
            self.assertIn('missing\t-\tK05909', text)
            self.assertIn('dbcan_hmm_candidates', text)
            product = make_cazy_product(merged, path, genome_id='001')
            result = specialise(merged, path / 'metabolism', genome_id='001', cazy_product=product)
            self.assertIn('assessed', result['assessment'].read_text())
            # An ambiguous alternative KO must not create polyphenolic evidence.
            _, pathways = read_tsv(path / 'kolach_cazy_pathways.tsv')
            self.assertEqual(next(row['called'] for row in pathways if row['pathway'] == 'CAZy-Polyphenolics'), 'False')
            with self.assertRaises(ValueError):
                merge_overview(merged, FIXTURES / 'dbcan_overview.tsv', path / 'again.tsv')

    def test_download_hash_verification_and_interruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('kolach.metabolism.dbcan.urlopen', side_effect=lambda *a, **k: io.BytesIO(b'reference')) as network:
                manifest = download_dbcan(tmp)
                self.assertEqual(network.call_count, 4)
                download_dbcan(tmp)
                self.assertEqual(network.call_count, 4)
                (Path(tmp) / DATABASE_FILES[0]).write_bytes(b'corrupt')
                with self.assertRaises(ValueError):
                    download_dbcan(tmp)
            self.assertEqual(json.loads(manifest.read_text())['release'], 'db_v5-2-9_5-5-2026')
        with tempfile.TemporaryDirectory() as tmp:
            with patch('kolach.metabolism.dbcan.urlopen', side_effect=OSError('interrupted')):
                with self.assertRaises(OSError):
                    download_dbcan(tmp)
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_runner_is_fresh_version_pinned_and_failure_safe(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            for name in DATABASE_FILES:
                (path / name).write_bytes(b'fixture reference')
            fasta = path / 'input.faa'
            fasta.write_text('>g0\nM\n')
            output = path / 'out'
            output.mkdir()
            (output / 'overview.tsv').write_text('stale')
            def successful(command, check):
                directory = Path(command[command.index('--output_dir') + 1])
                self.assertNotEqual(directory, output)
                self.assertEqual(command[command.index('--methods') + 1], 'diamond,hmm,dbCANsub')
                (directory / 'overview.tsv').write_bytes((FIXTURES / 'dbcan_overview.tsv').read_bytes())
            with patch('kolach.metabolism.dbcan.shutil.which', return_value='run_dbcan'), patch('kolach.metabolism.dbcan.subprocess.check_output', return_value='dbCAN version: 5.2.9'), patch('kolach.metabolism.dbcan.subprocess.run', side_effect=successful):
                result = run_dbcan(fasta, path, output, 2)
                self.assertNotEqual(result.read_text(), 'stale')
            before = result.read_bytes()
            with patch('kolach.metabolism.dbcan.shutil.which', return_value='run_dbcan'), patch('kolach.metabolism.dbcan.subprocess.check_output', return_value='dbCAN version: 5.2.9'), patch('kolach.metabolism.dbcan.subprocess.run', side_effect=subprocess.CalledProcessError(1, 'run_dbcan')):
                with self.assertRaises(subprocess.CalledProcessError):
                    run_dbcan(fasta, path, output)
            self.assertEqual(result.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
