"""Unit tests for the method wrappers (KOfam, eggNOG, DeepKOALA) without external tools."""
import argparse
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from urllib.error import URLError

from kolach import download
from kolach.methods import deepkoala, eggnog, kofam


def _profiles():
    return {
        "K00001": {"knum": "K00001", "threshold": 100.0, "score_type": "full", "definition": "one"},
        "K00002": {"knum": "K00002", "threshold": 100.0, "score_type": "domain", "definition": "two"},
        "K00003": {"knum": "K00003", "threshold": 200.0, "score_type": "full", "definition": "three"},
    }


def _hit(gene, ko, bit, ev, dbit=None, dev=None):
    return {"gene_id": gene, "ko": ko, "bit_score": bit, "e_value": ev,
            "domain_bit_score": bit if dbit is None else dbit,
            "domain_e_value": ev if dev is None else dev}


class TestKofamAssignHits(unittest.TestCase):
    def assign(self, hits, **kwargs):
        rows = kofam.assign_hits(hits, _profiles(), **kwargs)
        return [(row["gene_id"], row["ko"], row["assignment"]) for row in rows]

    def test_threshold_is_inclusive(self):
        self.assertEqual(self.assign([_hit("g", "K00001", 100.0, 1.0)]), [("g", "K00001", "threshold")])
        self.assertEqual(self.assign([_hit("g", "K00001", 99.9, 1.0)]), [])

    def test_domain_profiles_use_domain_scores(self):
        hits = [_hit("g", "K00002", 50.0, 1.0, dbit=150.0, dev=1e-30)]
        self.assertEqual(self.assign(hits), [("g", "K00002", "threshold")])
        hits = [_hit("g", "K00002", 150.0, 1e-30, dbit=50.0, dev=1.0)]
        self.assertEqual(self.assign(hits), [])

    def test_rescue_requires_strictly_greater_score_and_evalue(self):
        # 0.75 * 100 = 75: strictly greater is required for rescue.
        self.assertEqual(self.assign([_hit("g", "K00001", 75.0, 1e-10)]), [])
        self.assertEqual(self.assign([_hit("g", "K00001", 75.1, 1e-10)]), [("g", "K00001", "rescued")])
        # E-value cutoff is inclusive.
        self.assertEqual(self.assign([_hit("g", "K00001", 80.0, 1e-5)]), [("g", "K00001", "rescued")])
        self.assertEqual(self.assign([_hit("g", "K00001", 80.0, 2e-5)]), [])

    def test_rescue_only_for_unannotated_genes_with_unique_candidate(self):
        # Gene already annotated above threshold: no rescue.
        hits = [_hit("g", "K00001", 150.0, 1e-30), _hit("g", "K00003", 180.0, 1e-30)]
        self.assertEqual(self.assign(hits), [("g", "K00001", "threshold")])
        # Two distinct rescue candidates: ambiguous, no rescue.
        hits = [_hit("g", "K00001", 80.0, 1e-10), _hit("g", "K00003", 160.0, 1e-10)]
        self.assertEqual(self.assign(hits), [])
        # The same KO hit twice is still a unique candidate; best E-value row is kept.
        hits = [_hit("g", "K00001", 80.0, 1e-10), _hit("g", "K00001", 90.0, 1e-20)]
        rows = kofam.assign_hits(hits, _profiles())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["bit_score"], 90.0)

    def test_rescue_can_be_disabled_and_tuned(self):
        hits = [_hit("g", "K00001", 80.0, 1e-10)]
        self.assertEqual(self.assign(hits, rescue=False), [])
        self.assertEqual(self.assign(hits, fraction=0.9), [])
        self.assertEqual(self.assign(hits, evalue=1e-12), [])

    def test_invalid_parameters(self):
        with self.assertRaises(ValueError):
            kofam.assign_hits([], _profiles(), fraction=1.5)
        with self.assertRaises(ValueError):
            kofam.assign_hits([], _profiles(), evalue=-1)

    def test_output_sorted(self):
        hits = [_hit("b", "K00001", 150.0, 1e-30), _hit("a", "K00003", 250.0, 1e-30),
                _hit("a", "K00001", 150.0, 1e-30)]
        self.assertEqual(self.assign(hits), [("a", "K00001", "threshold"), ("a", "K00003", "threshold"),
                                             ("b", "K00001", "threshold")])


class TestKofamParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_parse_hits(self):
        table = self.dir / "hits.tbl"
        table.write_text(
            "# header\n\n"
            "gene1 - K00001 - 1.5e-30 105.2 0.0 2e-29 101.0 0.0 1.0 1 1 1 1 1 1 1 some target description\n"
        )
        self.assertEqual(list(kofam.parse_hits(table)), [{
            "gene_id": "gene1", "ko": "K00001", "e_value": 1.5e-30, "bit_score": 105.2,
            "domain_e_value": 2e-29, "domain_bit_score": 101.0,
        }])
        table.write_text("gene1 - K00001 - 1.5e-30\n")
        with self.assertRaises(ValueError):
            list(kofam.parse_hits(table))

    def test_count_proteins(self):
        fasta = self.dir / "p.faa"
        fasta.write_text(">a desc\nMK\nLL\n\n>b\nM\n")
        self.assertEqual(kofam.count_proteins(fasta), 2)
        for bad in (">a\nM\n>a\nM\n", ">a\n>b\nM\n", "MK\n>a\nM\n", ">a\nM\n>b\n", "", ">\nM\n"):
            fasta.write_text(bad)
            with self.assertRaises(ValueError, msg=bad):
                kofam.count_proteins(fasta)

    def test_read_profiles_skips_missing_thresholds(self):
        ko_list = self.dir / "ko_list"
        ko_list.write_text(
            "knum\tthreshold\tscore_type\tdefinition\n"
            "K00001\t100.5\tfull\tone\n"
            "K00002\t-\t-\ttwo\n"
        )
        profiles = kofam.read_profiles(ko_list)
        self.assertEqual(set(profiles), {"K00001"})
        self.assertEqual(profiles["K00001"]["threshold"], 100.5)

    def test_annotate_checks_hmmsearch_and_builds_command(self):
        db = self.dir / "db" / "kofam"
        db.mkdir(parents=True)
        (db / "ko_list").write_text("knum\tthreshold\tscore_type\tdefinition\nK00001\t100\tfull\tone\n")
        fasta = self.dir / "p.faa"
        fasta.write_text(">g1\nMK\n>g2\nMK\n")
        args = argparse.Namespace(database_dir=str(self.dir / "db"), protein_fasta=str(fasta), threads=2)

        with mock.patch("shutil.which", return_value=None):
            with self.assertRaisesRegex(FileNotFoundError, "hmmsearch"):
                kofam.annotate(args, io.StringIO())

        calls = []

        def fake_run(cmd, check):
            calls.append(cmd)
            Path(cmd[cmd.index("--tblout") + 1]).write_text(
                "g1 - K00001 - 1e-30 150.0 0.0 1e-30 150.0 0.0 1 1 1 1 1 1 1 1 d\n")

        out = io.StringIO()
        with mock.patch("shutil.which", return_value="/usr/bin/hmmsearch"), \
                mock.patch("subprocess.run", fake_run):
            kofam.annotate(args, out)
        cmd = calls[0]
        self.assertEqual(cmd[:7], ["hmmsearch", "--cpu", "2", "--noali", "-o", os.devnull, "--tblout"])
        self.assertEqual(cmd[8:12], ["-Z", "2", "--domZ", "2"])
        self.assertIn("g1\tK00001\tthreshold\tfull\t100.0\t150.0", out.getvalue())


class TestEggnog(unittest.TestCase):
    def test_optional(self):
        for value in (None, "", "None", "none", " NONE "):
            self.assertIsNone(eggnog._optional(value))
        for value in ("fast", 0, 2.5, "/tmp"):
            self.assertEqual(eggnog._optional(value), value)

    def test_package_data_dir_points_next_to_package(self):
        spec = mock.Mock(origin=str(Path("/env/lib/site/eggnogmapper/__init__.py")))
        with mock.patch("importlib.util.find_spec", return_value=spec):
            self.assertEqual(eggnog._package_data_dir(), Path("/env/lib/site/data").resolve())
        with mock.patch("importlib.util.find_spec", return_value=None), \
                mock.patch("sysconfig.get_paths", return_value={"purelib": "/purelib"}):
            self.assertEqual(eggnog._package_data_dir(), Path("/purelib") / "data")

    def test_annotate_builds_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "db" / "eggnog").mkdir(parents=True)
            (tmp / "db" / "eggnog" / "eggnog.db").write_text("")
            calls = []

            def fake_run(cmd, env):
                calls.append(cmd)
                out_dir = Path(cmd[cmd.index("--output_dir") + 1])
                (out_dir / "eggnog_res.emapper.annotations").write_text(
                    "## meta\n#query\tKEGG_ko\ng1\tko:K00001\n## end\n")
                return mock.Mock(returncode=0)

            args = argparse.Namespace(
                database_dir=str(tmp / "db"), protein_fasta=str(tmp / "p.faa"), threads=3,
                eggnog_mode="diamond", eggnog_sensmode="none", eggnog_temp_dir=None,
                eggnog_dmnd_block_size="2.5", eggnog_dmnd_index_chunks=4,
            )
            with mock.patch("shutil.which", return_value="/bin/emapper.py"), \
                    mock.patch("subprocess.run", fake_run), \
                    mock.patch.object(eggnog, "ensure_package_data_dir"):
                eggnog.annotate(args, tmp / "out.tsv")
            cmd = calls[0]
            self.assertNotIn("--dmnd_sensmode", cmd)
            self.assertEqual(cmd[cmd.index("--dmnd_block_size") + 1], "2.5")
            self.assertEqual(cmd[cmd.index("--dmnd_index_chunks") + 1], "4")
            self.assertEqual(cmd[cmd.index("--cpu") + 1], "3")
            self.assertEqual((tmp / "out.tsv").read_text(), "query\tKEGG_ko\ng1\tko:K00001\n")

    def test_annotate_requires_emapper(self):
        args = argparse.Namespace(database_dir="db", protein_fasta="p.faa")
        with mock.patch("shutil.which", return_value=None):
            with self.assertRaisesRegex(FileNotFoundError, "emapper.py"):
                eggnog.annotate(args, "out.tsv")


class TestDeepkoalaDownloads(unittest.TestCase):
    def test_download_atomic_keeps_previous_file_on_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "weights_full.pt"
            dest.write_bytes(b"previous")

            class Broken(io.BytesIO):
                def read(self, size=-1):
                    raise URLError("connection reset")

            with mock.patch.object(download, "urlopen", return_value=Broken()):
                with self.assertRaises(URLError):
                    download.download_atomic("https://example.invalid/x", dest)
            self.assertEqual(dest.read_bytes(), b"previous")
            self.assertEqual(sorted(p.name for p in Path(tmp).iterdir()), ["weights_full.pt"])

            with mock.patch.object(download, "urlopen", return_value=io.BytesIO(b"new weights")):
                digest = download.download_atomic("https://example.invalid/x", dest)
            self.assertEqual(dest.read_bytes(), b"new weights")
            self.assertEqual(len(digest), 64)

    def test_resolve_release(self):
        with mock.patch.object(deepkoala, "get_available_releases", return_value=["202501", "202608"]):
            self.assertEqual(deepkoala.resolve_release("latest"), "202608")
            self.assertEqual(deepkoala.resolve_release(202501), "202501")
            with self.assertRaisesRegex(ValueError, "not found"):
                deepkoala.resolve_release("209912")

    def test_check_dependencies(self):
        with mock.patch("importlib.util.find_spec", return_value=None):
            with self.assertRaisesRegex(ImportError, "deepkoala"):
                deepkoala.check_dependencies()


if __name__ == "__main__":
    unittest.main()
