"""Snakemake dry runs of the bundled workflows (skipped when snakemake is unavailable)."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

WORKFLOWS = Path(__file__).resolve().parents[1] / "src" / "kolach" / "workflows"
SRC = WORKFLOWS.parents[1]


@unittest.skipUnless(shutil.which("snakemake"), "snakemake not installed")
class TestWorkflowDryRun(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        # Path with characters that `--config key=value` parsing used to mangle.
        self.db = self.dir / "db: x,y"
        (self.db / "kofam").mkdir(parents=True)
        (self.db / "kofam" / "profiles.hmm").write_text("")
        (self.db / "kofam" / "ko_list").write_text("")
        (self.db / "eggnog").mkdir()
        (self.db / "eggnog" / ".download_complete").write_text("")
        self.fasta = self.dir / "proteins.faa"
        self.fasta.write_text(">g1\nMK\n")

    def tearDown(self):
        self.tmp.cleanup()

    def dry_run(self, snakefile, config):
        config_path = self.dir / "config.json"
        config_path.write_text(json.dumps(config))
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(SRC), os.environ.get("PYTHONPATH")])))
        env["XDG_CACHE_HOME"] = str(self.dir / "cache")
        result = subprocess.run(
            ["snakemake", "-n", "-s", str(WORKFLOWS / snakefile), "--cores", "1",
             "--configfile", str(config_path), "--directory", str(self.dir)],
            capture_output=True, text=True, env=env,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def test_annotate_all_methods_with_pathways(self):
        out = self.dry_run("annotate.smk", {
            "database_dir": str(self.db), "protein_fasta": str(self.fasta),
            "output_dir": str(self.dir / "out"), "databases": ["kofam", "eggnog", "deepkoala"],
            "threads": 1, "add_pathways": True, "skip_bitscore_heuristic": False,
            "detail": True, "deepkoala_release": "latest",
        })
        for rule in ("annotate_kofam", "annotate_eggnog", "annotate_deepkoala",
                     "integrate_annotations", "annotate_pathways"):
            self.assertIn(rule, out)

    def test_annotate_kofam_only(self):
        out = self.dry_run("annotate.smk", {
            "database_dir": str(self.db), "protein_fasta": str(self.fasta),
            "output_dir": str(self.dir / "out"), "databases": ["kofam"],
        })
        self.assertIn("annotate_kofam", out)
        self.assertNotIn("annotate_pathways", out)

    def test_download_kofam_eggnog_and_kegg(self):
        empty = self.dir / "empty_db"
        out = self.dry_run("download.smk", {"database_dir": str(empty), "databases": ["kofam", "eggnog"]})
        for rule in ("download_kofam", "download_eggnog", "download_kegg_hierarchy"):
            self.assertIn(rule, out)

    def test_download_dbcan_and_kegg(self):
        out = self.dry_run("download.smk", {
            "database_dir": str(self.dir / "empty_db"), "databases": ["dbcan"],
        })
        for target in ("download_dbcan", "download_kegg_hierarchy", "manifest.json",
                       "CAZy.dmnd", "dbCAN.hmm", "dbCAN-sub.hmm", "fam-substrate-mapping.tsv"):
            self.assertIn(target, out)


if __name__ == "__main__":
    unittest.main()
