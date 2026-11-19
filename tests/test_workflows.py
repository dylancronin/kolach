"""Snakemake dry runs of the bundled workflows (skipped when snakemake is unavailable)."""
import importlib.util
import json
import os
import re
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

    def dry_run(self, snakefile, config, extra=(), extra_env=None):
        config_path = self.dir / "config.json"
        config_path.write_text(json.dumps(config))
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(SRC), os.environ.get("PYTHONPATH")])))
        env["XDG_CACHE_HOME"] = str(self.dir / "cache")
        env.update(extra_env or {})
        cores = ["--local-cores", "1"] if extra else ["--cores", "1"]
        result = subprocess.run(
            ["snakemake", "-n", "-s", str(WORKFLOWS / snakefile), *cores, *extra,
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


    def test_download_deepkoala_uses_release_resolved_by_cli(self):
        # Under --executor slurm the CLI resolves the release once; submitted jobs
        # must not query GenomeNet again (a dead proxy makes any lookup fail).
        dead = "http://127.0.0.1:9"
        out = self.dry_run("download.smk", {
            "database_dir": str(self.dir / "empty_db"), "databases": ["deepkoala"],
            "deepkoala_release": "202608", "deepkoala_release_resolved": True,
        }, extra_env={"http_proxy": dead, "https_proxy": dead, "HTTP_PROXY": dead, "HTTPS_PROXY": dead})
        self.assertIn("deepkoala/202608/manifest.json", out)

    def test_download_graftm_and_kegg(self):
        out = self.dry_run("download.smk", {
            "database_dir": str(self.dir / "empty_db"), "databases": ["graftm"],
        })
        for target in ("download_graftm", "download_kegg_hierarchy", "graftm/manifest.json"):
            self.assertIn(target, out)

    @unittest.skipUnless(importlib.util.find_spec("snakemake_executor_plugin_slurm"),
                         "Slurm executor plugin not installed")
    def test_annotate_slurm_profile_resources_and_localrules(self):
        out = self.dry_run("annotate.smk", {
            "database_dir": str(self.db), "protein_fasta": str(self.fasta),
            "output_dir": str(self.dir / "out"), "databases": ["kofam", "eggnog", "deepkoala"],
            "add_pathways": True, "deepkoala_resources": {"gpu": 1, "slurm_partition": "gpu"},
        }, extra=["--executor", "slurm", "--jobs", "4", "--workflow-profile", str(WORKFLOWS / "profiles" / "slurm"),
                  "--default-resources", "slurm_account=TEST"])
        blocks = {}
        for block in re.split(r"^(?:local)?rule ", out, flags=re.M)[1:]:
            name, _, body = block.partition(":")
            blocks[name] = body.split("\n[", 1)[0]
        self.assertIn("mem_mb=64000", blocks["annotate_eggnog"])
        self.assertIn("runtime=720", blocks["annotate_eggnog"])
        self.assertIn("mem_mb=32000", blocks["annotate_kofam"])
        self.assertIn("gpu=1", blocks["annotate_deepkoala"])
        self.assertIn("slurm_partition=gpu", blocks["annotate_deepkoala"])
        self.assertNotIn("gpu=", blocks["annotate_kofam"])

    def test_light_rules_are_localrules(self):
        # Light steps run in the controller job; heavy ones must be submitted.
        local = set()
        for path in [WORKFLOWS / "annotate.smk", WORKFLOWS / "download.smk", *WORKFLOWS.glob("rules/*.smk")]:
            for line in path.read_text().splitlines():
                if line.startswith("localrules:"):
                    local.update(name.strip() for name in line.split(":", 1)[1].split(","))
        self.assertEqual(local, {"all", "integrate_annotations", "annotate_pathways", "specialise_genomes",
                                 "download_kegg_hierarchy"})

    def test_specialisations_require_dbcan(self):
        config_path = self.dir / "config.json"
        config_path.write_text(json.dumps({
            "database_dir": str(self.db), "protein_fasta": str(self.fasta),
            "output_dir": str(self.dir / "out"), "databases": ["kofam"], "add_specialisations": True,
        }))
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(filter(None, [str(SRC), os.environ.get("PYTHONPATH")])))
        result = subprocess.run(
            ["snakemake", "-n", "-s", str(WORKFLOWS / "annotate.smk"), "--cores", "1",
             "--configfile", str(config_path), "--directory", str(self.dir)],
            capture_output=True, text=True, env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires dbcan", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
