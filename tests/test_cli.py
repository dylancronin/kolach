"""Tests for kolach.cli argument parsing and workflow handoff."""
import contextlib
import io
import json
from pathlib import Path
import unittest
from unittest import mock

from kolach import cli


def parse(*argv):
    return cli.build_parser().parse_args(list(argv))


class TestParser(unittest.TestCase):
    def assert_rejected(self, *argv):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parse(*argv)

    def test_positive_ints(self):
        base = ["annotate", "--protein-fasta", "p.faa", "--database-dir", "db"]
        for flag in ("--threads", "--batch-size"):
            for bad in ("0", "-2", "x"):
                self.assert_rejected(*base, flag, bad)
        args = parse(*base, "--threads", "8", "--batch-size", "32")
        self.assertEqual((args.threads, args.batch_size), (8, 32))

    def test_shared_integration_options_and_defaults(self):
        annotate = parse("annotate", "--protein-fasta", "p.faa", "--database-dir", "db")
        integrate = parse("integrate", "--output-file", "out.tsv")
        for args in (annotate, integrate):
            self.assertEqual(args.heuristic_bitscore_fraction, 0.75)
            self.assertEqual(args.heuristic_e_value, 1e-5)
            self.assertEqual(args.eggnog_min_bitscore, 60.0)
            self.assertEqual(args.eggnog_max_evalue, 1e-5)
            self.assertEqual(args.eggnog_filter_multi, "disambiguate")
            self.assertEqual(args.conflict_strategy, "multiple")
            self.assertFalse(args.add_pathways)
        self.assert_rejected("integrate", "--output-file", "o", "--heuristic-bitscore-fraction", "1.5")
        self.assert_rejected("annotate", "--protein-fasta", "p", "--database-dir", "d", "--eggnog-max-evalue", "-1")

    def test_integrate_add_pathways_requires_database_dir(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            cli.main(["integrate", "--output-file", "o.tsv", "--add-pathways"])


class TestWorkflowConfig(unittest.TestCase):
    def test_annotate_config_is_typed_and_omits_unset_options(self):
        args = parse(
            "annotate", "--protein-fasta", "dir: with/odd,chars/p.faa", "--database-dir", "db",
            "--databases", "kofam", "eggnog", "--threads", "4", "--no-detail",
            "--deepkoala-release", "202608", "--eggnog-dmnd-block-size", "2.5",
        )
        config = cli.annotate_config(args)
        self.assertEqual(config["protein_fasta"], "dir: with/odd,chars/p.faa")
        self.assertEqual(config["databases"], ["kofam", "eggnog"])
        self.assertEqual(config["threads"], 4)
        self.assertIs(config["detail"], False)
        self.assertIs(config["skip_bitscore_heuristic"], False)
        self.assertEqual(config["deepkoala_release"], "202608")
        self.assertEqual(config["eggnog_dmnd_block_size"], 2.5)
        self.assertNotIn("eggnog_sensmode", config)
        self.assertNotIn("eggnog_temp_dir", config)
        self.assertEqual(json.loads(json.dumps(config)), config)

    def test_run_workflow_passes_configfile(self):
        seen = {}

        def fake_run(cmd):
            path = Path(cmd[cmd.index("--configfile") + 1])
            seen["cmd"] = cmd
            seen["config"] = json.loads(path.read_text(encoding="utf-8"))
            seen["path"] = path
            return mock.Mock(returncode=3)

        with mock.patch("shutil.which", return_value="/bin/snakemake"), \
                mock.patch("subprocess.run", fake_run):
            code = cli.run_workflow("annotate.smk", {"protein_fasta": "a: b", "threads": 2}, cores=2)
        self.assertEqual(code, 3)
        self.assertEqual(seen["config"], {"protein_fasta": "a: b", "threads": 2})
        self.assertNotIn("--config", seen["cmd"])
        self.assertEqual(seen["cmd"][seen["cmd"].index("--cores") + 1], "2")
        self.assertTrue(seen["cmd"][2].endswith("annotate.smk"))
        self.assertFalse(seen["path"].exists())

    def test_annotate_preflight_stops_before_workflow(self):
        with mock.patch("shutil.which", return_value=None), \
                mock.patch.object(cli, "run_workflow") as run, \
                contextlib.redirect_stderr(io.StringIO()) as err, \
                self.assertRaises(SystemExit) as exit_:
            cli.main(["annotate", "--protein-fasta", "p.faa", "--database-dir", "db", "--databases", "kofam"])
        self.assertEqual(exit_.exception.code, 1)
        self.assertIn("hmmsearch", err.getvalue())
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
