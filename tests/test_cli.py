"""Tests for kolach.cli argument parsing and workflow handoff."""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
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

    def test_annotate_rejects_duplicate_fasta_ids_before_workflow(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "p.faa"
            fasta.write_text(">g1 a\nM\n>g2\nM\n>g1 b\nM\n")
            with mock.patch.object(cli, "check_annotation_dependencies"), \
                    mock.patch.object(cli, "run_workflow") as run, \
                    contextlib.redirect_stderr(io.StringIO()) as err, \
                    self.assertRaises(SystemExit) as exit_:
                cli.main(["annotate", "--protein-fasta", str(fasta), "--database-dir", "db",
                          "--databases", "eggnog"])
        self.assertEqual(exit_.exception.code, 1)
        self.assertIn("1 duplicate gene ID(s): g1", err.getvalue())
        run.assert_not_called()


class TestSlurmExecutor(unittest.TestCase):
    BASE = ["annotate", "--protein-fasta", "p.faa", "--database-dir", "db"]

    def run_main(self, argv):
        """Run cli.main with the workflow mocked; return (run_workflow mock, stderr)."""
        with mock.patch.object(cli, "check_annotation_dependencies"), \
                mock.patch("kolach.integrate.read_fasta_ids", return_value=["g1"]), \
                mock.patch.object(cli, "run_workflow", return_value=0) as run, \
                contextlib.redirect_stderr(io.StringIO()) as err:
            cli.main(argv)
        return run, err.getvalue()

    def assert_main_rejected(self, argv, message):
        with mock.patch.object(cli, "run_workflow") as run, \
                contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit):
            cli.main(argv)
        self.assertIn(message, err.getvalue())
        run.assert_not_called()

    def test_local_is_default_and_unchanged(self):
        run, _ = self.run_main(self.BASE + ["--threads", "3"])
        self.assertEqual(run.call_args.kwargs["slurm_args"], None)
        self.assertEqual(run.call_args.kwargs["cores"], 3)

    def test_slurm_requires_account(self):
        self.assert_main_rejected(self.BASE + ["--executor", "slurm"], "requires --slurm-account")

    def test_slurm_options_require_slurm_executor(self):
        self.assert_main_rejected(self.BASE + ["--slurm-account", "A"], "require(s) --executor slurm")
        self.assert_main_rejected(["download", "--database-dir", "db", "--databases", "kofam",
                                   "--slurm-partition", "p"], "require(s) --executor slurm")

    def test_slurm_profile_must_contain_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assert_main_rejected(self.BASE + ["--executor", "slurm", "--slurm-account", "A",
                                                   "--slurm-profile", tmp], "containing config.yaml")

    def test_annotate_slurm_arguments(self):
        run, _ = self.run_main(self.BASE + ["--executor", "slurm", "--slurm-account", "PDS0325",
                                            "--slurm-partition", "cpu", "--jobs", "4", "--threads", "2",
                                            "--output-dir", "out"])
        args = run.call_args.kwargs["slurm_args"]
        self.assertEqual(args[:4], ["--executor", "slurm", "--jobs", "4"])
        self.assertEqual(Path(args[args.index("--workflow-profile") + 1]), cli.SLURM_PROFILE)
        self.assertTrue((cli.SLURM_PROFILE / "config.yaml").is_file())
        # Quoted so Snakemake, which evaluates resource values as Python, keeps them literal.
        self.assertEqual(args[args.index("--default-resources") + 1:],
                         ["slurm_account='PDS0325'", "slurm_partition='cpu'"])
        self.assertEqual(run.call_args.kwargs["cores"], 2)
        self.assertEqual(run.call_args.kwargs["config_dir"], Path("out") / ".kolach")

    def test_download_slurm_arguments(self):
        with mock.patch.object(cli, "run_workflow", return_value=0) as run:
            cli.main(["download", "--database-dir", "db", "--databases", "kofam", "eggnog",
                      "--executor", "slurm", "--slurm-account", "A"])
        self.assertIn("slurm_account='A'", run.call_args.kwargs["slurm_args"])
        self.assertEqual(run.call_args.kwargs["config_dir"], Path("db") / ".kolach")
        self.assertEqual(run.call_args.kwargs["cores"], 1)

    def test_user_profile_replaces_bundled(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "config.yaml").write_text("set-threads: {}\n")
            run, _ = self.run_main(self.BASE + ["--executor", "slurm", "--slurm-account", "A", "--slurm-profile", tmp])
        args = run.call_args.kwargs["slurm_args"]
        self.assertEqual(args[args.index("--workflow-profile") + 1], str(Path(tmp).resolve()))

    def test_deepkoala_gpu_only_with_cuda_on_slurm(self):
        slurm = ["--executor", "slurm", "--slurm-account", "A", "--slurm-partition", "cpu"]
        cases = [
            (["--databases", "deepkoala", "--device", "cuda", *slurm, "--slurm-gpu-partition", "gpu"],
             {"gpu": 1, "slurm_partition": "gpu"}),
            (["--databases", "deepkoala", "--device", "cuda", *slurm], {"gpu": 1, "slurm_partition": "cpu"}),
            (["--databases", "deepkoala", "--device", "auto", *slurm], None),
            (["--databases", "kofam", "--device", "cuda", *slurm], None),
            (["--databases", "deepkoala", "--device", "cuda"], None),
        ]
        for extra, expected in cases:
            with self.subTest(extra=extra):
                config = cli.annotate_config(parse(*self.BASE, *extra))
                self.assertEqual(config.get("deepkoala_resources"), expected)

    def test_slurm_values_stay_literal_strings(self):
        for value in ("all", "None", "0123", "1e3", "a'b"):
            args = parse(*self.BASE, "--executor", "slurm", "--slurm-account", value, "--slurm-partition", value)
            resources = dict(item.split("=", 1) for item in cli.slurm_arguments(args)[7:])
            self.assertEqual(eval(resources["slurm_account"]), value)
            self.assertEqual(eval(resources["slurm_partition"]), value)

    def test_user_profile_default_resources_are_kept(self):
        for body in ("default-resources:\n  mem_mb: 5000\n  slurm_extra: \"'--mail-type=NONE'\"\n  slurm_account: x\n",
                     "default-resources:\n  - mem_mb=5000\n  - \"slurm_extra='--mail-type=NONE'\"\n"):
            with self.subTest(body=body), tempfile.TemporaryDirectory() as tmp:
                (Path(tmp) / "config.yaml").write_text(body)
                args = parse(*self.BASE, "--executor", "slurm", "--slurm-account", "A", "--slurm-profile", tmp)
                slurm_args = cli.slurm_arguments(args)
                resources = slurm_args[slurm_args.index("--default-resources") + 1:]
                # The profile's defaults survive; --slurm-account wins over its account.
                self.assertEqual(resources, ["mem_mb=5000", "slurm_extra='--mail-type=NONE'", "slurm_account='A'"])

    @unittest.skipUnless(importlib.util.find_spec("snakemake_executor_plugin_slurm"), "Slurm plugin not installed")
    def test_profile_slurm_extra_with_managed_option_is_rejected(self):
        # The plugin only rejects these at submission, and the controller then hangs.
        bodies = [
            "default-resources:\n  slurm_extra: \"'--comment=x'\"\n",
            "set-resources:\n  annotate_kofam:\n    slurm_extra: \"'--qos=debug'\"\n",
            "set-resources:\n  - \"annotate_kofam:slurm_extra='--mem=4G'\"\n",
        ]
        for body in bodies:
            with self.subTest(body=body), tempfile.TemporaryDirectory() as tmp:
                (Path(tmp) / "config.yaml").write_text(body)
                self.assert_main_rejected(self.BASE + ["--executor", "slurm", "--slurm-account", "A",
                                                       "--slurm-profile", tmp], "which the Slurm executor manages")
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "config.yaml").write_text("default-resources:\n  slurm_extra: \"'--mail-type=NONE'\"\n")
            run, _ = self.run_main(self.BASE + ["--executor", "slurm", "--slurm-account", "A", "--slurm-profile", tmp])
            self.assertIn("slurm_extra='--mail-type=NONE'", run.call_args.kwargs["slurm_args"])

    def test_malformed_profile_is_rejected(self):
        for body in ("- just\n- a list\n", "default-resources:\n  - mem_mb\n", "a: [unclosed\n"):
            with self.subTest(body=body), tempfile.TemporaryDirectory() as tmp:
                (Path(tmp) / "config.yaml").write_text(body)
                self.assert_main_rejected(self.BASE + ["--executor", "slurm", "--slurm-account", "A",
                                                       "--slurm-profile", tmp], "invalid Slurm profile")

    def test_bundled_profile_has_no_default_resources(self):
        self.assertEqual(cli.profile_default_resources(cli.SLURM_PROFILE), {})

    def test_jobs_and_gpu_partition_need_their_context(self):
        self.assert_main_rejected(self.BASE + ["--jobs", "4"], "require(s) --executor slurm")
        slurm = ["--executor", "slurm", "--slurm-account", "A", "--slurm-gpu-partition", "gpu"]
        self.assert_main_rejected(self.BASE + ["--databases", "deepkoala", *slurm], "requires --device cuda")
        self.assert_main_rejected(self.BASE + ["--databases", "kofam", "--device", "cuda", *slurm],
                                  "requires --device cuda")

    def test_download_slurm_resolves_deepkoala_release_once(self):
        with mock.patch.object(cli, "run_workflow", return_value=0) as run, \
                mock.patch("kolach.methods.deepkoala.resolve_release", return_value="202608") as resolve:
            cli.main(["download", "--database-dir", "db", "--databases", "deepkoala",
                      "--executor", "slurm", "--slurm-account", "A"])
        resolve.assert_called_once_with("latest")
        config = run.call_args.args[1]
        self.assertEqual((config["deepkoala_release"], config["deepkoala_release_resolved"]), ("202608", True))
        with mock.patch.object(cli, "run_workflow", return_value=0) as run, \
                mock.patch("kolach.methods.deepkoala.resolve_release") as resolve:
            cli.main(["download", "--database-dir", "db", "--databases", "deepkoala"])
        resolve.assert_not_called()
        self.assertNotIn("deepkoala_release_resolved", run.call_args.args[1])

    def test_run_workflow_slurm_uses_unique_shared_configfile(self):
        seen = {}

        def fake_controller(cmd):
            path = Path(cmd[cmd.index("--configfile") + 1])
            seen.update(cmd=cmd, path=path, config=json.loads(path.read_text()))
            return 0

        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch("shutil.which", return_value="/bin/snakemake"), \
                mock.patch.object(cli, "run_controller", fake_controller):
            config_dir = Path(tmp) / "out" / ".kolach"
            slurm_args = ["--executor", "slurm", "--jobs", "2", "--default-resources", "slurm_account='A'"]
            code = cli.run_workflow("annotate.smk", {"threads": 2}, cores=3, slurm_args=slurm_args,
                                    config_dir=config_dir)
            first = seen["path"]
            cli.run_workflow("annotate.smk", {"threads": 2}, cores=3, slurm_args=slurm_args, config_dir=config_dir)
            self.assertEqual(code, 0)
            self.assertEqual(first.parent, config_dir.resolve())
            self.assertTrue(first.name.startswith("annotate_config_"))
            self.assertNotEqual(first, seen["path"])  # per-run name: a second run cannot overwrite it
            self.assertEqual(seen["config"], {"threads": 2})
            self.assertEqual(list(config_dir.iterdir()), [])  # removed once Snakemake returns
        cmd = seen["cmd"]
        self.assertNotIn("--cores", cmd)
        self.assertEqual(cmd[cmd.index("--local-cores") + 1], "3")
        # The configfile must follow the nargs='+' --default-resources values as its own flag.
        self.assertLess(cmd.index("slurm_account='A'"), cmd.index("--configfile"))

    def test_controller_turns_sigterm_into_sigint(self):
        # Snakemake cancels its Slurm jobs and unlocks on SIGINT but not on SIGTERM.
        child = ("import signal, sys, time\n"
                 "signal.signal(signal.SIGINT, lambda *a: sys.exit(42))\n"
                 "print('ready', flush=True)\n"
                 "time.sleep(30)\n")
        script = ("import os, signal, subprocess, sys, threading, time\n"
                  "from kolach import cli\n"
                  "threading.Timer(2, lambda: os.kill(os.getpid(), signal.SIGTERM)).start()\n"
                  f"sys.exit(cli.run_controller([sys.executable, '-c', {child!r}]))\n")
        env = dict(os.environ, PYTHONPATH=str(Path(cli.__file__).parents[1]))
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env, timeout=60)
        self.assertEqual(result.returncode, 42, result.stderr)


if __name__ == "__main__":
    unittest.main()
