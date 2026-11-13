"""dbCAN conversion, workflow, and KO-consensus regression checks."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from kolach.integrate import integrate_annotations
from kolach.methods import dbcan


def native_fixture(directory, empty=False):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    hmm_cols = ["HMM Name", "Target Name", "Target From", "Target To", "i-Evalue", "Coverage"]
    sub_cols = ["Subfam Name", "Target Name", "Target From", "Target To", "i-Evalue", "Coverage",
                "Subfam Composition", "Subfam EC", "Substrate"]
    dia_cols = ["Gene ID", "CAZy ID", "Gene Start", "Gene End", "E Value", "Bit Score"]
    overview_cols = ["Gene ID", "#ofTools", "Recommend Results"]
    tables = [
        ("dbCAN_hmm_results.tsv", hmm_cols,
         [["GH13.hmm", "001", "1", "100", "1e-30", "0.8"],
          ["GH5.hmm", "001", "110", "210", "1e-20", "0.7"]]),
        ("dbCANsub_hmm_results.tsv", sub_cols,
         [["GH13_e1", "001", "1", "100", "1e-30", "0.8", "100", "3.2.1.1:10", "starch"],
          ["GH5_e2", "001", "110", "210", "1e-20", "0.7", "20", "3.2.1.4:3", "cellulose"]]),
        ("diamond.out", dia_cols, [["001", "AA12345.1|GH13|CBM20", "1", "210", "1e-120", "400"]]),
        ("overview.tsv", overview_cols, [["001", "3", "GH13_e1|GH5_e2"]]),
    ]
    for name, columns, values in tables:
        dbcan.write_table(directory / name, columns, [] if empty else [dict(zip(columns, row)) for row in values])


class TestDbcan(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.native = self.root / "native"
        native_fixture(self.native)
        self.fasta = self.root / "proteins.faa"
        self.fasta.write_text(">001\n" + "A" * 220 + "\n>empty\nAAAA\n")

    def summary(self):
        rows, _ = dbcan.parse_outputs(self.native, ["001", "empty"])
        path = self.root / "dbcan.tsv"
        dbcan.write_table(path, ["gene_id", *dbcan.SUMMARY_COLUMNS], rows)
        return path

    def test_domains_and_correlated_pairs(self):
        rows, evidence = dbcan.parse_outputs(self.native, ["001", "empty"])
        self.assertEqual(len(evidence), 5)
        self.assertEqual(rows[0]["dbcan_id"], "GH13;GH5")
        self.assertEqual(rows[0]["dbcan_family_ec"], "GH13;EC:3.2.1.1|GH5;EC:3.2.1.4")
        self.assertNotIn("GH13;EC:3.2.1.4", rows[0]["dbcan_family_ec"])
        self.assertEqual(rows[0]["dbcan_cazy_id"], "CBM20;GH13")
        self.assertEqual(rows[0]["dbcan_recommended"], "GH13_e1|GH5_e2")
        self.assertEqual(rows[1]["dbcan_tools_count"], "0")
        self.assertEqual(rows[1]["dbcan_id"], "-")
        self.assertEqual(json.loads(evidence[0]["source_row"])["Target Name"], "001")

    def test_empty_native_tables_are_valid(self):
        native_fixture(self.native, empty=True)
        rows, evidence = dbcan.parse_outputs(self.native, ["001", "empty"])
        self.assertFalse(evidence)
        self.assertTrue(all(row["dbcan_tools_count"] == "0" for row in rows))

    def test_reversible_search_ids(self):
        mapping = dbcan.prepare_input(self.fasta, self.native, ["001", "empty"])
        self.assertEqual(mapping["kolach_dbcan_00000000"], "001")
        self.assertTrue((self.native / "input.faa").read_text().startswith(">kolach_dbcan_00000000\n"))
        for path in self.native.glob("*"):
            if path.name in {"input.faa", "gene_ids.tsv"}:
                continue
            path.write_text(path.read_text().replace("001", "kolach_dbcan_00000000"))
        rows, evidence = dbcan.parse_outputs(self.native, ["001", "empty"], mapping)
        self.assertEqual(rows[0]["gene_id"], "001")
        self.assertTrue(all(row["gene_id"] == "001" for row in evidence))
        with self.assertRaises(ValueError):
            dbcan.parse_outputs(self.native, ["001", "001"])

    def test_missing_and_malformed_results_fail(self):
        for text in ["bad\theader\n", "Gene ID\t#ofTools\tRecommend Results\n001\t3\n",
                     "Gene ID\t#ofTools\tRecommend Results\n001\tnan\tGH13\n",
                     "Gene ID\t#ofTools\tRecommend Results\n001\t2\tGH13\n"]:
            with self.subTest(text=text):
                (self.native / "overview.tsv").write_text(text)
                with self.assertRaises(ValueError):
                    dbcan.parse_outputs(self.native, ["001", "empty"])
        (self.native / "overview.tsv").unlink()
        with self.assertRaises(FileNotFoundError):
            dbcan.parse_outputs(self.native, ["001", "empty"])

    def test_unknown_gene_and_duplicate_overview_fail(self):
        with self.assertRaises(ValueError):
            dbcan.parse_outputs(self.native, ["other"])
        with open(self.native / "overview.tsv", "a") as handle:
            handle.write("001\t3\tGH13\n")
        with self.assertRaises(ValueError):
            dbcan.parse_outputs(self.native, ["001"])

    def test_numeric_validation(self):
        path = self.native / "dbCAN_hmm_results.tsv"
        text = path.read_text()
        for old, bad in [("0.8", "nan"), ("0.8", "1.5"), ("1e-30", "-1"), ("\t100\t", "\t0\t")]:
            with self.subTest(bad=bad):
                path.write_text(text.replace(old, bad))
                with self.assertRaises(ValueError):
                    dbcan.parse_outputs(self.native, ["001"])

    def test_integrate_preserves_ko_decisions_and_evidence(self):
        kofam = self.root / "kofam.tsv"
        kofam.write_text("gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\n"
                         "001\tK00001\trescued\t80\t1e-7\t100\tfull\n")
        baseline = integrate_annotations(protein_fasta=self.fasta, kofam_tsv=kofam)
        enriched = integrate_annotations(protein_fasta=self.fasta, kofam_tsv=kofam, dbcan_tsv=self.summary())
        pd.testing.assert_frame_equal(baseline, enriched[baseline.columns])
        pd.testing.assert_frame_equal(baseline.attrs["evidence"], enriched.attrs["evidence"])
        self.assertEqual(enriched.iloc[0]["gene_id"], "001")
        self.assertEqual(enriched.iloc[0]["dbcan_id"], "GH13;GH5")

    def test_dbcan_only_and_union_inventory(self):
        table = self.summary()
        result = integrate_annotations(dbcan_tsv=table)
        self.assertEqual(list(result.gene_id), ["001", "empty"])
        self.assertEqual(set(result.accepted_ko), {"-"})
        self.assertEqual(set(result.consensus_level), {"unannotated"})
        self.assertTrue(result.attrs["evidence"].empty)
        fasta = self.root / "other.faa"
        fasta.write_text(">other\nAAAA\n")
        with self.assertRaises(ValueError):
            integrate_annotations(protein_fasta=fasta, dbcan_tsv=table)
        kofam = self.root / "kofam.tsv"
        kofam.write_text("gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\n"
                         "other\tK00001\tthreshold\t110\t1e-20\t100\tfull\n")
        result = integrate_annotations(kofam_tsv=kofam, dbcan_tsv=table)
        self.assertEqual(list(result.gene_id), ["other", "001", "empty"])
        self.assertEqual(result.iloc[0]["dbcan_tools_count"], "0")

    def test_duplicate_summary_rejected(self):
        table = self.summary()
        with open(table, "a") as handle:
            handle.write(table.read_text().splitlines()[1] + "\n")
        with self.assertRaises(ValueError):
            integrate_annotations(dbcan_tsv=table)

    def test_annotate_command_and_provenance(self):
        output = self.root / "out" / "dbcan_annotations.tsv"
        calls = []

        def fake_run(command, log):
            calls.append(command)
            native_fixture(command[command.index("--output_dir") + 1])
            Path(log).touch()

        with patch.object(dbcan, "validate_database"), patch.object(dbcan, "tool_version", return_value="5.2.9"), \
                patch.object(dbcan, "sha256", return_value="checksum"), patch.object(dbcan, "_run", side_effect=fake_run), \
                patch.object(dbcan.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "diamond 2.2.6\n")):
            dbcan.annotate(argparse.Namespace(protein_fasta=self.fasta, database_dir=self.root, threads="4"), output)
            note = output.parent / "dbcan" / "user-notes.txt"
            note.write_text("preserve me")
            dbcan.annotate(argparse.Namespace(protein_fasta=self.fasta, database_dir=self.root, threads="4"), output)
            self.assertEqual(note.read_text(), "preserve me")
        command = calls[0]
        self.assertEqual(command[command.index("--threads") + 1], "4")
        self.assertEqual(command[command.index("--mode") + 1], "protein")
        self.assertEqual(command[command.index("--methods") + 1], "diamond,hmm,dbCANsub")
        self.assertTrue((output.parent / "dbcan" / "overview.tsv").exists())
        self.assertEqual(len(dbcan.read_summary(output)), 2)
        metadata = json.loads((output.parent / "dbcan_metadata.json").read_text())
        self.assertEqual(metadata["protein_fasta"]["sha256"], "checksum")

    def test_output_collision_guards_preserve_inputs(self):
        database = self.root / "dbcan"
        database.mkdir()
        for name in dbcan.ASSETS:
            (database / name).write_text("original database")
        cases = [self.root / "dbcan_annotations.tsv"]
        nested = database / "nested"
        nested.mkdir()
        cases.append(nested / "dbcan_annotations.tsv")
        for output in cases:
            with self.subTest(output=output), patch.object(dbcan, "validate_database") as validation:
                with self.assertRaisesRegex(ValueError, "overlap"):
                    dbcan.annotate(argparse.Namespace(protein_fasta=self.fasta, database_dir=self.root, threads=1), output)
                validation.assert_not_called()
        self.assertTrue(all((database / name).read_text() == "original database" for name in dbcan.ASSETS))
        out = self.root / "results"
        native = out / "dbcan"
        native.mkdir(parents=True)
        note = native / "keep.txt"
        note.write_text("user data")
        with self.assertRaisesRegex(ValueError, "nonempty"):
            dbcan.annotate(argparse.Namespace(protein_fasta=self.fasta, database_dir=self.root, threads=1), out / "dbcan_annotations.tsv")
        self.assertEqual(note.read_text(), "user data")
        other = self.root / "other"
        other.mkdir()
        (other / "dbcan").symlink_to(database, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symlink"):
            dbcan.annotate(argparse.Namespace(protein_fasta=self.fasta, database_dir=self.root, threads=1), other / "dbcan_annotations.tsv")
        # Guard input FASTAs even in a directory previously created by Kolach.
        (native / dbcan.OWNER_FILE).write_text(dbcan.OWNER_TEXT)
        fasta = native / "input.faa"
        fasta.write_text(">original\nAAAA\n")
        with self.assertRaisesRegex(ValueError, "input FASTA"):
            dbcan.annotate(argparse.Namespace(protein_fasta=fasta, database_dir=self.root, threads=1), out / "dbcan_annotations.tsv")
        self.assertEqual(fasta.read_text(), ">original\nAAAA\n")

    def test_sidecar_symlinks_and_hardlinks_cannot_overwrite_sources(self):
        database = self.root / "dbcan"
        database.mkdir()
        asset = database / "CAZy.dmnd"
        asset.write_text("original database")
        out = self.root / "results"
        out.mkdir()
        for filename in ("dbcan_evidence.tsv", "dbcan_metadata.json"):
            for hardlink in (False, True):
                with self.subTest(filename=filename, hardlink=hardlink):
                    link = out / filename
                    if hardlink:
                        os.link(asset, link)
                    else:
                        link.symlink_to(asset)
                    with self.assertRaisesRegex(ValueError, "alias"):
                        dbcan.annotate(argparse.Namespace(protein_fasta=self.fasta, database_dir=self.root, threads=1),
                                       out / "dbcan_annotations.tsv")
                    self.assertEqual(asset.read_text(), "original database")
                    link.unlink()

    def test_download_marker_only_after_validation(self):
        directory = self.root / "db"
        payload = "fixture"

        def fake_download(url, destination):
            self.assertEqual(url, f"{dbcan.BASE_URL}/{dbcan.ASSETS[destination.name]}")
            destination.write_text(payload)
            return dbcan.sha256(destination)

        with patch.object(dbcan, "download_atomic", side_effect=fake_download), \
                patch.object(dbcan, "validate_database"):
            dbcan.download_dbcan(directory)
        self.assertTrue((directory / ".download_complete").exists())
        manifest_text = (directory / "manifest.json").read_text()
        manifest = json.loads(manifest_text)
        self.assertEqual(manifest["release"], dbcan.RELEASE)
        self.assertEqual(manifest["compatible_tool_version"], dbcan.VERSION)
        for name in dbcan.ASSETS:
            self.assertEqual(manifest["files"][name]["sha256"], dbcan.sha256(directory / name))
        payload = "invalid replacement"
        with patch.object(dbcan, "download_atomic", side_effect=fake_download), \
                patch.object(dbcan, "validate_database", side_effect=ValueError("broken")):
            with self.assertRaises(ValueError):
                dbcan.download_dbcan(directory)
        self.assertFalse((directory / ".download_complete").exists())
        self.assertEqual((directory / "manifest.json").read_text(), manifest_text)
        self.assertTrue(all((directory / name).read_text() == "fixture" for name in dbcan.ASSETS))

    def test_interrupted_download_preserves_existing_bundle(self):
        from io import BytesIO
        from urllib.error import URLError

        directory = self.root / "db"
        directory.mkdir()
        for name in [*dbcan.ASSETS, "manifest.json", ".download_complete"]:
            (directory / name).write_text("previous")
        with patch("kolach.download.urlopen", side_effect=[BytesIO(b"new asset"), URLError("interrupted")]), \
                patch.object(dbcan, "validate_database") as validate:
            with self.assertRaises(URLError):
                dbcan.download_dbcan(directory)
        validate.assert_not_called()
        self.assertFalse((directory / ".download_complete").exists())
        self.assertTrue(all((directory / name).read_text() == "previous"
                            for name in [*dbcan.ASSETS, "manifest.json"]))
        self.assertEqual(sorted(p.name for p in directory.iterdir()), sorted([*dbcan.ASSETS, "manifest.json"]))

    def test_database_and_process_failures(self):
        with self.assertRaises(FileNotFoundError):
            dbcan.validate_database(self.root)
        with self.assertRaises(subprocess.CalledProcessError):
            dbcan._run([sys.executable, "-c", "raise SystemExit(1)"], self.root / "failed.log")
        with self.assertRaises(RuntimeError):
            dbcan._run([sys.executable, "-c", "print('ERROR: caught upstream error')"], self.root / "caught.log")
        with patch.object(dbcan.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "dbCAN version: 4.0\n")):
            with self.assertRaises(RuntimeError):
                dbcan.tool_version()

    def test_cli_forwarding(self):
        from kolach.cli import main
        for command in ["download", "annotate"]:
            argv = ["kolach", command, "--databases", "dbcan", "--database-dir", str(self.root)]
            if command == "annotate":
                argv += ["--protein-fasta", str(self.fasta)]
            with patch.object(sys, "argv", argv), patch("kolach.cli.run_workflow") as run, \
                    patch.object(dbcan, "check_dependencies"):
                run.return_value = 0
                main()
                self.assertEqual(run.call_args.args[0], f"{command}.smk")
                self.assertEqual(run.call_args.args[1]["databases"], ["dbcan"])
        output = self.root / "integrated.tsv"
        argv = ["kolach", "integrate", "--dbcan-table", str(self.summary()), "--output-file", str(output)]
        with patch.object(sys, "argv", argv):
            main()
        self.assertIn("dbcan_id", output.read_text().splitlines()[0])

    def test_cli_preflight_reports_missing_dbcan_executable(self):
        from kolach.cli import main
        import contextlib
        import io

        for missing in ("run_dbcan", "hmmstat", "diamond"):
            with self.subTest(missing=missing), \
                    patch.object(dbcan.shutil, "which", side_effect=lambda name: None if name == missing else f"/bin/{name}"), \
                    patch("kolach.cli.run_workflow") as run, \
                    contextlib.redirect_stderr(io.StringIO()) as error, \
                    self.assertRaises(SystemExit) as exit_:
                main(["annotate", "--databases", "dbcan", "--database-dir", str(self.root),
                      "--protein-fasta", str(self.fasta)])
            self.assertEqual(exit_.exception.code, 1)
            self.assertIn(missing, error.getvalue())
            run.assert_not_called()


@unittest.skipUnless(shutil.which("snakemake"), "Snakemake is optional")
class TestDbcanWorkflow(unittest.TestCase):
    """Run the real DAG with deterministic stand-ins for search executables."""

    setUp = TestDbcan.setUp
    def test_workflow_reruns_and_pathway(self):
        root = Path(__file__).resolve().parents[1]
        bindir = self.root / "bin"
        bindir.mkdir()
        script = """#!/usr/bin/env python3
import os, sys, pathlib, shutil
if pathlib.Path(sys.argv[0]).name != 'run_dbcan':
    raise SystemExit(0)
if sys.argv[1] == 'version':
    print('dbCAN version: 5.2.9')
else:
    out = pathlib.Path(sys.argv[sys.argv.index('--output_dir')+1])
    out.mkdir(parents=True, exist_ok=True)
    for source in pathlib.Path(os.environ['DBCAN_FIXTURE']).glob('*'):
        shutil.copyfile(source, out / source.name)
"""
        for name in ["run_dbcan", "hmmstat", "hmmsearch", "diamond"]:
            path = bindir / name
            path.write_text(script)
            path.chmod(0o755)
        database = self.root / "database"
        db = database / "dbcan"
        db.mkdir(parents=True)
        for name in dbcan.ASSETS:
            (db / name).write_text("fixture")
        (db / "fam-substrate-mapping.tsv").write_text(
            "Family\tSubstrate_curated\tEC_Number\nGH13\tstarch\t3.2.1.1\n")
        (database / "ko00001.json").write_text(json.dumps({"name": "root", "children": []}))
        out = self.root / "results"
        env = {**os.environ, "PATH": f"{bindir}:{os.environ['PATH']}",
               "DBCAN_FIXTURE": str(self.native), "PYTHONPATH": str(root / "src"),
               "XDG_CACHE_HOME": str(self.root / "cache")}
        cmd = [sys.executable, "-m", "kolach.cli", "annotate", "--databases", "dbcan",
               "--protein-fasta", str(self.fasta), "--database-dir", str(database),
               "--output-dir", str(out), "--add-pathways"]
        # cli.py is an executable module; Snakemake uses the same isolated cwd.
        result = subprocess.run(cmd, env=env, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        summary = out / "kolach_annotations.tsv"
        self.assertTrue(summary.exists())
        self.assertIn("dbcan_id", summary.read_text().splitlines()[0])
        timestamp = summary.stat().st_mtime_ns
        result = subprocess.run(cmd, env=env, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(summary.stat().st_mtime_ns, timestamp)
        (db / "dbCAN.hmm").write_text("changed fixture")
        result = subprocess.run(cmd, env=env, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertGreater(summary.stat().st_mtime_ns, timestamp)
        # Add a captured KOfam table; exercise the mixed-tool integration DAG
        # without running another annotation search.
        kofam_db = database / "kofam"
        kofam_db.mkdir()
        (kofam_db / "profiles.hmm").write_text("fixture")
        (kofam_db / "ko_list").write_text("knum\tdefinition\nK00001\tExample\n")
        (out / "kofam_annotations.tsv").write_text(
            "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\n"
            "001\tK00001\tthreshold\t110\t1e-20\t100\tfull\n")
        mixed = cmd.copy()
        mixed.insert(mixed.index("dbcan"), "kofam")
        result = subprocess.run(mixed, env=env, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with open(summary) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["accepted_ko"], "K00001")
        self.assertEqual(rows[0]["dbcan_id"], "GH13;GH5")
        # Every advertised deliverable must be repaired when removed, even
        # when the summary consumed by integration is still present.
        repair = [out / "dbcan_evidence.tsv", out / "dbcan_metadata.json", out / "dbcan" / "overview.tsv"]
        for path in repair:
            path.unlink()
        result = subprocess.run(mixed, env=env, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(all(path.exists() for path in repair))
        # A success marker does not hide a missing asset at DAG-build time.
        (db / ".download_complete").touch()
        (db / "dbCAN.hmm").unlink()
        fresh = cmd.copy()
        fresh[fresh.index("--output-dir") + 1] = str(self.root / "fresh")
        result = subprocess.run(fresh, env=env, cwd=self.root, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("MissingInputException", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
