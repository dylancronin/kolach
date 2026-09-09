"""Tests for eggNOG integration into kolach."""
import argparse
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from kolach.methods.eggnog import annotate
import kolach.cli as cli


class TestEggNOGIntegration(unittest.TestCase):

    @patch("kolach.methods.eggnog.shutil.which", return_value=None)
    def test_missing_executable_raises_error(self, mock_which):
        args = argparse.Namespace(
            database_dir="/fake/db",
            protein_fasta="/fake/proteins.faa",
        )
        with self.assertRaises(FileNotFoundError) as ctx:
            annotate(args, "/fake/out.tsv")
        self.assertIn("emapper.py", str(ctx.exception))

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    def test_missing_database_dir_raises_error(self, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(
                database_dir=tmp,
                protein_fasta="/fake/proteins.faa",
            )
            with self.assertRaises(FileNotFoundError) as ctx:
                annotate(args, f"{tmp}/out.tsv")
            self.assertIn("No eggNOG database directory found", str(ctx.exception))

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    def test_empty_database_dir_raises_error(self, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            eggnog_dir = Path(tmp) / "eggnog"
            eggnog_dir.mkdir()
            args = argparse.Namespace(
                database_dir=tmp,
                protein_fasta="/fake/proteins.faa",
            )
            with self.assertRaises(FileNotFoundError) as ctx:
                annotate(args, f"{tmp}/out.tsv")
            self.assertIn("is empty", str(ctx.exception))

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_annotate_success(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")

            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\n")
            output_tsv = tmp_path / "eggnog_annotations.tsv"

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]

                annotations_file = out_dir / f"{prefix}.emapper.annotations"
                content = (
                    "## emapper-2.1.12\n"
                    "## Execution config: ...\n"
                    "#query\tseed_ortholog\tevalue\tscore\tKEGG_ko\tDescription\n"
                    "seq1\t1234.PROT1\t1e-25\t120.5\tko:K00001,ko:K00002\talcohol dehydrogenase\n"
                    "## Done\n"
                )
                annotations_file.write_text(content)
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                threads=4,
                eggnog_mode="diamond",
                eggnog_sensmode="sensitive",
            )

            annotate(args, str(output_tsv))

            self.assertTrue(output_tsv.exists())
            with open(output_tsv) as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)

            self.assertEqual(len(rows), 2)
            self.assertEqual(
                rows[0],
                ["query", "seed_ortholog", "evalue", "score", "KEGG_ko", "Description"]
            )
            self.assertEqual(
                rows[1],
                ["seq1", "1234.PROT1", "1e-25", "120.5", "ko:K00001,ko:K00002", "alcohol dehydrogenase"]
            )

            call_cmd = mock_run.call_args[0][0]
            self.assertIn("-m", call_cmd)
            self.assertIn("diamond", call_cmd)
            self.assertIn("--cpu", call_cmd)
            self.assertIn("4", call_cmd)
            self.assertIn("--dmnd_sensmode", call_cmd)
            self.assertIn("sensitive", call_cmd)
            self.assertIn("--temp_dir", call_cmd)
            self.assertIn("--output_dir", call_cmd)

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_annotate_failure_raises_runtime_error(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            mock_res = MagicMock()
            mock_res.returncode = 1
            mock_res.stdout = "output error"
            mock_res.stderr = "Fatal error in diamond"
            mock_run.return_value = mock_res

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
            )

            with self.assertRaises(RuntimeError) as ctx:
                annotate(args, f"{tmp}/out.tsv")
            self.assertIn("emapper.py failed", str(ctx.exception))

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_sensmode_unspecified_omits_flag(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\nseq1\n")
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_sensmode=None,
            )
            annotate(args, f"{tmp}/out.tsv")
            call_cmd = mock_run.call_args[0][0]
            self.assertNotIn("--dmnd_sensmode", call_cmd)

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_sensmode_default_forwarded(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\nseq1\n")
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_sensmode="default",
            )
            annotate(args, f"{tmp}/out.tsv")
            call_cmd = mock_run.call_args[0][0]
            self.assertIn("--dmnd_sensmode", call_cmd)
            idx = call_cmd.index("--dmnd_sensmode")
            self.assertEqual(call_cmd[idx + 1], "default")

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_sensmode_explicit_fast_forwarded(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\nseq1\n")
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_sensmode="fast",
            )
            annotate(args, f"{tmp}/out.tsv")
            call_cmd = mock_run.call_args[0][0]
            self.assertIn("--dmnd_sensmode", call_cmd)
            idx = call_cmd.index("--dmnd_sensmode")
            self.assertEqual(call_cmd[idx + 1], "fast")

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_mmseqs_mode_excludes_diamond_flags(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\nseq1\n")
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_mode="mmseqs",
                eggnog_sensmode="sensitive",
                eggnog_dmnd_block_size=2.0,
                eggnog_dmnd_index_chunks=4,
            )
            annotate(args, f"{tmp}/out.tsv")
            call_cmd = mock_run.call_args[0][0]
            self.assertNotIn("--dmnd_sensmode", call_cmd)
            self.assertNotIn("--dmnd_block_size", call_cmd)
            self.assertNotIn("--dmnd_index_chunks", call_cmd)

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_temp_dir_custom(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            custom_scratch = tmp_path / "my_scratch"
            custom_scratch.mkdir()

            observed_tmp_dirs = []

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                temp_dir_idx = cmd.index("--temp_dir") + 1
                temp_dir = Path(cmd[temp_dir_idx])
                observed_tmp_dirs.append((out_dir, temp_dir))

                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\nseq1\n")
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_temp_dir=str(custom_scratch),
            )
            annotate(args, f"{tmp}/out.tsv")

            self.assertEqual(len(observed_tmp_dirs), 1)
            out_dir, temp_dir = observed_tmp_dirs[0]
            self.assertEqual(out_dir, temp_dir)
            # Must be placed inside custom_scratch
            self.assertEqual(out_dir.parent.resolve(), custom_scratch.resolve())
            # After annotate completes, the temp directory should be cleaned up
            self.assertFalse(out_dir.exists())

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    @patch("kolach.methods.eggnog.subprocess.run")
    def test_diamond_memory_controls_forwarded(self, mock_run, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            def fake_subprocess_run(cmd, env=None, **kwargs):
                out_dir_idx = cmd.index("--output_dir") + 1
                out_dir = Path(cmd[out_dir_idx])
                prefix_idx = cmd.index("-o") + 1
                prefix = cmd[prefix_idx]
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\nseq1\n")
                mock_res = MagicMock()
                mock_res.returncode = 0
                return mock_res

            mock_run.side_effect = fake_subprocess_run

            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_dmnd_block_size=2.5,
                eggnog_dmnd_index_chunks=4,
            )
            annotate(args, f"{tmp}/out.tsv")
            call_cmd = mock_run.call_args[0][0]
            self.assertIn("--dmnd_block_size", call_cmd)
            self.assertEqual(call_cmd[call_cmd.index("--dmnd_block_size") + 1], "2.5")
            self.assertIn("--dmnd_index_chunks", call_cmd)
            self.assertEqual(call_cmd[call_cmd.index("--dmnd_index_chunks") + 1], "4")

    @patch("kolach.methods.eggnog.shutil.which", return_value="/fake/emapper.py")
    def test_diamond_memory_controls_validation(self, mock_which):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "eggnog"
            db_dir.mkdir(parents=True)
            (db_dir / "eggnog.db").write_text("dummy")
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVT\n")

            args_invalid_block = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_dmnd_block_size=-1.0,
            )
            with self.assertRaises(ValueError):
                annotate(args_invalid_block, f"{tmp}/out.tsv")

            args_invalid_chunks = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                eggnog_dmnd_index_chunks=0,
            )
            with self.assertRaises(ValueError):
                annotate(args_invalid_chunks, f"{tmp}/out.tsv")

    def test_cli_dbmem_removed_raises_error(self):
        test_args = [
            "annotate",
            "--protein-fasta", "dummy.faa",
            "--database-dir", "dbs",
            "--eggnog-dbmem",
        ]
        with patch("sys.argv", ["kolach"] + test_args):
            with self.assertRaises(SystemExit):
                cli.main()

    def test_cli_defaults_sensitivity_to_none(self):
        # Test positive validators directly
        self.assertEqual(cli.positive_float("2.5"), 2.5)
        with self.assertRaises(argparse.ArgumentTypeError):
            cli.positive_float("-1.0")
        with self.assertRaises(argparse.ArgumentTypeError):
            cli.positive_float("0")

        self.assertEqual(cli.positive_int("4"), 4)
        with self.assertRaises(argparse.ArgumentTypeError):
            cli.positive_int("-2")
        with self.assertRaises(argparse.ArgumentTypeError):
            cli.positive_int("0")

    @patch("kolach.cli.subprocess.run")
    def test_cli_config_propagation_unspecified_sensitivity(self, mock_run):
        mock_res = MagicMock()
        mock_res.returncode = 0
        mock_run.return_value = mock_res

        test_args = [
            "kolach", "annotate",
            "--protein-fasta", "dummy.faa",
            "--database-dir", "dbs",
            "--databases", "eggnog",
        ]
        with patch("sys.argv", test_args):
            cli.main()

        cmd = mock_run.call_args[0][0]
        # Check that eggnog_sensmode is NOT in cmd
        sensmode_args = [arg for arg in cmd if arg.startswith("eggnog_sensmode=")]
        self.assertEqual(len(sensmode_args), 0)
        # Check that eggnog_temp_dir, block size, index chunks are also omitted
        self.assertFalse(any(arg.startswith("eggnog_temp_dir=") for arg in cmd))
        self.assertFalse(any(arg.startswith("eggnog_dmnd_block_size=") for arg in cmd))
        self.assertFalse(any(arg.startswith("eggnog_dmnd_index_chunks=") for arg in cmd))

    @patch("kolach.cli.subprocess.run")
    def test_cli_config_propagation_explicit_options(self, mock_run):
        mock_res = MagicMock()
        mock_res.returncode = 0
        mock_run.return_value = mock_res

        test_args = [
            "kolach", "annotate",
            "--protein-fasta", "dummy.faa",
            "--database-dir", "dbs",
            "--databases", "eggnog",
            "--eggnog-sensmode", "default",
            "--eggnog-temp-dir", "/scratch/tmp",
            "--eggnog-dmnd-block-size", "2.0",
            "--eggnog-dmnd-index-chunks", "4",
        ]
        with patch("sys.argv", test_args):
            cli.main()

        cmd = mock_run.call_args[0][0]
        self.assertIn("eggnog_sensmode=default", cmd)
        self.assertIn("eggnog_temp_dir=/scratch/tmp", cmd)
        self.assertIn("eggnog_dmnd_block_size=2.0", cmd)
        self.assertIn("eggnog_dmnd_index_chunks=4", cmd)


if __name__ == "__main__":
    unittest.main()
