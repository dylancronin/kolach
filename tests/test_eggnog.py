"""Tests for eggNOG integration into kolach."""
import argparse
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from kolach.methods.eggnog import annotate


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

            def fake_subprocess_run(cmd, env=None, capture_output=True, text=True):
                # Find output_dir argument
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
                eggnog_dbmem=True,
            )

            annotate(args, str(output_tsv))

            self.assertTrue(output_tsv.exists())
            with open(output_tsv) as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)

            # Check that ## comment lines were stripped and #query became query
            self.assertEqual(len(rows), 2)
            self.assertEqual(
                rows[0],
                ["query", "seed_ortholog", "evalue", "score", "KEGG_ko", "Description"]
            )
            self.assertEqual(
                rows[1],
                ["seq1", "1234.PROT1", "1e-25", "120.5", "ko:K00001,ko:K00002", "alcohol dehydrogenase"]
            )

            # Verify subprocess call args
            call_cmd = mock_run.call_args[0][0]
            self.assertIn("-m", call_cmd)
            self.assertIn("diamond", call_cmd)
            self.assertIn("--cpu", call_cmd)
            self.assertIn("4", call_cmd)
            self.assertIn("--sensmode", call_cmd)
            self.assertIn("sensitive", call_cmd)
            self.assertIn("--dbmem", call_cmd)

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


if __name__ == "__main__":
    unittest.main()
