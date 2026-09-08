"""Tests for DeepKOALA integration into kolach."""
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from kolach.methods.deepkoala import (
    get_available_releases,
    setup_resources_link,
    annotate,
)


class TestDeepKoalaIntegration(unittest.TestCase):

    @patch("kolach.methods.deepkoala.urlopen")
    def test_get_available_releases(self, mock_urlopen):
        mock_html = """
        <html>
        <a href="202502/">202502/</a>
        <a href="202510/">202510/</a>
        <a href="202608/">202608/</a>
        </html>
        """
        mock_resp = MagicMock()
        mock_resp.read.return_value = mock_html.encode("utf-8")
        mock_urlopen.return_value.__enter__.return_value = mock_resp

        releases = get_available_releases()
        self.assertEqual(releases, ["202502", "202510", "202608"])

    def test_setup_resources_link_graceful_without_pkg(self):
        with tempfile.TemporaryDirectory() as tmp:
            # When deepkoala is not installed yet or in mock, function should not crash
            setup_resources_link(tmp, "202608")

    def test_annotate_tsv_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            db_dir = tmp_path / "deepkoala"
            (db_dir / "202608").mkdir(parents=True)
            output_tsv = tmp_path / "output.tsv"
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">seq1\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\n")

            def fake_inference(input_path, output_path, **kwargs):
                with open(output_path, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["name", "predict_label"])
                    writer.writerow(["seq1", "K00001"])

            import argparse
            args = argparse.Namespace(
                database_dir=str(tmp_path),
                protein_fasta=str(fasta),
                deepkoala_release="202608",
                deepkoala_model="full",
                threads=1,
                device="cpu",
                batch_size=32,
                detail=False,
            )

            mock_infer_module = MagicMock()
            mock_infer_module.inference = fake_inference
            mock_pkg = MagicMock()
            mock_pkg.__file__ = "/mock/path/deepkoala/__init__.py"

            with patch.dict("sys.modules", {"deepkoala": mock_pkg, "deepkoala.infer": mock_infer_module}):
                annotate(args, str(output_tsv))

            self.assertTrue(output_tsv.exists())
            with open(output_tsv) as f:
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)
                self.assertEqual(rows[0], ["name", "predict_label"])
                self.assertEqual(rows[1], ["seq1", "K00001"])


if __name__ == "__main__":
    unittest.main()
