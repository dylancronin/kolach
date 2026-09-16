"""Unit tests for kolach.pathways KEGG hierarchy parsing and pathway annotation."""
import json
from pathlib import Path
import tempfile
import unittest
import pandas as pd

from kolach.pathways import parse_kegg_hierarchy, annotate_pathways


class TestPathways(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temp_dir.name)

        # Create a mock ko00001.json
        self.mock_hierarchy = {
            "name": "ko00001",
            "children": [
                {
                    "name": "09100 Metabolism",
                    "children": [
                        {
                            "name": "09101 Carbohydrate metabolism",
                            "children": [
                                {
                                    "name": "00010 Glycolysis / Gluconeogenesis [PATH:ko00010]",
                                    "children": [
                                        {"name": "K00844  HK, glk; hexokinase [EC:2.7.1.1]"},
                                        {"name": "K00001  E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]"},
                                    ],
                                }
                            ],
                        }
                    ],
                },
                {
                    "name": "09180 Brite Hierarchies",
                    "children": [
                        {
                            "name": "09183 Protein families: signaling and cellular processes",
                            "children": [
                                {
                                    "name": "02000 Transporters [BR:ko02000]",
                                    "children": [
                                        {"name": "K05559  mnhA; multicomponent K+:H+ antiporter subunit A"}
                                    ],
                                }
                            ],
                        }
                    ],
                },
            ],
        }

        self.json_path = self.tmp / "ko00001.json"
        with open(self.json_path, "w", encoding="utf-8") as f:
            json.dump(self.mock_hierarchy, f)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_parse_kegg_hierarchy(self):
        mapping = parse_kegg_hierarchy(self.json_path)

        self.assertIn("K00844", mapping)
        self.assertEqual(mapping["K00844"]["category"], "09100 Metabolism")
        self.assertEqual(mapping["K00844"]["subcategory"], "09101 Carbohydrate metabolism")
        self.assertEqual(mapping["K00844"]["pathway"], "00010 Glycolysis / Gluconeogenesis")

        self.assertIn("K05559", mapping)
        self.assertEqual(mapping["K05559"]["category"], "09180 Brite Hierarchies")
        self.assertEqual(mapping["K05559"]["subcategory"], "09183 Protein families: signaling and cellular processes")
        self.assertEqual(mapping["K05559"]["pathway"], "02000 Transporters")

    def test_annotate_pathways(self):
        tsv_content = (
            "gene_id\taccepted_ko\tdefinition\n"
            "gene_1\tK00844\thexokinase\n"
            "gene_2\tK05559\tantiporter\n"
            "gene_3\t-\t-\n"
            "gene_4\tK99999\tunknown\n"
        )
        in_tsv = self.tmp / "annotations.tsv"
        in_tsv.write_text(tsv_content, encoding="utf-8")

        out_tsv = self.tmp / "annotated.tsv"
        df = annotate_pathways(
            annotation_table=in_tsv,
            output_tsv=out_tsv,
            database_dir=self.tmp,
            ko_col="accepted_ko",
            verbose=False,
        )

        self.assertTrue(out_tsv.is_file())
        self.assertIn("kegg_brite_category", df.columns)
        self.assertIn("kegg_brite_subcategory", df.columns)
        self.assertIn("kegg_brite_pathway", df.columns)

        # Check column ordering: kegg_brite_category, kegg_brite_subcategory, kegg_brite_pathway right after definition
        cols = list(df.columns)
        def_idx = cols.index("definition")
        self.assertEqual(cols[def_idx + 1], "kegg_brite_category")
        self.assertEqual(cols[def_idx + 2], "kegg_brite_subcategory")
        self.assertEqual(cols[def_idx + 3], "kegg_brite_pathway")

        # Row 1 (K00844)
        row1 = df.loc[df["gene_id"] == "gene_1"].iloc[0]
        self.assertEqual(row1["kegg_brite_category"], "09100 Metabolism")
        self.assertEqual(row1["kegg_brite_subcategory"], "09101 Carbohydrate metabolism")
        self.assertEqual(row1["kegg_brite_pathway"], "00010 Glycolysis / Gluconeogenesis")

        # Row 3 (unannotated)
        row3 = df.loc[df["gene_id"] == "gene_3"].iloc[0]
        self.assertEqual(row3["kegg_brite_category"], "-")
        self.assertEqual(row3["kegg_brite_subcategory"], "-")
        self.assertEqual(row3["kegg_brite_pathway"], "-")

        # Row 4 (unknown KO not in hierarchy)
        row4 = df.loc[df["gene_id"] == "gene_4"].iloc[0]
        self.assertEqual(row4["kegg_brite_category"], "-")
        self.assertEqual(row4["kegg_brite_subcategory"], "-")
        self.assertEqual(row4["kegg_brite_pathway"], "-")

    def test_reannotate_existing_table(self):
        tsv_content = (
            "gene_id\taccepted_ko\tdefinition\tcategory\tsubcategory\tpathway\n"
            "gene_1\tK00844\thexokinase\told_cat\told_subcat\told_pathway\n"
        )
        in_tsv = self.tmp / "already_annotated.tsv"
        in_tsv.write_text(tsv_content, encoding="utf-8")

        out_tsv = self.tmp / "reannotated.tsv"
        df = annotate_pathways(
            annotation_table=in_tsv,
            output_tsv=out_tsv,
            database_dir=self.tmp,
            ko_col="accepted_ko",
            verbose=False,
        )

        row1 = df.loc[df["gene_id"] == "gene_1"].iloc[0]
        self.assertNotIn("category", df.columns)
        self.assertEqual(row1["kegg_brite_category"], "09100 Metabolism")
        self.assertEqual(row1["kegg_brite_subcategory"], "09101 Carbohydrate metabolism")
        self.assertEqual(row1["kegg_brite_pathway"], "00010 Glycolysis / Gluconeogenesis")

    def test_is_valid_hierarchy_json(self):
        from kolach.pathways import is_valid_hierarchy_json

        # Valid file
        self.assertTrue(is_valid_hierarchy_json(self.json_path))

        # Missing file
        self.assertFalse(is_valid_hierarchy_json(self.tmp / "nonexistent.json"))

        # Empty file
        empty = self.tmp / "empty.json"
        empty.write_text("", encoding="utf-8")
        self.assertFalse(is_valid_hierarchy_json(empty))

        # Corrupt / invalid json
        corrupt = self.tmp / "corrupt.json"
        corrupt.write_text("<html>403 Forbidden</html>", encoding="utf-8")
        self.assertFalse(is_valid_hierarchy_json(corrupt))

        # JSON without children
        no_children = self.tmp / "no_children.json"
        no_children.write_text('{"name": "test"}', encoding="utf-8")
        self.assertFalse(is_valid_hierarchy_json(no_children))

    def test_strict_ko_regex_validation(self):
        # Tokens like K1234, K123456, Kabcde should not match or get annotated
        tsv_content = (
            "gene_id\taccepted_ko\tdefinition\n"
            "gene_1\tK00844\thexokinase\n"
            "gene_2\tK0084\ttruncated_ko\n"
            "gene_3\tK008444\tlong_ko\n"
            "gene_4\tKABCDE\talphabetic\n"
        )
        in_tsv = self.tmp / "invalid_kos.tsv"
        in_tsv.write_text(tsv_content, encoding="utf-8")

        out_tsv = self.tmp / "invalid_annotated.tsv"
        df = annotate_pathways(
            annotation_table=in_tsv,
            output_tsv=out_tsv,
            database_dir=self.tmp,
            ko_col="accepted_ko",
            verbose=False,
        )

        row1 = df.loc[df["gene_id"] == "gene_1"].iloc[0]
        self.assertEqual(row1["kegg_brite_category"], "09100 Metabolism")

        for bad_id in ["gene_2", "gene_3", "gene_4"]:
            row = df.loc[df["gene_id"] == bad_id].iloc[0]
            self.assertEqual(row["kegg_brite_category"], "-")
            self.assertEqual(row["kegg_brite_subcategory"], "-")
            self.assertEqual(row["kegg_brite_pathway"], "-")


if __name__ == "__main__":
    unittest.main()
