"""Tests for the pandas-based annotation integration module in kolach."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from kolach.integrate import (
    parse_kos,
    format_kos,
    read_fasta_ids,
    load_kofam,
    load_deepkoala,
    load_eggnog,
    filter_and_disambiguate_eggnog,
    adjudicate_consensus,
    integrate_annotations,
)


class TestIntegrate(unittest.TestCase):

    def test_parse_kos(self):
        self.assertEqual(parse_kos("-"), set())
        self.assertEqual(parse_kos(""), set())
        self.assertEqual(parse_kos(None), set())
        self.assertEqual(parse_kos("ko:K00001"), {"K00001"})
        self.assertEqual(parse_kos("ko:K00001,ko:K00002"), {"K00001", "K00002"})
        self.assertEqual(parse_kos("K00001;K00003+K00005"), {"K00001", "K00003", "K00005"})
        self.assertEqual(format_kos({"K00002", "K00001"}), "K00001,K00002")
        self.assertEqual(format_kos(set()), "-")

    def test_load_kofam(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tmp:
            tmp.write(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "gene1\tK00001\tthreshold\tfull\t100.0\t250.5\t1e-50\t250.5\t1e-50\talcohol dehydrogenase\n"
                "gene2\tK00002\trescued\tdomain\t80.0\t70.0\t1e-08\t70.0\t1e-08\talcohol dehydrogenase (NADP+)\n"
                "gene2\tK00003\tthreshold\tdomain\t90.0\t120.0\t1e-20\t120.0\t1e-20\thomoserine dehydrogenase\n"
            )
            tmp_path = tmp.name

        df = load_kofam(tmp_path)
        self.assertEqual(len(df), 2)
        row1 = df[df["gene_id"] == "gene1"].iloc[0]
        self.assertEqual(row1["kofam_ko"], "K00001")
        self.assertEqual(row1["kofam_bit_score"], 250.5)
        self.assertEqual(row1["kofam_evalue"], 1e-50)
        self.assertEqual(row1["kofam_assignment"], "threshold")
        self.assertEqual(row1["kofam_threshold"], 100.0)

        # Multi-hit aggregation for gene2: threshold preferred, bit_score max
        row2 = df[df["gene_id"] == "gene2"].iloc[0]
        self.assertEqual(row2["kofam_ko"], "K00002,K00003")
        self.assertEqual(row2["kofam_assignment"], "threshold")
        self.assertEqual(row2["kofam_bit_score"], 120.0)

        Path(tmp_path).unlink(missing_ok=True)

    def test_load_deepkoala(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tmp:
            tmp.write(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "gene1\tK00001\t0.985\t0.500\t*\n"
                "gene2\t\t0.120\t0.600\t\n"
                "gene3\tK05565\t0.6004\t0.8674\t\n"
            )
            tmp_path = tmp.name

        df = load_deepkoala(tmp_path)
        self.assertEqual(len(df), 3)
        row1 = df[df["gene_id"] == "gene1"].iloc[0]
        self.assertEqual(row1["deepkoala_ko"], "K00001")
        self.assertAlmostEqual(row1["deepkoala_score"], 0.985)
        self.assertAlmostEqual(row1["deepkoala_threshold"], 0.500)

        # gene2 has empty predict_label
        row2 = df[df["gene_id"] == "gene2"].iloc[0]
        self.assertEqual(row2["deepkoala_ko"], "-")

        # gene3 has K05565 but probability 0.6004 < threshold 0.8674 -> must not be trusted!
        row3 = df[df["gene_id"] == "gene3"].iloc[0]
        self.assertEqual(row3["deepkoala_ko"], "-")
        self.assertAlmostEqual(row3["deepkoala_score"], 0.6004)
        self.assertAlmostEqual(row3["deepkoala_threshold"], 0.8674)

        Path(tmp_path).unlink(missing_ok=True)

    def test_load_eggnog_and_filtering(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tmp:
            tmp.write(
                "#query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "gene1\t123.1\t1e-40\t150.0\tko:K00001\talcohol dehydrogenase\n"
                "gene2\t123.2\t1e-05\t45.0\tko:K00002\tlow score hit\n"
                "gene3\t123.3\t1e-60\t200.0\tko:K00010,ko:K00020\tmulti hit\n"
            )
            tmp_path = tmp.name

        eggnog_df = load_eggnog(tmp_path)
        self.assertEqual(len(eggnog_df), 3)

        # Merge onto test base containing KOfam calls
        base = pd.DataFrame({
            "gene_id": ["gene1", "gene2", "gene3"],
            "kofam_ko": ["K00001", "-", "K00010"],
            "deepkoala_ko": ["-", "-", "-"],
        })

        filtered = filter_and_disambiguate_eggnog(
            eggnog_df, base, min_bitscore=60.0, filter_multi_mode="disambiguate"
        )

        row1 = filtered[filtered["gene_id"] == "gene1"].iloc[0]
        self.assertEqual(row1["eggnog_ko"], "K00001")
        self.assertEqual(row1["eggnog_bit_score"], 150.0)
        self.assertEqual(row1["eggnog_evalue"], 1e-40)

        # gene2 has bitscore 45 < 60 -> filtered to '-'
        row2 = filtered[filtered["gene_id"] == "gene2"].iloc[0]
        self.assertEqual(row2["eggnog_ko"], "-")
        self.assertEqual(row2["eggnog_bit_score"], 45.0)  # score retained!

        # gene3 multi-KO (K00010, K00020) disambiguated by KOfam (K00010)
        row3 = filtered[filtered["gene_id"] == "gene3"].iloc[0]
        self.assertEqual(row3["eggnog_ko"], "K00010")

        Path(tmp_path).unlink(missing_ok=True)

    def test_adjudicate_consensus_categories(self):
        data = pd.DataFrame({
            "gene_id": ["g_unanimous", "g_majority", "g_single", "g_conflict", "g_unannotated"],
            "kofam_ko": ["K00001", "K00002", "K00003", "K00004", "-"],
            "deepkoala_ko": ["K00001", "K00002", "-", "K00005", "-"],
            "eggnog_ko": ["K00001", "K00099", "-", "K00006", "-"],
            "kofam_definition": ["Def 1", "Def 2", "Def 3", "Def 4", "-"],
            "eggnog_description": ["-", "-", "-", "-", "-"],
        })

        active_tools = ["kofam", "deepkoala", "eggnog"]

        # 1. Multiple strategy (default: report as conflict and list all KOs)
        res_m = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="multiple")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "consensus_level"].values[0], "unanimous")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "ko"].values[0], "K00001")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "evidence"].values[0], "deepkoala,eggnog,kofam")

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "consensus_level"].values[0], "majority")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "ko"].values[0], "K00002")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "evidence"].values[0], "deepkoala,kofam")

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "ko"].values[0], "K00003")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "evidence"].values[0], "kofam")

        # Conflict: zero overlap -> level is 'conflict', and ko lists all candidate KOs comma-separated
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "ko"].values[0], "K00004,K00005,K00006")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "evidence"].values[0], "deepkoala,eggnog,kofam")

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "consensus_level"].values[0], "unannotated")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "ko"].values[0], "-")

        # 2. Priority strategy
        res_p = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="priority")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict_priority")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "ko"].values[0], "K00004")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "evidence"].values[0], "kofam")

        # 3. Drop conflict strategy
        res_d = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="drop")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict_dropped")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "ko"].values[0], "-")

    def test_end_to_end_integrate(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            fasta = tmp / "test.faa"
            fasta.write_text(">g1\nMKW\n>g2\nMLK\n>g3\nMST\n")

            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK01783\tthreshold\tfull\t248.0\t316.8\t4.1e-95\t316.7\t4.5e-95\tribulose-phosphate 3-epimerase\n"
            )

            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK01783\t0.991\t0.500\t*\n"
                "g2\tK04771\t0.850\t0.450\t*\n"
            )

            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-80\t250.0\tK01783\tribulose-phosphate 3-epimerase\n"
                "g2\ts2\t1e-40\t120.0\tK04771\tprotein translocase\n"
            )

            out_tsv = tmp / "kolach_annotations.tsv"

            df = integrate_annotations(
                protein_fasta=fasta,
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                eggnog_min_bitscore=60.0,
            )

            self.assertTrue(out_tsv.exists())
            self.assertEqual(len(df), 3)

            # Check that scores, evalues, and candidate KOs are retained in the output dataframe
            self.assertIn("kofam_bit_score", df.columns)
            self.assertIn("kofam_evalue", df.columns)
            self.assertIn("deepkoala_ko", df.columns)
            self.assertIn("deepkoala_candidate_ko", df.columns)
            self.assertIn("deepkoala_score", df.columns)
            self.assertIn("eggnog_ko", df.columns)
            self.assertIn("eggnog_candidate_ko", df.columns)
            self.assertIn("eggnog_bit_score", df.columns)
            self.assertIn("eggnog_evalue", df.columns)

            g1 = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(g1["ko"], "K01783")
            self.assertEqual(g1["consensus_level"], "unanimous")
            self.assertEqual(g1["deepkoala_candidate_ko"], "K01783")
            self.assertEqual(g1["eggnog_candidate_ko"], "K01783")
            self.assertEqual(g1["kofam_bit_score"], 316.8)
            self.assertEqual(g1["kofam_evalue"], 4.1e-95)
            self.assertAlmostEqual(g1["deepkoala_score"], 0.991)
            self.assertEqual(g1["eggnog_bit_score"], 250.0)
            self.assertEqual(g1["eggnog_evalue"], 1e-80)

            # Read exported TSV directly to check formatting
            tsv_read = pd.read_csv(out_tsv, sep="\t", dtype=str)
            self.assertEqual(len(tsv_read), 3)
            # Check g3 is unannotated
            g3_row = tsv_read[tsv_read["gene_id"] == "g3"].iloc[0]
            self.assertEqual(g3_row["ko"], "-")
            self.assertEqual(g3_row["consensus_level"], "unannotated")
            self.assertEqual(g3_row["kofam_bit_score"], "-")
            self.assertEqual(g3_row["deepkoala_candidate_ko"], "-")
            self.assertEqual(g3_row["eggnog_candidate_ko"], "-")


if __name__ == "__main__":
    unittest.main()
