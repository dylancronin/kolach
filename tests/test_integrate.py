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
    EVIDENCE_COLUMNS,
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

        # Multi-hit aggregation for gene2 preserves separate scores and statuses
        row2 = df[df["gene_id"] == "gene2"].iloc[0]
        self.assertEqual(row2["kofam_ko"], "K00002,K00003")
        self.assertIn("K00002:rescued", row2["kofam_assignment"])
        self.assertIn("K00003:threshold", row2["kofam_assignment"])
        self.assertIn("K00002:70.0", str(row2["kofam_bit_score"]))
        self.assertIn("K00003:120.0", str(row2["kofam_bit_score"]))

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
        self.assertEqual(row3["deepkoala_candidate_ko"], "K05565")
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

        base = pd.DataFrame({
            "gene_id": ["gene1", "gene2", "gene3"],
            "kofam_ko": ["K00001", "-", "K00010"],
            "deepkoala_ko": ["-", "-", "-"],
        })

        filtered = filter_and_disambiguate_eggnog(
            eggnog_df, base, min_bitscore=60.0, max_evalue=1e-5, filter_multi_mode="disambiguate"
        )

        row1 = filtered[filtered["gene_id"] == "gene1"].iloc[0]
        self.assertEqual(row1["eggnog_ko"], "K00001")
        self.assertEqual(row1["eggnog_bit_score"], 150.0)
        self.assertEqual(row1["eggnog_evalue"], 1e-40)

        # gene2 has bitscore 45 < 60 -> filtered to '-' (never promoted)
        row2 = filtered[filtered["gene_id"] == "gene2"].iloc[0]
        self.assertEqual(row2["eggnog_ko"], "-")
        self.assertEqual(row2["eggnog_candidate_ko"], "K00002")
        self.assertEqual(row2["eggnog_bit_score"], 45.0)

        # gene3 multi-KO disambiguated by KOfam (K00010)
        row3 = filtered[filtered["gene_id"] == "gene3"].iloc[0]
        self.assertEqual(row3["eggnog_ko"], "K00010")
        self.assertEqual(row3["eggnog_candidate_ko"], "K00010,K00020")

        Path(tmp_path).unlink(missing_ok=True)

    def test_adjudicate_consensus_categories(self):
        data = pd.DataFrame({
            "gene_id": ["g_unanimous", "g_majority", "g_single", "g_conflict", "g_unannotated", "g_multi_eggnog", "g_disambiguated"],
            "kofam_ko": ["K00001", "K00002", "K00003", "K00004", "-", "-", "K07979"],
            "deepkoala_ko": ["K00001", "K00002", "-", "K00005", "-", "-", "-"],
            "eggnog_ko": ["K00001", "K00099", "-", "K00006", "-", "K01447,K01448", "K07979"],
            "eggnog_candidate_ko": ["K00001", "K00099", "-", "K00006", "-", "K01447,K01448", "K00375,K07979"],
            "kofam_definition": ["Def 1", "Def 2", "Def 3", "Def 4", "-", "-", "GntR regulator"],
            "eggnog_description": ["-", "-", "-", "-", "-", "amidase", "secondary metabolite"],
        })

        active_tools = ["kofam", "deepkoala", "eggnog"]

        # 1. Multiple strategy (default: accepted_ko and ko set to '-', alternative_kos lists conflicting KOs)
        res_m = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="multiple")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "consensus_level"].values[0], "unanimous")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "accepted_ko"].values[0], "K00001")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "ko"].values[0], "K00001")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "alternative_kos"].values[0], "-")

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "consensus_level"].values[0], "majority")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "accepted_ko"].values[0], "K00002")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "ko"].values[0], "K00002")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "alternative_kos"].values[0], "-")

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "accepted_ko"].values[0], "K00003")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "ko"].values[0], "K00003")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "alternative_kos"].values[0], "-")

        # Solitary multi-KO hit: accepted_ko is strictly single-hit only ('-'), candidates in alternative_kos
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "accepted_ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "alternative_kos"].values[0], "K01447,K01448")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "definition"].values[0], "-")

        # Disambiguated eggNOG multi-hit: accepted_ko has the agreed single KO ('K07979'), dropped eggNOG hit ('K00375') in alternative_kos
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "consensus_level"].values[0], "majority")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "accepted_ko"].values[0], "K07979")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "ko"].values[0], "K07979")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "alternative_kos"].values[0], "K00375")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "definition"].values[0], "GntR regulator")

        # Conflict under default multiple: accepted_ko is '-', alternative_kos has sorted KOs
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "accepted_ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "alternative_kos"].values[0], "K00004,K00005,K00006")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "definition"].values[0], "-")
        self.assertIn("Def 4", res_m.loc[res_m["gene_id"] == "g_conflict", "alternative_definition"].values[0])

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "consensus_level"].values[0], "unannotated")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "accepted_ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "alternative_kos"].values[0], "-")

        # 2. Priority strategy
        res_p = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="priority")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict_priority")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "accepted_ko"].values[0], "K00004")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "ko"].values[0], "K00004")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "alternative_kos"].values[0], "K00005,K00006")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "definition"].values[0], "Def 4")

        # 3. Drop conflict strategy
        res_d = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="drop")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict_dropped")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "accepted_ko"].values[0], "-")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "ko"].values[0], "-")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "alternative_kos"].values[0], "K00004,K00005,K00006")
        self.assertEqual(res_d.loc[res_d["gene_id"] == "g_conflict", "definition"].values[0], "-")

    def test_kofam_confident_plus_subthreshold_eggnog_does_not_become_majority_or_unanimous(self):
        """A confident KOfam call plus matching below-threshold eggNOG evidence must remain single_tool_with_candidate."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tfull\t100.0\t250.0\t1e-50\t250.0\t1e-50\talcohol dehydrogenase\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-03\t45.0\tko:K00001\talcohol dehydrogenase\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
                eggnog_min_bitscore=60.0,
            )

            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "K00001")
            self.assertEqual(row["ko"], "K00001")
            self.assertNotEqual(row["consensus_level"], "majority")
            self.assertNotEqual(row["consensus_level"], "unanimous")
            self.assertEqual(row["consensus_level"], "single_tool_with_candidate")
            self.assertIn("eggnog(candidate)", row["evidence"])

            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            en_ev = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "eggnog")].iloc[0]
            self.assertEqual(en_ev["original_status"], "below_threshold")
            self.assertEqual(en_ev["cross_method_status"], "rescued")

    def test_deepkoala_confident_plus_subthreshold_eggnog_does_not_become_majority_or_unanimous(self):
        """A confident DeepKOALA call plus matching below-threshold eggNOG evidence must remain single_tool_with_candidate."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK00001\t0.950\t0.500\t*\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-03\t45.0\tko:K00001\talcohol dehydrogenase\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
                eggnog_min_bitscore=60.0,
            )

            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "K00001")
            self.assertNotEqual(row["consensus_level"], "majority")
            self.assertNotEqual(row["consensus_level"], "unanimous")
            self.assertEqual(row["consensus_level"], "single_tool_with_candidate")
            self.assertIn("deepkoala", row["evidence"])
            self.assertIn("eggnog(candidate)", row["evidence"])

    def test_genuine_threshold_passing_agreement(self):
        """Genuine independently threshold-passing agreement receives unanimous or majority."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tfull\t100.0\t250.0\t1e-50\t250.0\t1e-50\tdef1\n"
            )
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK00001\t0.950\t0.500\t*\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-40\t150.0\tko:K00001\tdef1\n"
            )
            out_tsv = tmp / "out.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
            )
            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["consensus_level"], "unanimous")
            self.assertEqual(row["accepted_ko"], "K00001")
            self.assertEqual(row["ko"], "K00001")

    def test_gene_with_threshold_and_heuristic_kofam_retains_separate_scores_and_statuses(self):
        """A gene with one threshold-passing KOfam KO and one heuristic-rescued KO retains separate scores and statuses."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g2\tK00002\trescued\tdomain\t80.0\t70.0\t1e-08\t70.0\t1e-08\tdef2\n"
                "g2\tK00003\tthreshold\tdomain\t90.0\t120.0\t1e-20\t120.0\t1e-20\tdef3\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            integrate_annotations(
                kofam_tsv=kf,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            self.assertEqual(len(ev_df), 2)

            k2 = ev_df[ev_df["ko"] == "K00002"].iloc[0]
            self.assertEqual(k2["original_status"], "heuristic_rescued")
            self.assertEqual(float(k2["bit_score"]), 70.0)
            self.assertEqual(float(k2["threshold"]), 80.0)

            k3 = ev_df[ev_df["ko"] == "K00003"].iloc[0]
            self.assertEqual(k3["original_status"], "threshold_passing")
            self.assertEqual(float(k3["bit_score"]), 120.0)
            self.assertEqual(float(k3["threshold"]), 90.0)

    def test_subthreshold_candidate_partial_overlap_does_not_leak_rescue(self):
        """A sub-threshold candidate that overlaps only one of several KOs does not confer rescued support on others."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tfull\t100.0\t200.0\t1e-50\t200.0\t1e-50\tdef1\n"
                "g1\tK00002\tthreshold\tfull\t100.0\t190.0\t1e-45\t190.0\t1e-45\tdef2\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-03\t45.0\tko:K00001\tdef1\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            integrate_annotations(
                kofam_tsv=kf,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
                eggnog_min_bitscore=60.0,
            )

            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            # K00001 has eggNOG candidate row
            en_k1 = ev_df[(ev_df["ko"] == "K00001") & (ev_df["method"] == "eggnog")]
            self.assertEqual(len(en_k1), 1)
            self.assertEqual(en_k1.iloc[0]["cross_method_status"], "rescued")

            # K00002 has NO eggNOG row
            en_k2 = ev_df[(ev_df["ko"] == "K00002") & (ev_df["method"] == "eggnog")]
            self.assertEqual(len(en_k2), 0)

    def test_annotation_and_standalone_integration_produce_evidence_table(self):
        """Verify kolach_evidence.tsv is automatically created alongside annotations."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK01783\tthreshold\tfull\t248.0\t316.8\t4.1e-95\t316.7\t4.5e-95\tribulose-phosphate 3-epimerase\n"
            )
            out_tsv = tmp / "kolach_annotations.tsv"

            integrate_annotations(
                kofam_tsv=kf,
                output_tsv=out_tsv,
            )

            expected_ev_path = tmp / "kolach_evidence.tsv"
            self.assertTrue(out_tsv.exists())
            self.assertTrue(expected_ev_path.exists())

            ev_df = pd.read_csv(expected_ev_path, sep="\t", dtype=str)
            for col in EVIDENCE_COLUMNS:
                self.assertIn(col, ev_df.columns)
            self.assertEqual(len(ev_df), 1)
            self.assertEqual(ev_df.iloc[0]["gene_id"], "g1")
            self.assertEqual(ev_df.iloc[0]["ko"], "K01783")
            self.assertEqual(ev_df.iloc[0]["method"], "kofam")
            self.assertEqual(ev_df.iloc[0]["original_status"], "threshold_passing")
            self.assertEqual(ev_df.iloc[0]["cross_method_status"], "accepted")

    def test_missing_methods_and_empty_results_handled_consistently(self):
        """Empty tables or missing methods produce consistent schema."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                protein_fasta=None,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            self.assertTrue(out_tsv.exists())
            self.assertTrue(ev_tsv.exists())
            self.assertTrue(df.empty)

            ev_read = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            for col in EVIDENCE_COLUMNS:
                self.assertIn(col, ev_read.columns)
            self.assertEqual(len(ev_read), 0)

    def test_end_to_end_integrate(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            fasta = tmp / "test.faa"
            fasta.write_text(">g1\nMKW\n>g2\nMLK\n>g3\nMST\n>g4\nMAA\n")

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
                "g4\ts4\t1e-45\t180.0\tK01447,K01448\tN-acetylmuramoyl-L-alanine amidase\n"
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
            self.assertEqual(len(df), 4)

            # Check that new and legacy columns exist
            self.assertIn("accepted_ko", df.columns)
            self.assertIn("ko", df.columns)
            self.assertIn("alternative_kos", df.columns)
            self.assertIn("alternative_definition", df.columns)
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
            self.assertEqual(g1["accepted_ko"], "K01783")
            self.assertEqual(g1["ko"], "K01783")
            self.assertEqual(g1["alternative_kos"], "-")
            self.assertEqual(g1["consensus_level"], "unanimous")
            self.assertEqual(g1["deepkoala_candidate_ko"], "K01783")
            self.assertEqual(g1["eggnog_candidate_ko"], "K01783")
            self.assertEqual(g1["kofam_bit_score"], 316.8)
            self.assertEqual(g1["kofam_evalue"], 4.1e-95)
            self.assertAlmostEqual(g1["deepkoala_score"], 0.991)
            self.assertEqual(g1["eggnog_bit_score"], 250.0)
            self.assertEqual(g1["eggnog_evalue"], 1e-80)

            # Check g4 (multi-KO hit): accepted_ko must be strictly single-hit ('-'), candidates in alternative_kos
            g4 = df[df["gene_id"] == "g4"].iloc[0]
            self.assertEqual(g4["accepted_ko"], "-")
            self.assertEqual(g4["ko"], "-")
            self.assertEqual(g4["alternative_kos"], "K01447,K01448")
            self.assertEqual(g4["consensus_level"], "single_tool")
            self.assertEqual(g4["definition"], "-")
            self.assertIn("N-acetylmuramoyl-L-alanine amidase", g4["alternative_definition"])

            # Read exported TSV directly to check formatting
            tsv_read = pd.read_csv(out_tsv, sep="\t", dtype=str)
            self.assertEqual(len(tsv_read), 4)
            # Check g3 is unannotated
            g3_row = tsv_read[tsv_read["gene_id"] == "g3"].iloc[0]
            self.assertEqual(g3_row["accepted_ko"], "-")
            self.assertEqual(g3_row["ko"], "-")
            self.assertEqual(g3_row["alternative_kos"], "-")
            self.assertEqual(g3_row["consensus_level"], "unannotated")
            self.assertEqual(g3_row["kofam_bit_score"], "-")
            self.assertEqual(g3_row["deepkoala_candidate_ko"], "-")
            self.assertEqual(g3_row["eggnog_candidate_ko"], "-")

            # Evidence table also generated
            ev_file = tmp / "kolach_evidence.tsv"
            self.assertTrue(ev_file.exists())


if __name__ == "__main__":
    unittest.main()
