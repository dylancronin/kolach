"""Tests for the pandas-based annotation integration module in kolach."""
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from kolach.integrate import (
    parse_kos,
    format_kos,
    read_fasta_ids,
    load_ko_definitions,
    load_kofam,
    load_deepkoala,
    load_eggnog,
    extract_kofam_records,
    extract_deepkoala_records,
    extract_eggnog_records,
    filter_and_disambiguate_eggnog,
    adjudicate_consensus,
    integrate_annotations,
    read_eggnog_tsv,
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
        self.assertEqual(row1["kofam_score_type"], "full")
        self.assertEqual(row1["kofam_bit_score"], 250.5)
        self.assertTrue(math.isclose(row1["kofam_evalue"], 1e-50, rel_tol=1e-12, abs_tol=0.0))
        self.assertEqual(row1["kofam_assignment"], "threshold")
        self.assertEqual(row1["kofam_threshold"], 100.0)

        # Multi-hit aggregation for gene2 preserves separate scores and statuses
        row2 = df[df["gene_id"] == "gene2"].iloc[0]
        self.assertEqual(row2["kofam_ko"], "K00002,K00003")
        self.assertEqual(row2["kofam_score_type"], "domain,domain")
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

        # gene3 has K05565 but probability 0.6004 < threshold 0.8674 -> below threshold candidate
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
            "kofam_assignment": ["threshold", "-", "threshold"],
            "deepkoala_ko": ["-", "-", "-"],
        })

        filtered = filter_and_disambiguate_eggnog(
            eggnog_df, base, min_bitscore=60.0, max_evalue=1e-5, filter_multi_mode="disambiguate"
        )

        row1 = filtered[filtered["gene_id"] == "gene1"].iloc[0]
        self.assertEqual(row1["eggnog_ko"], "K00001")
        self.assertEqual(row1["eggnog_bit_score"], 150.0)
        self.assertTrue(math.isclose(row1["eggnog_evalue"], 1e-40, rel_tol=1e-12, abs_tol=0.0))

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
            "kofam_assignment": ["threshold", "threshold", "threshold", "threshold", "-", "-", "threshold"],
            "deepkoala_ko": ["K00001", "K00002", "-", "K00005", "-", "-", "-"],
            "eggnog_ko": ["K00001", "K00099", "-", "K00006", "-", "K01447,K01448", "K07979"],
            "eggnog_candidate_ko": ["K00001", "K00099", "-", "K00006", "-", "K01447,K01448", "K00375,K07979"],
            "kofam_definition": ["Def 1", "Def 2", "Def 3", "Def 4", "-", "-", "GntR regulator"],
            "eggnog_description": ["-", "-", "-", "-", "-", "amidase", "secondary metabolite"],
        })

        active_tools = ["kofam", "deepkoala", "eggnog"]

        # 1. Multiple strategy (default)
        res_m = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="multiple")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "consensus_level"].values[0], "unanimous")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "accepted_ko"].values[0], "K00001")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unanimous", "alternative_kos"].values[0], "-")

        # g_majority has eggNOG calling K00099, so K00099 is preserved in alternative_kos
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "consensus_level"].values[0], "majority")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "accepted_ko"].values[0], "K00002")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_majority", "alternative_kos"].values[0], "K00099")

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "accepted_ko"].values[0], "K00003")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_single", "alternative_kos"].values[0], "-")

        # Solitary multi-KO hit: accepted_ko is strictly single-hit only ('-'), candidates in alternative_kos
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "accepted_ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "alternative_kos"].values[0], "K01447,K01448")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_multi_eggnog", "definition"].values[0], "-")

        # Disambiguated eggNOG multi-hit: accepted_ko has the agreed single KO ('K07979'), dropped eggNOG hit ('K00375') in alternative_kos
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "consensus_level"].values[0], "majority")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "accepted_ko"].values[0], "K07979")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "alternative_kos"].values[0], "K00375")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_disambiguated", "definition"].values[0], "GntR regulator")

        # Conflict under default multiple: accepted_ko is '-', alternative_kos has sorted KOs
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "accepted_ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "alternative_kos"].values[0], "K00004,K00005,K00006")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_conflict", "definition"].values[0], "-")
        self.assertIn("Def 4", res_m.loc[res_m["gene_id"] == "g_conflict", "alternative_definition"].values[0])

        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "consensus_level"].values[0], "unannotated")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "accepted_ko"].values[0], "-")
        self.assertEqual(res_m.loc[res_m["gene_id"] == "g_unannotated", "alternative_kos"].values[0], "-")

        # 2. Priority strategy
        res_p = adjudicate_consensus(data.copy(), active_tools, conflict_strategy="priority")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "consensus_level"].values[0], "conflict_priority")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "accepted_ko"].values[0], "K00004")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "alternative_kos"].values[0], "K00005,K00006")
        self.assertEqual(res_p.loc[res_p["gene_id"] == "g_conflict", "definition"].values[0], "Def 4")

    def test_priority_strategy_eggnog_over_deepkoala(self):
        """When KOfam is absent, priority strategy selects eggNOG over DeepKOALA."""
        data = pd.DataFrame({
            "gene_id": ["g1"],
            "kofam_ko": ["-"],
            "deepkoala_ko": ["K00002"],
            "eggnog_ko": ["K00003"],
            "eggnog_candidate_ko": ["K00003"],
        })
        active_tools = ["deepkoala", "eggnog"]
        res = adjudicate_consensus(data, active_tools, conflict_strategy="priority")
        self.assertEqual(res.loc[res["gene_id"] == "g1", "consensus_level"].values[0], "conflict_priority")
        self.assertEqual(res.loc[res["gene_id"] == "g1", "accepted_ko"].values[0], "K00003")
        self.assertEqual(res.loc[res["gene_id"] == "g1", "alternative_kos"].values[0], "K00002")

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
            self.assertEqual(en_k1.iloc[0]["cross_method_status"], "below_threshold_unrescued")

            # K00002 has NO eggNOG row
            en_k2 = ev_df[(ev_df["ko"] == "K00002") & (ev_df["method"] == "eggnog")]
            self.assertEqual(len(en_k2), 0)

    # -------------------------------------------------------------------------
    # Focused Tests for 4 Integration Requirements & Regression Cases
    # -------------------------------------------------------------------------

    def test_definition_resolution_winning_differs_from_kofam_regression(self):
        """Regression case: KOfam calls K00001 ('Function A'); DeepKOALA and eggNOG call K00002.

        Winning KO K00002 must never receive Function A merely because KOfam supplied that description.
        Minority call K00001 must be preserved in alternative_kos.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tfull\t100.0\t200.0\t1e-50\t200.0\t1e-50\tFunction A\n"
            )
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK00002\t0.990\t0.500\t*\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-50\t150.0\tK00002\teggnog seed description\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "K00002")
            self.assertEqual(row["alternative_kos"], "K00001")
            self.assertEqual(row["consensus_level"], "majority")

            # Definition of winning K00002 must NOT be Function A!
            self.assertNotEqual(row["definition"], "Function A")
            self.assertEqual(row["definition"], "-")

            # Alternative definition for K00001 must be Function A
            self.assertEqual(row["alternative_definition"], "Function A")

    def test_missing_ko_specific_definitions_return_dash(self):
        """When no matching KO definition is found in ko_list or evidence records, definition is '-'."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK99999\t0.950\t0.500\t*\n"
            )
            out_tsv = tmp / "out.tsv"

            df = integrate_annotations(
                deepkoala_tsv=dk,
                output_tsv=out_tsv,
            )

            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "K99999")
            self.assertEqual(row["definition"], "-")
            self.assertEqual(row["alternative_definition"], "-")

    def test_master_ko_list_definition_preferred_when_available(self):
        """Master ko_list definition takes precedence."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            db_dir = tmp / "databases"
            kofam_db = db_dir / "kofam"
            kofam_db.mkdir(parents=True)
            ko_list = kofam_db / "ko_list"
            ko_list.write_text("K00002\t100.0\tMaster Definition For K00002\n")

            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK00002\t0.950\t0.500\t*\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-50\t150.0\tK00002\teggnog seed description\n"
            )
            out_tsv = tmp / "out.tsv"

            df = integrate_annotations(
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                database_dir=db_dir,
                output_tsv=out_tsv,
            )

            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "K00002")
            self.assertEqual(row["definition"], "Master Definition For K00002")

    def test_missing_and_unknown_acceptance_status(self):
        """Missing or unrecognized status information must not default to threshold_passing.

        Unknown evidence must not count as a confident vote or silently enter candidate rescue.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            # DeepKOALA with missing score and threshold, empty annotate
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g_unknown\tK00001\t\t\t\n"
            )
            # KOfam with missing assignment and missing scores
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g_unknown\tK00001\t-\t-\t\t\t\t\t\t-\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            row = df[df["gene_id"] == "g_unknown"].iloc[0]
            self.assertEqual(row["accepted_ko"], "-")
            self.assertEqual(row["consensus_level"], "unannotated")
            self.assertEqual(row["evidence"], "-")

            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            self.assertEqual(len(ev_df), 2)
            for _, ev_row in ev_df.iterrows():
                self.assertEqual(ev_row["original_status"], "unknown")
                self.assertEqual(ev_row["cross_method_status"], "unknown")

    def test_mixed_complete_and_missing_score_rows(self):
        """DeepKOALA and eggNOG tables with mixed complete and missing score rows evaluate each row correctly."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK00001\t0.950\t0.500\t*\n"
                "g2\tK00002\t0.300\t0.500\t\n"
                "g3\tK00003\t\t\t\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-50\t150.0\tK00001\tdef1\n"
                "g2\ts2\t1e-02\t40.0\tK00002\tdef2\n"
                "g3\ts3\t\t\tK00003\tdef3\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)

            # g1: both methods threshold_passing -> unanimous
            g1_dk = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "deepkoala")].iloc[0]
            g1_en = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "eggnog")].iloc[0]
            self.assertEqual(g1_dk["original_status"], "threshold_passing")
            self.assertEqual(g1_en["original_status"], "threshold_passing")
            self.assertEqual(df[df["gene_id"] == "g1"].iloc[0]["consensus_level"], "unanimous")

            # g2: both below_threshold -> orthogonal_dual_candidate
            g2_dk = ev_df[(ev_df["gene_id"] == "g2") & (ev_df["method"] == "deepkoala")].iloc[0]
            g2_en = ev_df[(ev_df["gene_id"] == "g2") & (ev_df["method"] == "eggnog")].iloc[0]
            self.assertEqual(g2_dk["original_status"], "below_threshold")
            self.assertEqual(g2_en["original_status"], "below_threshold")
            self.assertEqual(df[df["gene_id"] == "g2"].iloc[0]["consensus_level"], "orthogonal_dual_candidate")

            # g3: both missing scores -> unknown -> unannotated
            g3_dk = ev_df[(ev_df["gene_id"] == "g3") & (ev_df["method"] == "deepkoala")].iloc[0]
            g3_en = ev_df[(ev_df["gene_id"] == "g3") & (ev_df["method"] == "eggnog")].iloc[0]
            self.assertEqual(g3_dk["original_status"], "unknown")
            self.assertEqual(g3_en["original_status"], "unknown")
            self.assertEqual(df[df["gene_id"] == "g3"].iloc[0]["consensus_level"], "unannotated")

    def test_minority_alternatives_retained_from_every_method(self):
        """Unselected threshold-passing calls from every method are preserved in alternative_kos."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)

            # Test A: DeepKOALA minority
            kf = tmp / "kf.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "gA\tK00001\tthreshold\tfull\t100.0\t200.0\t1e-50\t200.0\t1e-50\tdefA\n"
            )
            dk = tmp / "dk.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "gA\tK00003\t0.950\t0.500\t*\n"
            )
            en = tmp / "en.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "gA\ts1\t1e-50\t150.0\tK00001\tdefA\n"
            )
            out_tsv = tmp / "outA.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
            )
            rowA = df[df["gene_id"] == "gA"].iloc[0]
            self.assertEqual(rowA["accepted_ko"], "K00001")
            self.assertEqual(rowA["consensus_level"], "majority")
            # DeepKOALA minority K00003 preserved in alternative_kos
            self.assertEqual(rowA["alternative_kos"], "K00003")

    def test_one_resolved_accepted_ko_with_different_alternatives(self):
        """One resolved accepted KO with minority call and disambiguated dropped candidate."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kf.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tfull\t100.0\t200.0\t1e-50\t200.0\t1e-50\tdef1\n"
            )
            dk = tmp / "dk.tsv"
            dk.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g1\tK00002\t0.950\t0.500\t*\n"
            )
            en = tmp / "en.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-50\t150.0\tK00002,K00003\tdef_multi\n"
            )
            out_tsv = tmp / "out.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
            )
            row = df[df["gene_id"] == "g1"].iloc[0]
            # K00002 agreed by DeepKOALA and eggNOG -> majority
            self.assertEqual(row["accepted_ko"], "K00002")
            self.assertEqual(row["consensus_level"], "majority")
            # Both minority KOfam hit (K00001) and dropped eggNOG hit (K00003) in alternative_kos
            self.assertEqual(row["alternative_kos"], "K00001,K00003")

    def test_unresolved_multiple_candidates_produce_dash_accepted_ko(self):
        """Unresolved multi-KO calls produce accepted_ko = '-' and candidates in alternative_kos."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            # Solitary multi-KO eggNOG hit
            en = tmp / "en.tsv"
            en.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-50\t150.0\tK01447,K01448\tN-acetylmuramoyl-L-alanine amidase\n"
            )
            out_tsv = tmp / "out.tsv"

            df = integrate_annotations(
                eggnog_tsv=en,
                output_tsv=out_tsv,
            )
            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "-")
            self.assertEqual(row["alternative_kos"], "K01447,K01448")
            self.assertEqual(row["definition"], "-")
            self.assertEqual(row["alternative_definition"], "-")

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
            self.assertEqual(ev_df.iloc[0]["consensus_level"], "single_tool")
            self.assertEqual(ev_df.iloc[0]["original_status"], "threshold_passing")
            self.assertEqual(ev_df.iloc[0]["bit_score_threshold"], "248.0")
            self.assertEqual(ev_df.iloc[0]["e_value_threshold"], "-")
            self.assertEqual(ev_df.iloc[0]["eggnog_shared_seed_hit"], "False")
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

            # Check that new columns exist and redundant ko is absent
            self.assertIn("accepted_ko", df.columns)
            self.assertNotIn("ko", df.columns)
            self.assertIn("alternative_kos", df.columns)
            self.assertIn("alternative_definition", df.columns)
            self.assertIn("kofam_score_type", df.columns)
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
            self.assertEqual(g1["alternative_kos"], "-")
            self.assertEqual(g1["consensus_level"], "unanimous")
            self.assertEqual(g1["deepkoala_candidate_ko"], "K01783")
            self.assertEqual(g1["eggnog_candidate_ko"], "K01783")
            self.assertEqual(g1["kofam_score_type"], "full")
            self.assertEqual(g1["kofam_bit_score"], 316.8)
            self.assertTrue(math.isclose(g1["kofam_evalue"], 4.1e-95, rel_tol=1e-12, abs_tol=0.0))
            self.assertAlmostEqual(g1["deepkoala_score"], 0.991)
            self.assertEqual(g1["eggnog_bit_score"], 250.0)
            self.assertTrue(math.isclose(g1["eggnog_evalue"], 1e-80, rel_tol=1e-12, abs_tol=0.0))

            # Check g4 (multi-KO hit): accepted_ko must be strictly single-hit ('-'), candidates in alternative_kos
            g4 = df[df["gene_id"] == "g4"].iloc[0]
            self.assertEqual(g4["accepted_ko"], "-")
            self.assertEqual(g4["alternative_kos"], "K01447,K01448")
            self.assertEqual(g4["consensus_level"], "single_tool")
            self.assertEqual(g4["definition"], "-")
            # eggNOG generic description is NOT authoritative KO definition for multi-KO
            self.assertEqual(g4["alternative_definition"], "-")

            # Read exported TSV directly to check formatting
            tsv_read = pd.read_csv(out_tsv, sep="\t", dtype=str)
            self.assertEqual(len(tsv_read), 4)
            # Check g3 is unannotated
            g3_row = tsv_read[tsv_read["gene_id"] == "g3"].iloc[0]
            self.assertEqual(g3_row["accepted_ko"], "-")
            self.assertEqual(g3_row["alternative_kos"], "-")
            self.assertEqual(g3_row["consensus_level"], "unannotated")
            self.assertEqual(g3_row["kofam_bit_score"], "-")
            self.assertEqual(g3_row["deepkoala_candidate_ko"], "-")
            self.assertEqual(g3_row["eggnog_candidate_ko"], "-")

            # Evidence table also generated
            ev_file = tmp / "kolach_evidence.tsv"
            self.assertTrue(ev_file.exists())

    # =========================================================================
    # Targeted Regression Tests (Requirements 1-11)
    # =========================================================================

    def test_regression_1_full_and_domain_profiles_differ_substantially(self):
        """1. Full and domain profiles whose full/domain scores differ substantially."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g_full\tK00001\t-\tfull\t100.0\t250.0\t1e-50\t40.0\t1e-03\talcohol dehydrogenase\n"
                "g_domain\tK00002\t-\tdomain\t120.0\t350.0\t1e-80\t80.0\t1e-10\talcohol dehydrogenase (NADP+)\n"
            )
            records = extract_kofam_records(kf)
            self.assertEqual(len(records), 2)

            r_full = next(r for r in records if r["gene_id"] == "g_full")
            self.assertEqual(r_full["score_type"], "full")
            self.assertEqual(r_full["original_status"], "threshold_passing")
            self.assertEqual(r_full["bit_score"], 250.0)
            self.assertEqual(r_full["domain_bit_score"], 40.0)

            r_dom = next(r for r in records if r["gene_id"] == "g_domain")
            self.assertEqual(r_dom["score_type"], "domain")
            self.assertEqual(r_dom["original_status"], "below_threshold")
            self.assertEqual(r_dom["bit_score"], 350.0)
            self.assertEqual(r_dom["domain_bit_score"], 80.0)

    def test_regression_2_domain_profile_missing_domain_score_remains_unknown(self):
        """2. A domain profile with a missing domain score and a high full score remains unknown when status is inferred."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\t-\tdomain\t100.0\t450.0\t1e-90\t\t\tmissing domain score profile\n"
            )
            records = extract_kofam_records(kf)
            self.assertEqual(len(records), 1)
            r = records[0]
            self.assertEqual(r["original_status"], "unknown")
            self.assertEqual(r["bit_score"], 450.0)
            self.assertTrue(pd.isna(r["domain_bit_score"]))

    def test_regression_3_duplicate_domain_hits_selected_by_domain_score(self):
        """3. Duplicate domain hits being selected by domain score rather than full score."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            # Hit 1 has higher full score (400 > 200), but Hit 2 has higher domain score (180 > 90).
            # For domain profile, Hit 2 must be selected.
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tdomain\t100.0\t400.0\t1e-80\t90.0\t1e-15\thit1\n"
                "g1\tK00001\tthreshold\tdomain\t100.0\t200.0\t1e-40\t180.0\t1e-35\thit2\n"
            )
            records = extract_kofam_records(kf)
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0]["definition"], "hit2")
            self.assertEqual(records[0]["domain_bit_score"], 180.0)

    def test_regression_4_summary_scores_follow_score_type_evidence_unchanged(self):
        """4. Summary scores and E-values following score_type, while original evidence metrics remain unchanged."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g_dom\tK00001\tthreshold\tdomain\t100.0\t350.0\t1e-90\t150.5\t2.5e-30\tdomain protein\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            # Summary table asserts
            self.assertIn("kofam_score_type", df.columns)
            row = df[df["gene_id"] == "g_dom"].iloc[0]
            self.assertEqual(row["kofam_score_type"], "domain")
            self.assertEqual(row["kofam_bit_score"], 150.5)
            self.assertTrue(math.isclose(row["kofam_evalue"], 2.5e-30, rel_tol=1e-12, abs_tol=0.0))
            self.assertEqual(row["kofam_threshold"], 100.0)

            # Evidence table asserts: original full and domain metrics are unchanged
            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            ev_row = ev_df[ev_df["gene_id"] == "g_dom"].iloc[0]
            self.assertEqual(float(ev_row["bit_score"]), 350.0)
            self.assertTrue(math.isclose(float(ev_row["e_value"]), 1e-90, rel_tol=1e-12, abs_tol=0.0))
            self.assertEqual(float(ev_row["domain_bit_score"]), 150.5)
            self.assertTrue(math.isclose(float(ev_row["domain_e_value"]), 2.5e-30, rel_tol=1e-12, abs_tol=0.0))
            self.assertEqual(ev_row["score_type"], "domain")

    def test_regression_5_rescued_domain_hit_retains_threshold_and_label_displays_domain_score(self):
        """5. A rescued domain hit retaining its original KO threshold and rescued label while displaying its domain score."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00002\trescued\tdomain\t120.0\t180.0\t1e-25\t95.0\t5e-12\trescued domain KO\n"
            )
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"

            df = integrate_annotations(
                kofam_tsv=kf,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            row = df[df["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["accepted_ko"], "K00002")
            self.assertEqual(row["kofam_assignment"], "rescued")
            self.assertEqual(row["kofam_threshold"], 120.0)
            self.assertEqual(row["kofam_score_type"], "domain")
            self.assertEqual(row["kofam_bit_score"], 95.0)
            self.assertTrue(math.isclose(row["kofam_evalue"], 5e-12, rel_tol=1e-12, abs_tol=0.0))
            self.assertEqual(row["consensus_level"], "single_tool")
            self.assertEqual(row["evidence"], "kofam(rescued)")

            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)
            self.assertEqual(len(ev_df), 1)
            self.assertEqual(ev_df.iloc[0]["consensus_level"], "single_tool")
            self.assertEqual(ev_df.iloc[0]["original_status"], "heuristic_rescued")
            self.assertEqual(ev_df.iloc[0]["bit_score_threshold"], "90.0")  # 0.75 * 120.0
            self.assertEqual(float(ev_df.iloc[0]["e_value_threshold"]), 1e-5)
            self.assertEqual(ev_df.iloc[0]["eggnog_shared_seed_hit"], "False")

    def test_regression_6_native_and_normalized_eggnog_produce_identical_records(self):
        """6. Native and normalized eggNOG files producing identical records, retaining the first protein."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            native_file = tmp / "test.emapper.annotations"
            native_file.write_text(
                "## emapper-2.1.12\n"
                "## Command: emapper.py -i input.fa ...\n"
                "#query\tseed_ortholog\tevalue\tscore\tKEGG_ko\tDescription\n"
                "first_prot\ts1\t1e-50\t180.0\tko:K00001\tfirst description\n"
                "second_prot\ts2\t1e-30\t120.0\tko:K00002\tsecond description\n"
            )
            normalized_file = tmp / "normalized.tsv"
            normalized_file.write_text(
                "query\tseed_ortholog\tevalue\tscore\tKEGG_ko\tDescription\n"
                "first_prot\ts1\t1e-50\t180.0\tko:K00001\tfirst description\n"
                "second_prot\ts2\t1e-30\t120.0\tko:K00002\tsecond description\n"
            )

            rec_native = extract_eggnog_records(native_file)
            rec_norm = extract_eggnog_records(normalized_file)

            self.assertEqual(len(rec_native), 2)
            self.assertEqual(len(rec_norm), 2)
            self.assertEqual(rec_native[0]["gene_id"], "first_prot")
            self.assertEqual(rec_norm[0]["gene_id"], "first_prot")
            self.assertEqual(rec_native, rec_norm)

    def test_regression_7_header_only_and_malformed_eggnog(self):
        """7. Header-only eggNOG input and malformed headers."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)

            # Valid header-only native file
            header_only = tmp / "header_only.tsv"
            header_only.write_text(
                "## emapper metadata\n"
                "#query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
            )
            rec = extract_eggnog_records(header_only)
            self.assertEqual(rec, [])
            df_en = load_eggnog(header_only)
            self.assertTrue(df_en.empty)
            self.assertIn("gene_id", df_en.columns)

            # Malformed header: missing query
            missing_query = tmp / "missing_query.tsv"
            missing_query.write_text(
                "#seed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "s1\t1e-40\t150.0\tK00001\tdef\n"
            )
            with self.assertRaises(ValueError):
                extract_eggnog_records(missing_query)

            # Malformed header: missing KO
            missing_ko = tmp / "missing_ko.tsv"
            missing_ko.write_text(
                "#query\tseed\tevalue\tscore\tCOG\tDescription\n"
                "prot1\ts1\t1e-40\t150.0\tCOG0001\tdef\n"
            )
            with self.assertRaises(ValueError):
                extract_eggnog_records(missing_ko)

    def test_regression_8_eggnog_missing_metrics_unknown_and_zero_evalue_valid(self):
        """8. Missing/malformed/nonfinite eggNOG metrics remaining unknown, while E-value zero remains valid."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            en_file = tmp / "eggnog.tsv"
            en_file.write_text(
                "#query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g_zero\ts1\t0.0\t120.0\tko:K00001\tzero evalue\n"
                "g_missing_ev\ts2\t\t150.0\tko:K00002\tmissing evalue\n"
                "g_missing_bs\ts3\t1e-50\t\tko:K00003\tmissing bitscore\n"
                "g_inf\ts4\t1e-50\tinf\tko:K00004\tinf bitscore\n"
                "g_below\ts5\t1e-10\t45.0\tko:K00005\tbelow threshold\n"
            )

            records = extract_eggnog_records(en_file, min_bitscore=60.0, max_evalue=1e-5)
            rec_map = {r["gene_id"]: r for r in records}

            self.assertEqual(rec_map["g_zero"]["original_status"], "threshold_passing")
            self.assertEqual(rec_map["g_zero"]["e_value"], 0.0)

            self.assertEqual(rec_map["g_missing_ev"]["original_status"], "unknown")
            self.assertEqual(rec_map["g_missing_bs"]["original_status"], "unknown")
            self.assertEqual(rec_map["g_inf"]["original_status"], "unknown")
            self.assertEqual(rec_map["g_below"]["original_status"], "below_threshold")

            # Run integrate_annotations to confirm unknowns stay unknown without promotion
            out_tsv = tmp / "out.tsv"
            ev_tsv = tmp / "ev.tsv"
            df = integrate_annotations(eggnog_tsv=en_file, output_tsv=out_tsv, evidence_tsv=ev_tsv)

            self.assertEqual(df[df["gene_id"] == "g_zero"].iloc[0]["accepted_ko"], "K00001")
            self.assertEqual(df[df["gene_id"] == "g_missing_ev"].iloc[0]["consensus_level"], "unannotated")
            self.assertEqual(df[df["gene_id"] == "g_missing_bs"].iloc[0]["consensus_level"], "unannotated")
            self.assertEqual(df[df["gene_id"] == "g_inf"].iloc[0]["consensus_level"], "unannotated")

    def test_regression_9_explicit_missing_input_paths_raise_and_omitted_supported(self):
        """9. Explicit missing input paths raising errors and omitted optional paths remaining supported."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            out_tsv = tmp / "out.tsv"

            # 1. Nonexistent explicitly supplied KOfam file
            with self.assertRaises(FileNotFoundError):
                integrate_annotations(kofam_tsv=tmp / "nonexistent_kofam.tsv", output_tsv=out_tsv)

            # 2. Nonexistent explicitly supplied FASTA file
            with self.assertRaises(FileNotFoundError):
                integrate_annotations(protein_fasta=tmp / "nonexistent.faa", output_tsv=out_tsv)

            # 3. Completely omitted optional paths (None) work as empty runs
            df_empty = integrate_annotations(output_tsv=out_tsv)
            self.assertTrue(df_empty.empty)
            self.assertIn("accepted_ko", df_empty.columns)
            self.assertIn("kofam_score_type", df_empty.columns)

    def test_regression_10_legacy_adjudicate_consensus_single_tool_branch(self):
        """10. Legacy single-tool consensus branches completing without a NameError or column-length mismatch."""
        data = pd.DataFrame({
            "gene_id": ["g1", "g2"],
            "kofam_ko": ["K00001", "K00002"],
            "kofam_assignment": ["threshold", "rescued"],
            "deepkoala_ko": ["-", "-"],
            "eggnog_ko": ["-", "-"],
            "eggnog_candidate_ko": ["-", "-"],
        })
        active_tools = ["kofam", "deepkoala"]

        # Run adjudicate_consensus directly; must complete without NameError and length must match
        res = adjudicate_consensus(data, active_tools)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(res["consensus_level"]), 2)
        self.assertEqual(res.loc[res["gene_id"] == "g1", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res.loc[res["gene_id"] == "g2", "consensus_level"].values[0], "single_tool")
        self.assertEqual(res.loc[res["gene_id"] == "g2", "evidence"].values[0], "kofam(rescued)")

    def test_regression_11_accepted_ko_single_valued_and_alternatives_preserved(self):
        """11. Accepted KOs remaining single-valued and unresolved alternatives being preserved."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)

            # Solitary multi-KO hit
            en_tsv = tmp / "eggnog.tsv"
            en_tsv.write_text(
                "query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g_multi\ts1\t1e-50\t150.0\tko:K01447,ko:K01448\tmulti KO amidase\n"
            )
            df = integrate_annotations(eggnog_tsv=en_tsv)
            row = df[df["gene_id"] == "g_multi"].iloc[0]
            self.assertEqual(row["accepted_ko"], "-")
            self.assertEqual(row["alternative_kos"], "K01447,K01448")

            # Conflicting calls: KOfam calls K00001, DeepKOALA calls K00002
            kf_tsv = tmp / "kf.tsv"
            kf_tsv.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g_conflict\tK00001\tthreshold\tfull\t100.0\t200.0\t1e-50\t200.0\t1e-50\tdef1\n"
            )
            dk_tsv = tmp / "dk.tsv"
            dk_tsv.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "g_conflict\tK00002\t0.95\t0.50\t*\n"
            )
            df_conf = integrate_annotations(kofam_tsv=kf_tsv, deepkoala_tsv=dk_tsv, conflict_strategy="multiple")
            row_c = df_conf[df_conf["gene_id"] == "g_conflict"].iloc[0]
            self.assertEqual(row_c["accepted_ko"], "-")
            self.assertEqual(row_c["alternative_kos"], "K00001,K00002")
            self.assertEqual(row_c["consensus_level"], "conflict")

    def test_regression_12_kofam_missing_assignment_not_promoted_to_threshold_passing(self):
        """12. Missing KOfam status/assignment must not imply threshold passage in legacy helpers."""
        # In adjudicate_consensus, kofam_ko without assignment column or with "-" assignment
        data = pd.DataFrame({
            "gene_id": ["g1"],
            "kofam_ko": ["K00001"],
            # No kofam_assignment column
            "deepkoala_ko": ["-"],
            "eggnog_ko": ["-"],
            "eggnog_candidate_ko": ["-"],
        })
        active_tools = ["kofam"]
        res = adjudicate_consensus(data, active_tools)
        # Without kofam_assignment, K00001 must NOT be promoted to accepted threshold-passing call
        self.assertEqual(res.loc[res["gene_id"] == "g1", "accepted_ko"].values[0], "-")
        self.assertEqual(res.loc[res["gene_id"] == "g1", "consensus_level"].values[0], "unannotated")

        # In filter_and_disambiguate_eggnog:
        # eggnog has a multi-KO hit K00010,K00020.
        # If merged_df has kofam_ko="K00010" but NO kofam_assignment, KOfam cannot disambiguate as a trusted call
        eggnog_df = pd.DataFrame({
            "gene_id": ["g1"],
            "eggnog_raw_ko": ["K00010,K00020"],
            "eggnog_bit_score": [100.0],
            "eggnog_evalue": [1e-20],
            "eggnog_description": ["desc"],
        })
        base_df = pd.DataFrame({
            "gene_id": ["g1"],
            "kofam_ko": ["K00010"],
            "deepkoala_ko": ["-"],
            "deepkoala_candidate_ko": ["-"],
        })
        filtered = filter_and_disambiguate_eggnog(eggnog_df, base_df, filter_multi_mode="strict")
        # In strict mode without trusted other call, cannot disambiguate -> "-"
        self.assertEqual(filtered.loc[filtered["gene_id"] == "g1", "eggnog_ko"].values[0], "-")

    def test_regression_13_eggnog_duplicate_inf_score_does_not_displace_valid_finite_hit(self):
        """13. eggNOG duplicate with inf score must not displace a valid finite threshold-passing hit."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            tsv = tmp / "eggnog.tsv"
            # Row 1: inf score (unknown status)
            # Row 2: valid finite score 150.0, evalue 1e-50 (threshold_passing)
            tsv.write_text(
                "#query\tseed\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-50\tinf\tko:K00001\tinf hit\n"
                "g1\ts2\t1e-50\t150.0\tko:K00001\tvalid hit\n"
            )
            records = extract_eggnog_records(tsv, min_bitscore=60.0, max_evalue=1e-5)
            self.assertEqual(len(records), 1)
            rec = records[0]
            self.assertEqual(rec["bit_score"], 150.0)
            self.assertEqual(rec["original_status"], "threshold_passing")
            self.assertEqual(rec["eggnog_description"], "valid hit")

            # Also verify in load_eggnog
            df_en = load_eggnog(tsv)
            row = df_en[df_en["gene_id"] == "g1"].iloc[0]
            self.assertEqual(row["eggnog_bit_score"], 150.0)
            self.assertEqual(row["eggnog_candidate_ko"], "K00001")

    def test_regression_14_eggnog_pre_header_single_hash_comments_skipped(self):
        """14. Leading single-# comments before the recognized #query header are skipped."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            tsv = tmp / "eggnog.tsv"
            tsv.write_text(
                "## emapper-2.1.12\n"
                "# Command: emapper.py -i input.fa ...\n"
                "# Generated on 2026-09-14\n"
                "# Author: bioinfo\n"
                "#query\tseed_ortholog\tevalue\tscore\tKEGG_ko\tDescription\n"
                "prot1\ts1\t1e-50\t180.0\tko:K00001\tfirst protein\n"
                "prot2\ts2\t1e-30\t120.0\tko:K00002\tsecond protein\n"
            )
            df = read_eggnog_tsv(tsv)
            self.assertEqual(len(df), 2)
            self.assertIn("query", df.columns)
            self.assertEqual(df.iloc[0]["query"], "prot1")

            records = extract_eggnog_records(tsv)
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0]["gene_id"], "prot1")
            self.assertEqual(records[0]["ko"], "K00001")

            # Also check normalized TSV with pre-header single-# comment
            norm_tsv = tmp / "norm_eggnog.tsv"
            norm_tsv.write_text(
                "# Leading comment line\n"
                "# Another comment\n"
                "query\tseed_ortholog\tevalue\tscore\tKEGG_ko\tDescription\n"
                "prot1\ts1\t1e-50\t180.0\tko:K00001\tfirst protein\n"
            )
            df_norm = read_eggnog_tsv(norm_tsv)
            self.assertEqual(len(df_norm), 1)
            self.assertEqual(df_norm.iloc[0]["query"], "prot1")

    def test_evidence_table_full_schema_and_values(self):
        """Verify kolach_evidence.tsv contains consensus_level, correct thresholds, and eggnog_shared_seed_hit."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            kf = tmp / "kofam.tsv"
            kf.write_text(
                "gene_id\tko\tassignment\tscore_type\tthreshold\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tdefinition\n"
                "g1\tK00001\tthreshold\tfull\t200.0\t250.0\t1e-60\t250.0\t1e-60\talcohol dehydrogenase\n"
                "g2\tK00002\trescued\tfull\t100.0\t80.0\t1e-6\t80.0\t1e-6\talcohol dehydrogenase (NADP+)\n"
            )
            dk = tmp / "deepkoala.tsv"
            dk.write_text(
                "gene_id\tko\tdeepkoala_score\tdeepkoala_threshold\tannotate\n"
                "g1\tK00001\t0.85\t0.50\t*\n"
            )
            en = tmp / "eggnog.tsv"
            en.write_text(
                "#query\tseed_ortholog\tevalue\tscore\tKEGG_ko\tDescription\n"
                "g1\ts1\t1e-40\t150.0\tko:K00001,ko:K00009\tADH complex\n"
                "g3\ts3\t1e-20\t120.0\tko:K00003\thomoserine dehydrogenase\n"
            )
            out_tsv = tmp / "kolach_annotations.tsv"
            ev_tsv = tmp / "kolach_evidence.tsv"

            integrate_annotations(
                kofam_tsv=kf,
                deepkoala_tsv=dk,
                eggnog_tsv=en,
                output_tsv=out_tsv,
                evidence_tsv=ev_tsv,
            )

            self.assertTrue(ev_tsv.exists())
            ev_df = pd.read_csv(ev_tsv, sep="\t", dtype=str)

            for col in EVIDENCE_COLUMNS:
                self.assertIn(col, ev_df.columns)

            # Check g1 KOfam primary hit
            row_g1_kf = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "kofam")].iloc[0]
            self.assertEqual(row_g1_kf["consensus_level"], "unanimous")
            self.assertEqual(row_g1_kf["original_status"], "threshold_passing")
            self.assertEqual(row_g1_kf["bit_score_threshold"], "200.0")
            self.assertEqual(row_g1_kf["e_value_threshold"], "-")
            self.assertEqual(row_g1_kf["eggnog_shared_seed_hit"], "False")
            self.assertEqual(row_g1_kf["cross_method_status"], "accepted")

            # Check g1 DeepKOALA hit
            row_g1_dk = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "deepkoala")].iloc[0]
            self.assertEqual(row_g1_dk["consensus_level"], "unanimous")
            self.assertEqual(row_g1_dk["original_status"], "threshold_passing")
            self.assertEqual(row_g1_dk["bit_score_threshold"], "-")
            self.assertEqual(row_g1_dk["e_value_threshold"], "-")
            self.assertEqual(row_g1_dk["deepkoala_probability"], "0.85")
            self.assertEqual(row_g1_dk["deepkoala_threshold"], "0.5")
            self.assertEqual(row_g1_dk["eggnog_shared_seed_hit"], "False")
            self.assertEqual(row_g1_dk["cross_method_status"], "accepted")

            # Check g1 eggNOG multi-KO retained hit
            row_g1_en1 = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "eggnog") & (ev_df["ko"] == "K00001")].iloc[0]
            self.assertEqual(row_g1_en1["consensus_level"], "unanimous")
            self.assertEqual(row_g1_en1["original_status"], "threshold_passing")
            self.assertEqual(row_g1_en1["bit_score_threshold"], "60.0")
            self.assertEqual(float(row_g1_en1["e_value_threshold"]), 1e-5)
            self.assertEqual(row_g1_en1["eggnog_shared_seed_hit"], "True")
            self.assertEqual(row_g1_en1["cross_method_status"], "disambiguated_retained")

            # Check g1 eggNOG multi-KO dropped hit
            row_g1_en2 = ev_df[(ev_df["gene_id"] == "g1") & (ev_df["method"] == "eggnog") & (ev_df["ko"] == "K00009")].iloc[0]
            self.assertEqual(row_g1_en2["consensus_level"], "unanimous")
            self.assertEqual(row_g1_en2["bit_score_threshold"], "60.0")
            self.assertEqual(float(row_g1_en2["e_value_threshold"]), 1e-5)
            self.assertEqual(row_g1_en2["eggnog_shared_seed_hit"], "True")
            self.assertEqual(row_g1_en2["cross_method_status"], "disambiguated_dropped")

            # Check g2 KOfam rescued hit
            row_g2_kf = ev_df[(ev_df["gene_id"] == "g2") & (ev_df["method"] == "kofam")].iloc[0]
            self.assertEqual(row_g2_kf["consensus_level"], "single_tool")
            self.assertEqual(row_g2_kf["original_status"], "heuristic_rescued")
            self.assertEqual(row_g2_kf["bit_score_threshold"], "75.0")  # 0.75 * 100.0
            self.assertEqual(float(row_g2_kf["e_value_threshold"]), 1e-5)
            self.assertEqual(row_g2_kf["eggnog_shared_seed_hit"], "False")
            self.assertEqual(row_g2_kf["cross_method_status"], "heuristic_rescued")

            # Check g3 eggNOG single-KO hit
            row_g3_en = ev_df[(ev_df["gene_id"] == "g3") & (ev_df["method"] == "eggnog")].iloc[0]
            self.assertEqual(row_g3_en["consensus_level"], "single_tool")
            self.assertEqual(row_g3_en["original_status"], "threshold_passing")
            self.assertEqual(row_g3_en["bit_score_threshold"], "60.0")
            self.assertEqual(float(row_g3_en["e_value_threshold"]), 1e-5)
            self.assertEqual(row_g3_en["eggnog_shared_seed_hit"], "False")
            self.assertEqual(row_g3_en["cross_method_status"], "accepted")


if __name__ == "__main__":
    unittest.main()

