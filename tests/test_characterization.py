"""Characterization test suite for kolach.integrate.

Covers consensus adjudication, multi-tool agreement/conflict, candidate rescue,
eggNOG multi-KO disambiguation, FASTA ordering, score selection, definition
resolution, and table exports.
"""
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from kolach.integrate import (
    integrate_annotations,
    select_kofam_scores,
    parse_kos,
    format_kos,
    get_single_ko_definition,
    resolve_definition_for_kos,
    load_ko_definitions,
    read_fasta_ids,
)


class TestCharacterization(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_tsv(self, filename: str, content: str) -> Path:
        p = self.tmp / filename
        p.write_text(content.strip() + "\n", encoding="utf-8")
        return p

    def test_select_kofam_scores(self):
        # full score_type
        bs, ev, st = select_kofam_scores("full", 100.0, 1e-10, 50.0, 1e-5)
        self.assertEqual((bs, ev, st), (100.0, 1e-10, "full"))

        # domain score_type
        bs, ev, st = select_kofam_scores("domain", 100.0, 1e-10, 50.0, 1e-5)
        self.assertEqual((bs, ev, st), (50.0, 1e-5, "domain"))

        # missing domain score when domain required
        bs, ev, st = select_kofam_scores("domain", 100.0, 1e-10, np.nan, np.nan)
        self.assertTrue(np.isnan(bs))
        self.assertTrue(np.isnan(ev))
        self.assertEqual(st, "domain")

        # unrecognized score_type
        bs, ev, st = select_kofam_scores("other", 100.0, 1e-10, 50.0, 1e-5)
        self.assertTrue(np.isnan(bs))
        self.assertTrue(np.isnan(ev))
        self.assertEqual(st, "-")

    def test_parse_and_format_kos(self):
        self.assertEqual(parse_kos("K00001,K00002"), {"K00001", "K00002"})
        self.assertEqual(parse_kos("K00001; K00002"), {"K00001", "K00002"})
        self.assertEqual(parse_kos("-"), set())
        self.assertEqual(parse_kos(None), set())
        self.assertEqual(format_kos({"K00002", "K00001"}), "K00001,K00002")
        self.assertEqual(format_kos(set()), "-")

    def test_definitions(self):
        ko_defs = {"K00001": "alcohol dehydrogenase", "K00002": "alcohol dehydrogenase (NADP+)"}
        self.assertEqual(get_single_ko_definition("K00001", ko_defs), "alcohol dehydrogenase")
        self.assertEqual(get_single_ko_definition("K99999", ko_defs), "-")
        self.assertEqual(
            resolve_definition_for_kos({"K00001", "K00002"}, ko_defs),
            "K00001: alcohol dehydrogenase; K00002: alcohol dehydrogenase (NADP+)",
        )

    def test_unanimous_3_tools(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\tdeepkoala_annotate\ng1\tK00001\t0.95\t0.5\t*\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00001\t200.0\t1e-30\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk, eggnog_tsv=en)
        self.assertEqual(len(res), 1)
        row = res.iloc[0]
        self.assertEqual(row["gene_id"], "g1")
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["alternative_kos"], "-")
        self.assertEqual(row["consensus_level"], "unanimous")
        self.assertEqual(row["evidence"], "deepkoala,eggnog,kofam")

        ev = res.attrs["evidence"]
        self.assertEqual(len(ev), 3)
        self.assertTrue((ev["cross_method_status"] == "accepted").all())

    def test_majority_3_tools(self):
        # 2 agree on K00001, eggnog calls K00002
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\tdeepkoala_annotate\ng1\tK00001\t0.95\t0.5\t*\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00002\t200.0\t1e-30\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk, eggnog_tsv=en)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["alternative_kos"], "K00002")
        self.assertEqual(row["consensus_level"], "majority")
        self.assertEqual(row["evidence"], "deepkoala,kofam")

        ev = res.attrs["evidence"]
        k1_ev = ev[ev["ko"] == "K00001"]
        self.assertTrue((k1_ev["cross_method_status"] == "accepted").all())
        k2_ev = ev[ev["ko"] == "K00002"]
        self.assertEqual(k2_ev.iloc[0]["cross_method_status"], "alternative")

    def test_conflict_multiple_3_tools(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\tdeepkoala_annotate\ng1\tK00002\t0.95\t0.5\t*\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00003\t200.0\t1e-30\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk, eggnog_tsv=en, conflict_strategy="multiple")
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "-")
        self.assertEqual(row["alternative_kos"], "K00001,K00002,K00003")
        self.assertEqual(row["consensus_level"], "conflict")
        self.assertEqual(row["evidence"], "deepkoala,eggnog,kofam")

        ev = res.attrs["evidence"]
        self.assertTrue((ev["cross_method_status"] == "conflict").all())

    def test_conflict_priority_3_tools(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\tdeepkoala_annotate\ng1\tK00002\t0.95\t0.5\t*\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00003\t200.0\t1e-30\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk, eggnog_tsv=en, conflict_strategy="priority")
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["alternative_kos"], "K00002,K00003")
        self.assertEqual(row["consensus_level"], "conflict_priority")
        self.assertEqual(row["evidence"], "kofam")

    def test_conflict_priority_eggnog_over_deepkoala(self):
        # KOfam did not call, eggNOG K00003, DeepKOALA K00002
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\tdeepkoala_annotate\ng1\tK00002\t0.95\t0.5\t*\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00003\t200.0\t1e-30\n")

        res = integrate_annotations(deepkoala_tsv=dk, eggnog_tsv=en, conflict_strategy="priority")
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00003")
        self.assertEqual(row["alternative_kos"], "K00002")
        self.assertEqual(row["consensus_level"], "conflict_priority")
        self.assertEqual(row["evidence"], "eggnog")

    def test_single_tool_without_candidate(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        res = integrate_annotations(kofam_tsv=kf)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "single_tool")
        self.assertEqual(row["evidence"], "kofam")

    def test_single_tool_with_candidate(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\ng1\tK00001\t0.3\t0.5\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "single_tool_with_candidate")
        self.assertEqual(row["evidence"], "kofam,deepkoala(candidate)")

        ev = res.attrs["evidence"]
        dk_ev = ev[ev["method"] == "deepkoala"].iloc[0]
        self.assertEqual(dk_ev["cross_method_status"], "rescued")

    def test_zero_calls_kofam_rescue_alone(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\trescued\t80.0\t1e-6\t100.0\tfull\n")
        res = integrate_annotations(kofam_tsv=kf)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "single_tool")
        self.assertEqual(row["evidence"], "kofam(rescued)")

        ev = res.attrs["evidence"]
        self.assertEqual(ev.iloc[0]["cross_method_status"], "heuristic_rescued")

    def test_zero_calls_kofam_rescue_with_eggnog_candidate(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\trescued\t80.0\t1e-6\t100.0\tfull\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00001\t40.0\t1e-2\n")

        res = integrate_annotations(kofam_tsv=kf, eggnog_tsv=en)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "dual_candidate")
        self.assertEqual(row["evidence"], "eggnog(candidate),kofam(rescued)")

    def test_zero_calls_kofam_rescue_with_deepkoala_candidate(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\trescued\t80.0\t1e-6\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\ng1\tK00001\t0.3\t0.5\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "dual_candidate")
        self.assertEqual(row["evidence"], "deepkoala(candidate),kofam(rescued)")

    def test_zero_calls_kofam_rescue_with_both_candidates(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\trescued\t80.0\t1e-6\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\ng1\tK00001\t0.3\t0.5\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00001\t40.0\t1e-2\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk, eggnog_tsv=en)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "dual_candidate")
        self.assertEqual(row["evidence"], "deepkoala(candidate),eggnog(candidate),kofam(rescued)")

    def test_zero_calls_kofam_rescue_precedence_over_dk_en_candidates(self):
        # User specified: Kofam rescue in any form first (even alone) over DK & EN candidates
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\trescued\t80.0\t1e-6\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\ng1\tK00002\t0.3\t0.5\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00002\t40.0\t1e-2\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk, eggnog_tsv=en)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "single_tool")
        self.assertEqual(row["evidence"], "kofam(rescued)")

    def test_zero_calls_dk_and_en_candidate_agree(self):
        # No kofam rescue: DeepKOALA & eggNOG sub-threshold candidates agree
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\ng1\tK00002\t0.3\t0.5\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00002\t40.0\t1e-2\n")

        res = integrate_annotations(deepkoala_tsv=dk, eggnog_tsv=en)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00002")
        self.assertEqual(row["consensus_level"], "dual_candidate")
        self.assertEqual(row["evidence"], "deepkoala(candidate),eggnog(candidate)")

    def test_zero_calls_unannotated(self):
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\ng1\tK00001\t0.3\t0.5\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00002\t40.0\t1e-2\n")

        res = integrate_annotations(deepkoala_tsv=dk, eggnog_tsv=en)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "-")
        self.assertEqual(row["alternative_kos"], "-")
        self.assertEqual(row["consensus_level"], "unannotated")
        self.assertEqual(row["evidence"], "-")

    def test_eggnog_disambiguation(self):
        # eggNOG multi-KO K00001,K00002 disambiguated against KOfam calling K00001
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00001,K00002\t200.0\t1e-30\n")

        res = integrate_annotations(kofam_tsv=kf, eggnog_tsv=en, eggnog_filter_multi="disambiguate")
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["alternative_kos"], "K00002")
        self.assertEqual(row["consensus_level"], "unanimous")
        self.assertEqual(row["evidence"], "eggnog,kofam")

        ev = res.attrs["evidence"]
        k2_rec = ev[(ev["method"] == "eggnog") & (ev["ko"] == "K00002")].iloc[0]
        self.assertEqual(k2_rec["cross_method_status"], "disambiguated_dropped")

    def test_eggnog_filter_multi_none(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        en = self._write_tsv("en.tsv", "#query\tkegg_ko\tscore\tevalue\ng1\tK00001,K00002\t200.0\t1e-30\n")

        res = integrate_annotations(kofam_tsv=kf, eggnog_tsv=en, eggnog_filter_multi="none")
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        # Under "none", K00002 is not disambiguated away, so it stays in alternative_kos
        self.assertEqual(row["alternative_kos"], "K00002")

    def test_fasta_gene_universe_and_ordering(self):
        fasta = self._write_tsv("prots.faa", ">g2 some desc\nMKLL\n>g1 other desc\nMKLF\n>g3 third\nMKLA\n")
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")

        res = integrate_annotations(protein_fasta=fasta, kofam_tsv=kf)
        self.assertEqual(list(res["gene_id"]), ["g2", "g1", "g3"])
        g2_row = res[res["gene_id"] == "g2"].iloc[0]
        self.assertEqual(g2_row["consensus_level"], "unannotated")
        self.assertEqual(g2_row["accepted_ko"], "-")

    def test_export_tsv_and_evidence(self):
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\n")
        out_tsv = self.tmp / "output.tsv"

        res = integrate_annotations(kofam_tsv=kf, output_tsv=out_tsv)
        self.assertTrue(out_tsv.is_file())
        default_ev = self.tmp / "kolach_evidence.tsv"
        self.assertTrue(default_ev.is_file())

        out_lines = out_tsv.read_text(encoding="utf-8").splitlines()
        self.assertIn("accepted_ko", out_lines[0])
        self.assertIn("g1\tK00001", out_lines[1])

        ev_lines = default_ev.read_text(encoding="utf-8").splitlines()
        self.assertIn("cross_method_status", ev_lines[0])
        self.assertIn("accepted", ev_lines[1])

    def test_missing_files_raises_filenotfound(self):
        with self.assertRaises(FileNotFoundError):
            integrate_annotations(kofam_tsv=self.tmp / "nonexistent.tsv")

    def test_master_ko_list_precedence(self):
        db_dir = self.tmp / "db"
        (db_dir / "kofam").mkdir(parents=True)
        (db_dir / "kofam" / "ko_list").write_text("K00001\tMaster alcohol dehydrogenase\n", encoding="utf-8")

        kf = self._write_tsv(
            "kf.tsv",
            "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\tdefinition\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\tKOfam profile definition\n",
        )
        res = integrate_annotations(kofam_tsv=kf, database_dir=db_dir)
        row = res.iloc[0]
        self.assertEqual(row["definition"], "Master alcohol dehydrogenase")

    def test_deepkoala_simple_two_column(self):
        dk = self._write_tsv("dk.tsv", "name\tpredict_label\ng1\tK00001\n")
        res = integrate_annotations(deepkoala_tsv=dk)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "K00001")
        self.assertEqual(row["consensus_level"], "single_tool")
        ev = res.attrs["evidence"]
        self.assertEqual(ev.iloc[0]["score_type"], "upstream_filtered")

    def test_multi_ko_agreement_alternative_routing(self):
        # Tools agree on multi-KO: accepted_ko must be '-' according to strict single-KO policy
        kf = self._write_tsv("kf.tsv", "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\ng1\tK00001\tthreshold\t150.0\t1e-20\t100.0\tfull\ng1\tK00002\tthreshold\t160.0\t1e-20\t100.0\tfull\n")
        dk = self._write_tsv("dk.tsv", "gene_id\tpredict_label\tdeepkoala_score\tdeepkoala_threshold\tdeepkoala_annotate\ng1\tK00001,K00002\t0.95\t0.5\t*\n")

        res = integrate_annotations(kofam_tsv=kf, deepkoala_tsv=dk)
        row = res.iloc[0]
        self.assertEqual(row["accepted_ko"], "-")
        self.assertEqual(row["alternative_kos"], "K00001,K00002")
        self.assertEqual(row["consensus_level"], "unanimous")


if __name__ == "__main__":
    unittest.main()

