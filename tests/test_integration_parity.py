"""Compare the refactor with upstream db5a5ce, including exhaustive call categories."""
import gzip
import inspect
import itertools
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from kolach import integrate as current
from tests.reference import integrate_db5a5ce as previous

TOOLS = ("kofam", "deepkoala", "eggnog")
KOS = ("K00001", "K00002", "K00003")
STATUSES = ("threshold_passing", "below_threshold", "unknown")


def subsets(values):
    return [
        set(value for value, enabled in zip(values, mask) if enabled)
        for mask in itertools.product((False, True), repeat=len(values))
    ]


def records_for(tool, statuses, kos=KOS):
    records = []
    for ko, status in zip(kos, statuses):
        if status is None:
            continue
        records.append(previous.make_evidence_record(
            "gene", ko, tool, status, bit_score=100.0, e_value=0.0,
            domain_bit_score=80.0, domain_e_value=1e-10,
            score_type="full", threshold=60.0, assignment="threshold",
            deepkoala_probability=0.9, deepkoala_threshold=0.5,
            eggnog_sources=[{"group": ko, "is_multi": False, "status": status}],
        ))
    return records


def previous_decision(gene, filtering, strategy):
    # Build the oracle's indexes independently. Reusing the refactor's indexes
    # here would hide a categorization error in both expected and actual results.
    tool_records = {tool: [] for tool in gene.active_tools}
    for record in gene.records:
        if record["method"] in tool_records:
            tool_records[record["method"]].append(record)
    confident, candidates, rescued = {}, {}, {}
    for tool, records in tool_records.items():
        for status, target in [
            ("threshold_passing", confident), ("below_threshold", candidates),
            ("heuristic_rescued", rescued),
        ]:
            kos = {record["ko"] for record in records if record["original_status"] == status}
            if kos:
                target[tool] = kos
    retained, dropped, narrowed, unresolved = previous.disambiguate_eggnog(
        set(confident.get("eggnog", set())), tool_records, gene.active_tools,
        confident, candidates, rescued, filtering,
    )
    if retained:
        confident["eggnog"] = retained
    decision = previous.adjudicate_gene(
        gene.active_tools, confident, candidates, rescued, dropped, strategy,
    )
    accepted, alternatives, level, evidence = decision
    audit = previous.annotate_evidence_records(
        gene.gene_id, gene.records, tool_records, gene.active_tools,
        {accepted} if accepted != "-" else set(), alternatives, accepted, level,
        dropped, narrowed, unresolved,
    )
    return decision, (retained, dropped, narrowed, unresolved), audit


def assert_gene_parity(gene, filtering, strategy, case):
    expected, groups, expected_audit = previous_decision(gene, filtering, strategy)
    current._disambiguate_gene(gene, filtering)
    current._resolve_gene(gene, strategy)
    actual = gene.accepted_ko, gene.alternatives, gene.consensus_level, gene.evidence
    actual_groups = gene.eggnog_retained, gene.eggnog_dropped, gene.narrowed_groups, gene.unresolved_groups
    actual_audit = current._build_audit(gene)
    # Both versions reuse the same record metrics; compare every audit field.
    if actual != expected or actual_groups != groups or actual_audit != expected_audit:
        raise AssertionError(
            f"case={case!r}, filtering={filtering}, strategy={strategy}\n"
            f"decision: {actual!r} != {expected!r}\n"
            f"groups: {actual_groups!r} != {groups!r}\n"
            f"audit: {actual_audit!r} != {expected_audit!r}"
        )


def build_gene(records, active=TOOLS):
    gene = current.GeneAnnotation("gene", records, list(active))
    current._categorize_gene(gene)
    return gene


class ExhaustiveCategories(unittest.TestCase):
    def test_every_three_ko_category_assignment(self):
        profiles = []
        for tool in TOOLS:
            statuses = (None,) + STATUSES
            if tool == "kofam":
                statuses += ("heuristic_rescued",)
            # None means inactive; an all-absent tuple means active with no calls.
            choices = [(None, [])]
            for combination in itertools.product(statuses, repeat=len(KOS)):
                choices.append((combination, records_for(tool, combination)))
            profiles.append(choices)
        count = 0
        for configuration in itertools.product(*profiles):
            active = [tool for tool, (statuses, _) in zip(TOOLS, configuration) if statuses is not None]
            if not active:
                continue
            records = []
            for _, tool_records in configuration:
                records.extend(tool_records)
            case = tuple(statuses for statuses, _ in configuration)
            for filtering, strategy in itertools.product(("none", "disambiguate"), ("multiple", "priority")):
                assert_gene_parity(build_gene(records, active), filtering, strategy, case)
                count += 1
        self.assertEqual(count, 2_129_396)

    def test_every_two_row_eggnog_source_configuration(self):
        choices = [(set(), None)]
        for kos in subsets(KOS)[1:]:
            for status in STATUSES:
                choices.append((kos, status))
        count = 0
        for first, second in itertools.product(choices, repeat=2):
            sources = {}
            for index, (kos, status) in enumerate((first, second)):
                for ko in kos:
                    sources.setdefault(ko, []).append({
                        "group": index, "is_multi": len(kos) > 1, "status": status,
                    })
            eggnog = []
            rank = {"threshold_passing": 2, "below_threshold": 1, "unknown": 0}
            for ko, ko_sources in sources.items():
                status = max((source["status"] for source in ko_sources), key=rank.get)
                record = records_for("eggnog", [status], [ko])[0]
                record["eggnog_sources"] = ko_sources
                record["eggnog_shared_seed_hit"] = any(source["is_multi"] for source in ko_sources)
                eggnog.append(record)
            for trusted, candidates in itertools.product(subsets(KOS), repeat=2):
                records = eggnog + records_for(
                    "kofam", ["threshold_passing"] * len(trusted), sorted(trusted),
                ) + records_for(
                    "deepkoala", ["below_threshold"] * len(candidates), sorted(candidates),
                )
                for filtering, strategy in itertools.product(("none", "disambiguate"), ("multiple", "priority")):
                    case = first, second, trusted, candidates
                    assert_gene_parity(build_gene(records), filtering, strategy, case)
                    count += 1
        self.assertEqual(count, 123_904)

    def test_larger_and_repeated_groups(self):
        kos = tuple(f"K{index:05}" for index in range(1, 7))
        records = records_for("eggnog", ["threshold_passing"] * len(kos), kos)
        for record in records:
            record["eggnog_sources"] = [
                {"group": group, "is_multi": True, "status": "threshold_passing"}
                for group in range(4)
            ]
        records.extend(records_for("kofam", ["heuristic_rescued"], kos[:1]))
        records.extend(records_for("deepkoala", ["threshold_passing"] * 2, kos[1:3]))
        for filtering, strategy in itertools.product(("none", "disambiguate"), ("multiple", "priority")):
            assert_gene_parity(build_gene(records), filtering, strategy, "six KOs, four repeated groups")


class PipelineParity(unittest.TestCase):
    def compare(self, **kwargs):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            old_path, new_path = directory / "old/summary.tsv", directory / "new/summary.tsv"
            expected = previous.integrate_annotations(**kwargs, output_tsv=old_path)
            actual = current.integrate_annotations(**kwargs, output_tsv=new_path)
            pd.testing.assert_frame_equal(actual, expected)
            pd.testing.assert_frame_equal(actual.attrs["evidence"], expected.attrs["evidence"])
            self.assertEqual(new_path.read_bytes(), old_path.read_bytes())
            self.assertEqual(
                (new_path.parent / "kolach_evidence.tsv").read_bytes(),
                (old_path.parent / "kolach_evidence.tsv").read_bytes(),
            )
        return actual

    def test_all_example_tool_and_policy_combinations(self):
        examples = Path(__file__).resolve().parents[1] / "examples"
        for active in subsets(TOOLS)[1:]:
            for filtering, strategy in itertools.product(("none", "disambiguate"), ("multiple", "priority")):
                with self.subTest(active=active, filtering=filtering, strategy=strategy):
                    kwargs = {tool + "_tsv": examples / (tool + "_annotations.tsv") for tool in active}
                    self.compare(**kwargs, eggnog_filter_multi=filtering, conflict_strategy=strategy)

    def test_public_signatures(self):
        for name, function in inspect.getmembers(previous, inspect.isfunction):
            if not name.startswith("_"):
                # NaN defaults do not compare equal, even for identical signatures.
                self.assertEqual(str(inspect.signature(getattr(current, name))), str(inspect.signature(function)), name)

    def test_empty_and_fasta_ordering(self):
        self.compare()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            fasta = directory / "genes.faa.gz"
            with gzip.open(fasta, "wt") as stream:
                stream.write(">b description\nM\n>a\nM\n>b duplicate\nM\n")
            kofam = directory / "kofam.tsv"
            for content in ("", "gene_id\tko\n"):
                kofam.write_text(content)
                self.compare(protein_fasta=fasta, kofam_tsv=kofam)
            kofam.write_text("gene_id\tko\tassignment\nextra\tK00001\tthreshold\n")
            self.compare(protein_fasta=fasta, kofam_tsv=kofam)

    def test_metric_boundaries_duplicates_and_definition_precedence(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            kofam = directory / "kofam.tsv"
            kofam.write_text(
                "gene_id\tko\tassignment\tbit_score\te_value\tdomain_bit_score\tdomain_e_value\tthreshold\tscore_type\tdefinition\n"
                "gene\tK00001\tbelow\t200\t0\t80\t0\t100\tfull\tweak\n"
                "gene\tK00001\tthreshold\t100\t1e-5\t50\t1e-3\t100\tfull\tstrong\n"
                "gene\tK00001\tthreshold\t100\t1e-5\t50\t1e-3\t100\tfull\ttied\n"
                "gene\tK00002\t-\t300\t0\t-\t-\t100\tdomain\tmissing domain\n"
                "gene\tK00003\trescued\t80\t1e-6\t40\t1e-2\t100\tfull\trescue\n"
                "gene\tK00004\t-\tbad\tinf\t-\t-\t100\tfull\tunknown\n"
            )
            deepkoala = directory / "deepkoala.tsv"
            deepkoala.write_text(
                "name\tpredict_label\tprobability\tthreshold\tannotate\n"
                "gene\tK00001\t0.4\t0.5\t*\n"
                "gene\tK00001\t0.9\t0.95\t-\n"
                "gene\tK00002\t0.5\t0.5\t-\n"
                "gene\tK00003\tbad\t0.5\t-\n"
            )
            eggnog = directory / "eggnog.tsv.gz"
            with gzip.open(eggnog, "wt") as stream:
                stream.write(
                    "## metadata\n# comment\n#query\tkegg_ko\tscore\tevalue\tdescription\n"
                    "gene\tK00001,K00002\t60\t1e-5\tgeneric\n"
                    "gene\tK00001\t90\t1e-3\tweak duplicate\n"
                    "gene\tK00002\t60\t0\tsingleton protection\n"
                    "gene\tK00003\tbad\t0\tunknown\n# footer\n"
                )
            (directory / "kofam").mkdir()
            (directory / "kofam/ko_list").write_text("K00001\tMaster definition\n")
            (directory / "ko_list").write_text("K00001\tLower precedence\n")
            for filtering, strategy in itertools.product(("none", "disambiguate"), ("multiple", "priority")):
                result = self.compare(
                    kofam_tsv=kofam, deepkoala_tsv=deepkoala, eggnog_tsv=eggnog,
                    database_dir=directory, eggnog_filter_multi=filtering,
                    conflict_strategy=strategy, heuristic_bitscore_fraction=0.8,
                    heuristic_e_value=1e-6,
                )
                self.assertIn("Master definition", result.iloc[0]["definition"])

    def test_validation_parity(self):
        cases = [
            {"conflict_strategy": "bad"}, {"eggnog_filter_multi": "bad"},
            {"eggnog_min_bitscore": -1}, {"eggnog_min_bitscore": np.inf},
            {"eggnog_max_evalue": -1}, {"heuristic_e_value": np.nan},
            {"heuristic_bitscore_fraction": 1.1}, {"eggnog_min_bitscore": "bad"},
            {"protein_fasta": "missing-fasta"}, {"deepkoala_tsv": "missing-table"},
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.tsv"
            for tool, content in [
                ("kofam", "gene_id\n"), ("deepkoala", "wrong\tko\n"),
                ("eggnog", "#query\tscore\n"),
                ("kofam", "gene_id\tko\n \tK00001\n"),
                ("deepkoala", "name\tpredict_label\n \tK00001\n"),
                ("eggnog", "#query\tkegg_ko\n \tK00001\n"),
            ]:
                path.write_text(content)
                self.assert_same_error({tool + "_tsv": path})
        for kwargs in cases:
            self.assert_same_error(kwargs)

    def assert_same_error(self, kwargs):
        with self.assertRaises(Exception) as old:
            previous.integrate_annotations(**kwargs)
        with self.assertRaises(type(old.exception)) as new:
            current.integrate_annotations(**kwargs)
        self.assertEqual(str(new.exception), str(old.exception))


if __name__ == "__main__":
    unittest.main()
