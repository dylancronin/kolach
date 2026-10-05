"""Behavioral tests, including 1,485 pinned AnnoGuild pathway comparisons.

Core feature tests require only Python's standard library. The optional workflow
test uses installed Snakemake; existing integration and BRITE tests still use
Kolach's normal pandas/numpy dependencies.
"""
import ast
import csv
import io
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from kolach.metabolism.evaluate import (
    Expression, PathwayRule, evaluate_pathway, load_rules,
    parse_expression, signature_present,
)
from kolach.metabolism.specialise import (
    add_markers, classify_specialisation, ko_predict_methanogen,
    parse_bool, read_annotations, read_product, specialise,
)


ROOT = Path(__file__).resolve().parents[1]


class TestEvaluator(unittest.TestCase):
    def test_all_legacy_pathway_calls_and_coverage(self):
        fixture = json.loads((ROOT / "tests/fixtures/annoguild_parity.json").read_text())
        rules = load_rules()
        self.assertEqual(len(rules), 55)
        for profile in fixture["profiles"]:
            for rule in rules:
                with self.subTest(profile=profile["name"], pathway=rule.pathway_id):
                    result = evaluate_pathway(rule, set(profile["kos"]))
                    self.assertEqual(result["called"], profile["product"][rule.pathway_id])
                    for metric, expected in profile["metrics"][rule.pathway_id].items():
                        if metric == "signature_present" and rule.pathway_id == "nitrogen_redox-comammox":
                            # The old graph merges sibling OR groups under an AND.
                            # Correct strict signature semantics are tested below.
                            continue
                        self.assertAlmostEqual(result[metric], expected)

    def test_boolean_signature_and_precedence(self):
        expression = parse_expression(" (K00001 + K00002), K00003 ")
        self.assertFalse(signature_present(expression, {"K00001"}))
        self.assertTrue(signature_present(expression, {"K00001", "K00002"}))
        self.assertTrue(signature_present(expression, {"K00003"}))

    def test_invalid_definitions(self):
        for definition in ["", "()", "K00001+", "K00001,,K00002", "(K00001", "K00001)", "K00001 K00002", "K00001/2"]:
            with self.subTest(definition=definition), self.assertRaises(ValueError):
                parse_expression(definition)

    def test_partial_reaction_and_absent_reaction_denominator(self):
        # Two of three reactions have some evidence (2/3). Only those reactions
        # contribute to KO coverage: two observed subunits out of three (2/3).
        rule = PathwayRule("example", {"1": parse_expression("K00001+K00002"),
                                      "2": parse_expression("K00003"),
                                      "3": parse_expression("K00004")}, (("1", "2", "3"),), None)
        result = evaluate_pathway(rule, {"K00001", "K00003"}, 2/3, 2/3)
        self.assertTrue(result["called"])
        self.assertAlmostEqual(result["reaction_coverage"], 2/3)
        self.assertAlmostEqual(result["ko_coverage"], 2/3)
        self.assertEqual(result["missing_features"], "K00002;K00004")
        self.assertFalse(evaluate_pathway(rule, {"K00001", "K00003"}, 0.7, 0.6)["called"])

    def test_repeated_kos_count_per_reaction(self):
        rule = PathwayRule("example", {"1": parse_expression("K00001+K00002"),
                                      "2": parse_expression("K00001+K00003+K00004")}, (("1", "2"),), None)
        self.assertAlmostEqual(evaluate_pathway(rule, {"K00001"})["ko_coverage"], 2/5)

    def test_signature_blocks_otherwise_complete_pathway(self):
        rule = PathwayRule("example", {"1": parse_expression("K00001")}, (("1",),), parse_expression("K00002"))
        self.assertFalse(evaluate_pathway(rule, {"K00001"})["called"])

    def test_signature_sibling_or_groups_are_independent(self):
        expression = parse_expression("(K00001,K00002)+(K00003,K00004)")
        self.assertFalse(signature_present(expression, {"K00001"}))
        self.assertTrue(signature_present(expression, {"K00001", "K00003"}))

    def test_legacy_independent_route_maxima_are_preserved(self):
        rule = PathwayRule("example", {"1": parse_expression("K00001+K00002+K00003"),
                                      "2": parse_expression("K00004+K00005+K00006"),
                                      "3": parse_expression("K00007"),
                                      "4": parse_expression("K00008")}, (("1", "2"), ("3", "4")), None)
        result = evaluate_pathway(rule, {"K00001", "K00004", "K00007"})
        self.assertTrue(result["called"])
        self.assertEqual(result["reaction_coverage"], 1)
        self.assertEqual(result["ko_coverage"], 1)
        self.assertEqual(result["ko_coverage_route"], "3+4")


class TestClassifier(unittest.TestCase):
    def test_primary_guilds(self):
        cases = [
            ({}, set(), ""),
            ({"hydrogenotrophic_methanogenesis-all": True}, set(), "methanogen"),
            ({"methylotrophic_methanogenesis-h2_dependent": True}, set(), "methanogen"),
            ({"methylotrophic_methanogenesis-h2_independent": True}, set(), "methanogen"),
            ({"carbon_redox-dissimilatory_methanotrophy": True}, {"K109440"}, "methanotroph"),
            ({"Wood_Ljungdahl-acetogen": True}, set(), "homoacetogen"),
            ({"CAZy-A": True, "CAZy-B": True, "CAZy-C": True, "a_degradation": True,
              "b_degradation": True, "c_degradation": True}, set(), "generalist"),
            ({"ethanol_fermentation-all": True, "lactate_fermentation-all": True}, set(), "fermenter"),
            ({"CAZy-Starch": True}, set(), "macromolecule_degrader"),
            ({"fructose_degradation-all": True}, set(), "monomer_degrader"),
        ]
        for calls, kos, expected in cases:
            with self.subTest(expected=expected):
                self.assertEqual(classify_specialisation(calls, kos)["specialisation"], expected)

    def test_corrected_methyltransferase_literals(self):
        for ko in ["K04480", "K14082"]:
            self.assertTrue(ko_predict_methanogen({"K00399", "K16183", ko}))
        self.assertFalse(ko_predict_methanogen({"K00399", "K16183"}))

    def test_methanotrophy_exclusion_covers_pathway_branch(self):
        calls = {"hydrogenotrophic_methanogenesis-all": True,
                 "carbon_redox-dissimilatory_methanotrophy": True}
        self.assertEqual(classify_specialisation(calls, {"K109440"})["specialisation"], "methanotroph")

    def test_nitrogen_and_sulfur_priority(self):
        calls = {"nitrogen_redox-comammox": True, "nitrogen_redox-nitrogen_fixation": True,
                 "sulfur_redox-dissimilatory_sulfate_reduction": True,
                 "sulfur_redox-dissimilatory_sulfur_oxidation": True}
        labels = classify_specialisation(calls, set())
        self.assertEqual(labels["nitrogen_specialisation"], "comammox")
        self.assertEqual(labels["sulfur_specialisation"], "sulfate_reducer")

    def test_explicit_booleans(self):
        self.assertFalse(parse_bool("False"))
        self.assertTrue(parse_bool("True"))
        self.assertEqual(classify_specialisation({"CAZy-A": "False"}, set())["specialisation"], "")
        with self.assertRaises(ValueError):
            parse_bool("0.6")

    def test_against_corrected_original_classifier(self):
        # Execute the original classifier functions/constants, correcting only
        # the documented comma, pathway names, and precedence defects. pd.Series
        # is a dictionary-returning adapter because no pandas operation is used
        # inside these functions; all original control flow remains unchanged.
        source = (ROOT / "tests/fixtures/legacy_calculate_specialisation.py").read_text()
        source = source.replace('"K04480"                            #', '"K04480",                           #')
        source = source.replace('"methylotrophic_methanogenesis-h2_dep"', '"methylotrophic_methanogenesis-h2_dependent"')
        source = source.replace('"methylotrophic_methanogenesis-h2_indep"', '"methylotrophic_methanogenesis-h2_independent"')
        source = source.replace('any(p in called for p in [', '(any(p in called for p in [')
        source = source.replace('or row["ko_methanogen"]', 'or row["ko_methanogen"])')
        nodes = [node for node in ast.parse(source).body
                 if isinstance(node, ast.Assign) or isinstance(node, ast.FunctionDef)
                 and node.name in {"ko_predict_methanogen", "ko_predict_fdh", "classify_specialisation"}]
        namespace = {"pd": SimpleNamespace(Series=lambda values: values)}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), "<corrected AnnoGuild classifier>", "exec"), namespace)
        fixture = json.loads((ROOT / "tests/fixtures/annoguild_parity.json").read_text())
        for profile in fixture["profiles"]:
            kos = set(profile["kos"])
            row = {**profile["product"], "ko_methanogen": namespace["ko_predict_methanogen"](kos),
                   "ko_fdh": namespace["ko_predict_fdh"](kos)}
            expected = {key: value if value is not None else "" for key, value in namespace["classify_specialisation"](row).items()}
            with self.subTest(profile=profile["name"]):
                self.assertEqual(classify_specialisation(profile["product"], kos), expected)


class TestInputsAndOutputs(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.tmp = Path(self.directory.name)

    def write(self, name, text):
        path = self.tmp / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_golden_example_outputs(self):
        paths = specialise(ROOT / "examples/specialisation/annotations.tsv", self.tmp,
                           markers=ROOT / "examples/specialisation/markers.tsv",
                           cazy_product=ROOT / "examples/specialisation/cazy_product.tsv")
        for name in ["product", "specialisations"]:
            expected = ROOT / "examples/specialisation/expected" / paths[name].name
            self.assertEqual(paths[name].read_bytes(), expected.read_bytes())

    def test_single_genome_and_ambiguous_alternatives(self):
        table = self.write("a.tsv", "gene_id\taccepted_ko\talternative_kos\ng1\t-\tK00844\n")
        paths = specialise(table, self.tmp / "out", genome_id="MAG")
        _, rows = read_product(paths["product"])
        self.assertFalse(rows["MAG"]["fructose_degradation-all"])
        self.assertEqual(list(rows), ["MAG"])
        with self.assertRaises(ValueError):
            read_annotations(table)

    def test_unannotated_and_inventory_only_genomes(self):
        table = self.write("a.tsv", "gene_id\taccepted_ko\ng1\t-\n")
        mapping = self.write("map.tsv", "gene_id\tgenome\ng1\tMAG1\ng2\tMAG2\n")
        paths = specialise(table, self.tmp / "out", gene_genome_map=mapping)
        self.assertEqual(list(read_product(paths["product"])[1]), ["MAG1", "MAG2"])
        empty = self.write("empty.tsv", "gene_id\taccepted_ko\n")
        paths = specialise(empty, self.tmp / "empty", genome_id="MAG")
        self.assertEqual(list(read_product(paths["product"])[1]), ["MAG"])

    def test_missing_duplicate_and_invalid_input(self):
        cases = ["gene_id\taccepted_ko\tgenome\ng1\tK00844\tMAG\ng1\t-\tMAG\n",
                 "gene_id\taccepted_ko\tgenome\ng1\tK0084400\tMAG\n",
                 "gene_id\taccepted_ko\tgenome\ng1\tK00844,K00399\tMAG\n",
                 "gene_id\taccepted_ko\tgenome\ng1\t-\t\n",
                 "gene_id\taccepted_ko\tgenome\ng1\t-\n"]
        for text in cases:
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_annotations(self.write("a.tsv", text))

    def test_mapping_checks(self):
        table = self.write("a.tsv", "gene_id\taccepted_ko\ng1\t-\n")
        for text in ["gene_id\tgenome\ng2\tMAG\n", "gene_id\tgenome\ng1\tMAG\ng1\tOTHER\n"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_annotations(table, gene_genome_map=self.write("map.tsv", text))

    def test_markers_replace_parent_and_avoid_reverse_dsr_substring_bug(self):
        genomes = {"MAG": {"g1": {"K10944"}, "g2": {"K11180"}}}
        add_markers(genomes, self.write("m.tsv", "gene_id\tfeature\ng1\tmarker:amoA\ng2\trdsrA\n"))
        self.assertEqual(genomes["MAG"]["g1"], {"K109441"})
        self.assertEqual(genomes["MAG"]["g2"], {"K111801"})

    def test_markers_reject_conflicts_unknown_genes_and_ambiguous_ids(self):
        for text in ["gene_id\tfeature\ng1\tpmoA\ng1\tamoA\n",
                     "gene_id\tfeature\nunknown\tpmoA\n", "gene_id\tfeature\ng1\tunknown\n"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                add_markers({"MAG": {"g1": {"K10944"}}}, self.write("m.tsv", text))
        with self.assertRaises(ValueError):
            add_markers({"MAG": {"g1": set()}, "OTHER": {"g1": set()}}, self.write("m.tsv", "gene_id\tfeature\ng1\tpmoA\n"))

    def test_import_existing_product_and_cazy(self):
        annotation = self.write("a.tsv", "gene_id\tko_id\tgenome\ng1\tK109440\tMAG\n")
        product = self.write("p.tsv", "#genome\tcarbon_redox-dissimilatory_methanotrophy\tCAZy-A\nMAG\tTrue\tFalse\n")
        paths = specialise(annotation, self.tmp / "out", ko_column="ko_id", product_table=product)
        self.assertIn("MAG\tmethanotroph\tNone detected\t\n", paths["specialisations"].read_text())
        self.assertEqual(len(paths["pathways"].read_text().splitlines()), 1)
        bad_cazy = self.write("c.tsv", "genome\tCAZy-A\nOTHER\tTrue\n")
        with self.assertRaises(ValueError):
            specialise(annotation, self.tmp / "bad", ko_column="ko_id", cazy_product=bad_cazy)

    def test_unassessed_evidence_is_separate_from_legacy_schema(self):
        paths = specialise(self.write("a.tsv", "gene_id\taccepted_ko\ng1\t-\n"), self.tmp / "out", genome_id="MAG")
        self.assertEqual(paths["specialisations"].read_text().splitlines()[0],
                         "genome\tspecialisation\tnitrogen_specialisation\tsulfur_specialisation")
        self.assertIn("unassessed", paths["assessment"].read_text())
        self.assertIn("unresolved_refinement", paths["pathways"].read_text())

    def test_custom_rules_and_threshold_validation(self):
        reactions = self.write("r.tsv", "pathway\treaction_id\tdefinition\nexample\t1\tK00844\n")
        pathways = self.write("p.tsv", "pathway\tsubpathway\treaction\tsignature_definition\nexample\tall\t1\t\n")
        paths = specialise(self.write("a.tsv", "gene_id\taccepted_ko\ng1\tK00844\n"), self.tmp / "out",
                           genome_id="MAG", reactions_file=reactions, pathways_file=pathways)
        self.assertIn("MAG\tTrue", paths["product"].read_text())
        for threshold in [-1, 2, float("nan"), float("inf")]:
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                specialise("unused", self.tmp, reaction_threshold=threshold)
        pathways.write_text("pathway\tsubpathway\treaction\tsignature_definition\nexample\tall\t2\t\n")
        with self.assertRaises(ValueError):
            load_rules(reactions, pathways)

    def test_cli_example_and_validation(self):
        from kolach.cli import main
        argv = ["kolach", "specialise", "--annotation-table", str(ROOT / "examples/specialisation/annotations.tsv"),
                "--output-dir", str(self.tmp / "cli")]
        with patch.object(sys, "argv", argv), patch("sys.stdout", new_callable=io.StringIO):
            main()
        self.assertTrue((self.tmp / "cli/kolach_specialisations.tsv").is_file())
        for argv in [["kolach", "annotate", "--protein-fasta", "a.faa", "--database-dir", "db", "--add-specialisations"],
                     ["kolach", "specialise", "--annotation-table", "a.tsv", "--output-dir", "out", "--ko-threshold", "nan"]]:
            with patch.object(sys, "argv", argv), patch("sys.stderr", new_callable=io.StringIO), self.assertRaises(SystemExit) as error:
                main()
            self.assertEqual(error.exception.code, 2)

    def test_annotate_cli_passes_workflow_configuration(self):
        from kolach.cli import main
        argv = ["kolach", "annotate", "--protein-fasta", "a.faa", "--database-dir", "db",
                "--add-specialisations", "--genome-id", "MAG", "--markers", "m.tsv"]
        with patch.object(sys, "argv", argv), patch("subprocess.run", return_value=SimpleNamespace(returncode=0)) as run:
            main()
        command = run.call_args.args[0]
        self.assertIn("add_specialisations=True", command)
        self.assertIn('genome_id=["MAG"]', command)
        self.assertIn('markers="m.tsv"', command)

    def test_annotate_preserves_numeric_genome_id(self):
        from kolach.cli import main
        argv = ["kolach", "annotate", "--protein-fasta", "a.faa", "--database-dir", "db",
                "--add-specialisations", "--genome-id", "001"]
        with patch.object(sys, "argv", argv), patch("subprocess.run", return_value=SimpleNamespace(returncode=0)) as run:
            main()
        config_value = next(arg.split("=", 1)[1] for arg in run.call_args.args[0]
                            if arg.startswith("genome_id="))
        self.assertEqual(json.loads(config_value), ["001"])


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is not installed")
class TestSpecialisationWorkflow(unittest.TestCase):
    def test_numeric_ids_and_grouping_changes_in_real_workflow(self):
        """Use captured annotations; exercise integration and report jobs offline."""
        from kolach.cli import main
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            out = work / "out"
            out.mkdir()
            (work / "proteins.faa").write_text(">g1\nMAAAA\n")
            (out / "kofam_annotations.tsv").write_text(
                "gene_id\tko\tassignment\tbit_score\te_value\tthreshold\tscore_type\n"
                "g1\tK10944\tthreshold\t150\t1e-20\t100\tfull\n")
            (work / "markers.tsv").write_text("gene_id\tfeature\ng1\tpmoA\n")

            def command(genome):
                argv = ["kolach", "annotate", "--protein-fasta", str(work / "proteins.faa"),
                        "--database-dir", str(work / "db"), "--output-dir", str(out),
                        "--add-specialisations", "--genome-id", genome, "--markers", "markers.tsv"]
                with patch.object(sys, "argv", argv), patch("subprocess.run", return_value=SimpleNamespace(returncode=0)) as run:
                    main()
                # Limit execution to existing annotation fixtures and downstream jobs.
                return [sys.executable, "-m", "snakemake", *run.call_args.args[0][1:],
                        "--directory", str(work), "--allowed-rules",
                        "all", "integrate_annotations", "specialise_genomes"]

            def execute(args):
                result = subprocess.run(args, capture_output=True, text=True, timeout=90)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                return result.stdout + result.stderr

            execute(command("001"))
            product = out / "kolach_product_refined.tsv"
            self.assertEqual(list(read_product(product)[1]), ["001"])
            metadata = json.loads((out / "kolach_specialisation_metadata.json").read_text())
            self.assertEqual(metadata["inputs"]["markers"]["path"], str(work / "markers.tsv"))
            self.assertIn("Nothing to be done", execute(command("001") + ["--dry-run"]))
            self.assertIn("Params have changed", execute(command("002") + ["--dry-run"]))
            execute(command("002"))
            self.assertEqual(list(read_product(product)[1]), ["002"])


if __name__ == "__main__":
    unittest.main()
