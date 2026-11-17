"""Verify captured fixtures against AnnoGuild with its real dependencies.

Usage:
  python scripts/verify_annoguild_parity.py --annoguild-dir ../AnnoGuild

The old evaluator imports pandas, numpy, networkx, click, graphviz and pytest.
Use an AnnoGuild checkout pinned to the commit recorded in the fixture. This
script performs no downloads and writes no files; any mismatch exits nonzero.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annoguild-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    fixture = json.loads((root / "tests/fixtures/annoguild_parity.json").read_text())
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.annoguild_dir, text=True).strip()
    if commit != fixture["source_commit"]:
        parser.error(f"Check out AnnoGuild commit {fixture['source_commit']} before verification")

    import pandas as pd

    # Load the actual upstream code, including its original third-party imports.
    source = args.annoguild_dir / "scripts/paths_graph_parser.py"
    spec = importlib.util.spec_from_file_location("annoguild_reference", source)
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    pathways = pd.read_csv(args.annoguild_dir / "resources/pathways_refined.tsv", sep="\t", dtype=str).fillna("")
    reactions = pd.read_csv(args.annoguild_dir / "resources/custom_input_modules/EMERGE_pathways_module.tsv", sep="\t", dtype=str)
    reactions["number"] = reactions["module_name"].str.split("number:").str[1]
    reactions = reactions.set_index(["complex", "number"])
    import re
    all_kos = set(re.findall(r"K\d{5,6}", " ".join(reactions["definition"]) + " " + " ".join(pathways["signature_definition"])))
    compared = 0
    for profile in fixture["profiles"]:
        features = set(profile["kos"])
        dist = pd.Series({ko: int(ko in features) for ko in all_kos | features}, name=profile["name"])
        for _, pathway in pathways.iterrows():
            tree = legacy.RuleParser(dist)
            tree.make_path_graph(pathway["reaction"], reactions.loc[pathway["pathway"], "definition"])
            signature_ok = True
            if pathway["signature_definition"]:
                tree.make_signiture_graph(pathway["signature_definition"])
                signature_ok = bool(tree.check_signature())
            name = pathway["pathway"] + "-" + pathway["subpathway"]
            called = bool(signature_ok and tree.check_presance() and tree.check_ko_percent())
            assert called == profile["product"][name], (profile["name"], name, "called")
            try:
                ko_coverage = float(tree.get_ko_percent())
            except ValueError:
                ko_coverage = 0.0  # Legacy diagnostic getter fails on empty options.
            expected = profile["metrics"][name]
            assert abs(float(tree.get_presance()) - expected["reaction_coverage"]) < 1e-12
            assert abs(ko_coverage - expected["ko_coverage"]) < 1e-12
            assert signature_ok == expected["signature_present"]
            compared += 1
    print(f"Verified {compared} captured pathway calls and diagnostics against original AnnoGuild libraries.")


if __name__ == "__main__":
    main()
