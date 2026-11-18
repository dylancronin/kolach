"""Compare Kolach specialisation with AnnoGuild's SynComm run.

Logic parity (same inputs as AnnoGuild):
  L1  Kolach's evaluator on AnnoGuild's own KO table (after its GraftM step)
      versus EMERGE_20260226_product_refined.tsv (55 pathways x 6 genomes).
  L2  Kolach's classifier on AnnoGuild's product_refined_with_CAZy versus
      EMERGE_20260226_genome_specialisations.tsv.
End to end (Kolach's own evidence), when kolach_output/ exists:
  E1  GraftM markers: Kolach (gated, exact lookup) versus AnnoGuild's
      substring-based KO replacement.
  E2  Accepted KO per gene versus AnnoGuild's ko_id.
  E3  Pathway calls and guild labels, with differing evidence listed.

Usage: PYTHONPATH=<checkout>/src python compare_syncomm.py
Writes syncomm_comparison.json and syncomm_*.tsv next to this script.
"""
import csv
import json
from pathlib import Path
import sys
import tempfile

from kolach.metabolism.specialise import (
    MARKER_IDS, classify_specialisation, parse_kos, read_product, specialise,
)

HERE = Path(__file__).resolve().parent
ANNOGUILD = Path("/fs/scratch/PDS0325/Cronin/syncomm/annoguild")
KOLACH = HERE / "syncomm" / "kolach_output"
LEGACY_MARKERS = {value: name for name, value in MARKER_IDS.items()}


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def legacy_annotations(tmp):
    """AnnoGuild's post-GraftM KO table as gene_id/genome/ko_id."""
    # AnnoGuild's DRAM/MicrobeAnnotator merge leaves some genes on two rows;
    # its distillation uses KO presence, so merge each gene's KOs.
    merged = {}
    for r in read_rows(ANNOGUILD / "all_annotations_filled_graftm.tsv"):
        entry = merged.setdefault(r["Unnamed: 0"], {"Unnamed: 0": r["Unnamed: 0"], "fasta": r["fasta"], "kos": set()})
        entry["kos"] |= parse_kos(r["ko_id"])
    rows = [{**e, "ko_id": ",".join(sorted(e.pop("kos"))) or "-"} for e in merged.values()]
    path = Path(tmp) / "legacy.tsv"
    write_rows(path, ["gene_id", "genome", "ko_id"],
               [{"gene_id": r["Unnamed: 0"], "genome": r["fasta"], "ko_id": r["ko_id"]} for r in rows])
    return path, rows


def legacy_features(rows):
    features = {}
    for row in rows:
        features.setdefault(row["fasta"], set()).update(parse_kos(row["ko_id"]))
    return features


def logic_parity(tmp, summary):
    path, rows = legacy_annotations(tmp)
    paths = specialise(path, Path(tmp) / "L1", genome_column="genome", ko_column="ko_id")
    _, ours = read_product(paths["product"])
    names, theirs = read_product(ANNOGUILD / "EMERGE_20260226_product_refined.tsv")
    mismatches = [{"genome": g, "pathway": n, "annoguild": theirs[g][n], "kolach": ours[g][n]}
                  for g in theirs for n in names if ours[g][n] != theirs[g][n]]
    summary["L1_pathways"] = {"genomes": len(theirs), "pathways": len(names),
                              "compared": len(theirs) * len(names), "mismatches": len(mismatches)}
    write_rows(HERE / "syncomm_L1_mismatches.tsv", ["genome", "pathway", "annoguild", "kolach"], mismatches)

    features = legacy_features(rows)
    _, with_cazy = read_product(ANNOGUILD / "EMERGE_20260226_product_refined_with_CAZy.tsv")
    expected = {r["genome"]: r for r in read_rows(ANNOGUILD / "EMERGE_20260226_genome_specialisations.tsv")}
    rows_out = []
    for genome, calls in with_cazy.items():
        labels = classify_specialisation(calls, features[genome])
        for column, value in labels.items():
            if value != expected[genome][column]:
                rows_out.append({"genome": genome, "label": column,
                                 "annoguild": expected[genome][column], "kolach": value})
    summary["L2_labels"] = {"genomes": len(with_cazy), "compared": len(with_cazy) * 3, "mismatches": len(rows_out)}
    write_rows(HERE / "syncomm_L2_mismatches.tsv", ["genome", "label", "annoguild", "kolach"], rows_out)
    return rows, features, with_cazy, expected


def end_to_end(legacy_rows, legacy_feature_sets, with_cazy, expected, summary):
    gene_genome = {r["Unnamed: 0"]: r["fasta"] for r in legacy_rows}
    annotations_by_gene = {r["gene_id"]: r for r in read_rows(KOLACH / "kolach_annotations.tsv")}
    legacy_ko = {r["Unnamed: 0"]: parse_kos(r["ko_id"]) for r in legacy_rows}

    # E1: GraftM. AnnoGuild marks KOs replaced by GraftM with database == graftm
    # (including mmoX, K16157). Compare genes refined by either side.
    theirs = {}
    for r in read_rows(ANNOGUILD / "all_annotations_filled_graftm.tsv"):
        if r["database"] == "graftm":
            theirs.setdefault(r["Unnamed: 0"], set()).update(
                LEGACY_MARKERS[k] for k in parse_kos(r["ko_id"]) if k in LEGACY_MARKERS)
    ours_rows = read_rows(KOLACH / "kolach_graftm_markers.tsv")
    ours = {}
    for row in ours_rows:
        if row["status"] == "refined":
            ours.setdefault(row["gene_id"], set()).add(row["marker"])
    marker_rows = []
    for gene in sorted(set(theirs) | set(ours)):
        detail = [r for r in ours_rows if r["gene_id"] == gene]
        a, k = theirs.get(gene, set()), ours.get(gene, set())
        marker_rows.append({
            "gene_id": gene, "genome": gene_genome.get(gene, ""),
            "annoguild_markers": ";".join(sorted(a)) or "-", "kolach_markers": ";".join(sorted(k)) or "-",
            "outcome": "same" if a == k else "annoguild_only" if not k else "kolach_only" if not a else "different",
            "kolach_status": ";".join(f"{r['package']}:{r['status']}:{r['placement'] or '-'}:{r['gate_source']}"
                                      for r in detail) or "not_gated",
            "annoguild_kos": ",".join(sorted(legacy_ko.get(gene, ()))),
            "kolach_accepted_ko": annotations_by_gene.get(gene, {}).get("accepted_ko", ""),
            "kolach_alternative_kos": annotations_by_gene.get(gene, {}).get("alternative_kos", ""),
        })
    outcomes = {}
    for row in marker_rows:
        outcomes[row["outcome"]] = outcomes.get(row["outcome"], 0) + 1
    summary["E1_graftm"] = {"genes_refined_by_either": len(marker_rows), "outcomes": outcomes,
                            "kolach_status_counts": {s: sum(r["status"] == s for r in ours_rows)
                                                     for s in ("refined", "excluded", "unresolved", "superseded", "not_placed")}}
    write_rows(HERE / "syncomm_E1_graftm.tsv",
               list(marker_rows[0]) if marker_rows else ["gene_id"], marker_rows)

    # E2: accepted KO vs AnnoGuild's KO set (both before marker refinement is
    # not recoverable for AnnoGuild, so compare on shared six-digit KOs only).
    annotations = annotations_by_gene
    counts = {"both_same": 0, "both_different": 0, "kolach_only": 0, "annoguild_only": 0, "neither": 0}
    for gene, kos in legacy_ko.items():
        six = {k[:6] if len(k) == 7 else k for k in kos}
        accepted = annotations[gene]["accepted_ko"]
        accepted = set() if accepted in ("", "-") else {accepted}
        if accepted and six:
            counts["both_same" if accepted <= six else "both_different"] += 1
        elif accepted:
            counts["kolach_only"] += 1
        elif six:
            counts["annoguild_only"] += 1
        else:
            counts["neither"] += 1
    summary["E2_ko"] = counts

    # E3: pathway calls and guild labels.
    names, ours_product = read_product(KOLACH / "kolach_product_refined.tsv")
    pathway_names = [n for n in names if not n.startswith("CAZy")]
    detail = {(r["genome"], r["pathway"]): r for r in read_rows(KOLACH / "kolach_metabolic_pathways.tsv")}
    diffs = []
    for genome in with_cazy:
        for name in pathway_names:
            if ours_product[genome][name] != with_cazy[genome][name]:
                row = detail[genome, name]
                diffs.append({"genome": genome, "pathway": name, "annoguild": with_cazy[genome][name],
                              "kolach": ours_product[genome][name],
                              "kolach_reaction_coverage": row["reaction_coverage"],
                              "kolach_ko_coverage": row["ko_coverage"],
                              "kolach_supporting": row["supporting_features"],
                              "annoguild_supporting": ";".join(sorted(
                                  set(row["supporting_features"].split(";") + row["missing_features"].split(";"))
                                  & legacy_feature_sets[genome] - {""}))})
    summary["E3_pathways"] = {"compared": len(with_cazy) * len(pathway_names), "differences": len(diffs)}
    write_rows(HERE / "syncomm_E3_pathways.tsv",
               ["genome", "pathway", "annoguild", "kolach", "kolach_reaction_coverage", "kolach_ko_coverage",
                "kolach_supporting", "annoguild_supporting"], diffs)
    ours_labels = {r["genome"]: r for r in read_rows(KOLACH / "kolach_specialisations.tsv")}
    cazy_theirs = {g: sorted(n for n, v in calls.items() if n.startswith("CAZy") and v) for g, calls in with_cazy.items()}
    cazy_ours = {g: sorted(n for n, v in calls.items() if n.startswith("CAZy") and v) for g, calls in ours_product.items()}
    labels = []
    for genome in sorted(expected):
        labels.append({"genome": genome,
                       **{f"annoguild_{c}": expected[genome][c] for c in ("specialisation", "nitrogen_specialisation", "sulfur_specialisation")},
                       **{f"kolach_{c}": ours_labels[genome][c] for c in ("specialisation", "nitrogen_specialisation", "sulfur_specialisation")},
                       "annoguild_cazy": ";".join(cazy_theirs[genome]), "kolach_cazy": ";".join(cazy_ours[genome])})
    write_rows(HERE / "syncomm_E3_labels.tsv", list(labels[0]), labels)
    summary["E3_labels"] = labels


def main():
    summary = {}
    with tempfile.TemporaryDirectory() as tmp:
        rows, features, with_cazy, expected = logic_parity(tmp, summary)
    if (KOLACH / "kolach_specialisations.tsv").is_file():
        end_to_end(rows, features, with_cazy, expected, summary)
    (HERE / "syncomm_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
    json.dump(summary, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
