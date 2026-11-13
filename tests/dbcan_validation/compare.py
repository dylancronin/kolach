"""Compare actual DRAM2 searches, normalized evidence, and its original combiner."""
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
import sys

from kolach.integrate import integrate_annotations, read_fasta_ids
from kolach.methods.dbcan import sha256

ROOT = Path(__file__).resolve().parent


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def canonical(rows):
    return Counter(json.dumps(row, sort_keys=True, separators=(",", ":")) for row in rows)


def main():
    kolach = ROOT / "kolach_output"
    reference = ROOT / "reference_output"
    native = kolach / "dbcan"
    mapping = {r["search_id"]: r["gene_id"] for r in read(native / "gene_ids.tsv")}
    report = {"reference_commit": "4536cd93c97c07d783fab9e9ef99e56b271f3782", "native": {}}
    evidence = read(kolach / "dbcan_evidence.tsv")
    for name, method in [("dbCAN_hmm_results.tsv", "hmm"),
                         ("dbCANsub_hmm_results.tsv", "dbCANsub"), ("diamond.out", "diamond")]:
        ours = read(native / name)
        theirs = read(reference / f"genome_{name}")
        expected, actual = canonical(theirs), canonical(ours)
        missing, additional = expected - actual, actual - expected
        report["native"][method] = {"kolach_rows": len(ours), "dram2_rows": len(theirs),
                                     "missing": sum(missing.values()), "additional": sum(additional.values())}
        assert not missing and not additional, (method, missing, additional)
        # Every reference domain/alignment must survive in Kolach's audit table.
        audited = Counter(r["source_row"] for r in evidence if r["method"] == method)
        assert audited == actual, (method, "lost normalized evidence")
    assert canonical(read(native / "overview.tsv")) == canonical(read(reference / "genome_overview.tsv"))
    report["upstream_overview_equal"] = True
    summary_path = kolach / "dbcan_overview.tsv"
    if not summary_path.exists():
        summary_path = kolach / "dbcan_annotations.tsv"  # Original validation artifacts.
    summary = {r["gene_id"]: r for r in read(summary_path)}
    for row in read(reference / "genome_overview.tsv"):
        gene = mapping[row["Gene ID"]]
        assert summary[gene]["dbcan_recommended"] == row["Recommend Results"]
        assert summary[gene]["dbcan_tools_count"] == row["#ofTools"]
    fasta = ROOT.parent / "anvio_validation" / "proteins.faa"
    ids = read_fasta_ids(fasta)
    assert list(summary) == ids
    report["gene_inventory"] = len(ids)
    report["genes_with_hits"] = sum(r["dbcan_tools_count"] != "0" for r in summary.values())
    report["genes_with_recommendations"] = sum(r["dbcan_recommended"] != "-" for r in summary.values())

    # Invoke DRAM2's real table combiner, without gene-location extraction. Its
    # summary is compared separately because first-value aggregation is lossy.
    combined = ROOT / "dram2_combined.tsv"
    subprocess.run([sys.executable, str(ROOT / "upstream" / "combine_annotations.py"),
                    "--dbcan_dir", str(reference), "--output", str(combined)], check=True)
    extra = []
    for row in read(combined):
        gene = mapping[row["query_id"]]
        for column in ("dbcan_id", "dbcan_sub_id", "dbcan_sub_composition", "dbcan_sub_ec", "dbcan_sub_substrate"):
            retained = {part for part in row.get(column, "").split(";") if part not in ("", "-")}
            ours = {part for part in summary[gene][column].split(";") if part not in ("", "-")}
            assert retained <= ours, (gene, column, retained - ours)
            for value in sorted(ours - retained):
                extra.append({"gene_id": gene, "column": column, "additional_retained_value": value})
    from kolach.methods.dbcan import write_table
    write_table(ROOT / "additional_domains.tsv", ["gene_id", "column", "additional_retained_value"], extra)
    report["additional_summary_values_vs_dram2"] = len(extra)

    # Reuse all completed KO searches; compare every existing column and audit row.
    validation = ROOT.parent / "anvio_validation"
    kwargs = dict(protein_fasta=fasta, kofam_tsv=validation / "matched_output/kofam_annotations.tsv",
                  eggnog_tsv=validation / "kolach_output/eggnog_annotations.tsv",
                  deepkoala_tsv=validation / "kolach_output/deepkoala_annotations.tsv")
    baseline = integrate_annotations(**kwargs)
    enriched = integrate_annotations(**kwargs, dbcan_tsv=summary_path,
                                     dbcan_evidence_tsv=kolach / "dbcan_evidence.tsv",
                                     output_tsv=kolach / "kolach_annotations_with_ko.tsv")
    import pandas as pd
    pd.testing.assert_frame_equal(baseline, enriched[baseline.columns])
    combined_evidence = enriched.attrs["evidence"]
    ko_evidence = combined_evidence[combined_evidence.method != "dbcan"][baseline.attrs["evidence"].columns]
    pd.testing.assert_frame_equal(baseline.attrs["evidence"], ko_evidence.reset_index(drop=True), check_dtype=False)
    domains = combined_evidence[combined_evidence.method == "dbcan"]
    assert canonical(domains["source_row"].tolist()) == canonical([row["source_row"] for row in evidence])
    import re
    pairs = {(r.gene_id, ko) for r in enriched.itertuples()
             for ko in re.findall(r"K\d{5}", r.accepted_ko + ";" + r.alternative_kos)}
    reference_pairs = {(r["gene_id"], r["ko"]) for r in read(validation / "anvio_gene_ko.tsv")}
    assert reference_pairs <= pairs
    report["ko_columns_and_evidence_unchanged"] = True
    report["retained_anvio_pairs"] = len(reference_pairs)
    report["versions_and_database"] = json.loads((kolach / "dbcan_metadata.json").read_text())
    report["reference_sources"] = json.loads((ROOT / "reference_sources.json").read_text())
    report["protein_fasta_sha256"] = sha256(fasta)
    (ROOT / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ("versions_and_database", "reference_sources")}, indent=2))


if __name__ == "__main__":
    main()
