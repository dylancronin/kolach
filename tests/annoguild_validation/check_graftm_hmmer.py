"""Check GraftM 0.15.1 placements with HMMER 3.4 against HMMER 3.2.1 results.

A. Temporal-paper test sets: each package's EMERGE sequences (in combined.faa)
   versus the recorded GraftM-analysis/<pkg>/emerge_read_tax.tsv.
B. SynComm: AnnoGuild's all_genes.faa (all proteins, every package) versus its
   recorded graftm_results/<pkg>_graftm/all_genes/all_genes_read_tax.tsv.

Writes graftm_hmmer_check.json and graftm_hmmer_mismatches.tsv here.
"""
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
TEMPORAL = Path("/fs/scratch/PDS0325/Cronin/temporal-paper/temporal_git/temporal_paper_202505/GraftM-analysis")
SYNCOMM = Path("/fs/scratch/PDS0325/Cronin/syncomm/annoguild/graftm_processing")
PACKAGES = HERE / "db" / "graftm"
THREADS = sys.argv[1] if len(sys.argv) > 1 else "4"
# Optional label, e.g. "shared_env", to keep results from another environment apart.
LABEL = f"_{sys.argv[2]}" if len(sys.argv) > 2 else ""
WORK = HERE / f"graftm_hmmer_work{LABEL}"


def read_tax(path):
    rows = {}
    if path.is_file():
        for line in path.read_text().splitlines():
            if line.strip():
                name, placement = line.split("\t", 1)
                rows[name] = placement.strip()
    return rows


def graft(query, package, out):
    subprocess.run(["graftM", "graft", "--forward", str(query), "--graftm_package", str(PACKAGES / f"{package}.gpkg"),
                    "--output_directory", str(out), "--threads", THREADS, "--force"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return read_tax(out / query.stem / f"{query.stem}_read_tax.tsv")


def subset_fasta(source, names, target):
    keep, out = False, []
    for line in source.read_text().splitlines():
        if line.startswith(">"):
            keep = line[1:].split()[0] in names
        if keep:
            out.append(line)
    target.write_text("\n".join(out) + "\n")


def compare(label, package, expected, observed, mismatches):
    names = set(expected) | set(observed)
    same = sum(expected.get(n) == observed.get(n) for n in names)
    for n in sorted(names):
        if expected.get(n) != observed.get(n):
            mismatches.append([label, package, n, expected.get(n, "-"), observed.get(n, "-")])
    return {"expected": len(expected), "observed": len(observed), "identical": same, "compared": len(names)}


def main():
    WORK.mkdir(exist_ok=True)
    summary, mismatches = {"A_temporal": {}, "B_syncomm": {}}, []
    for package in ["amoA_pmoA", "amoB_pmoB", "amoC_pmoC", "dsrA", "dsrB", "narG_nxrA", "narH_nxrB", "norZ_norB"]:
        expected = read_tax(TEMPORAL / package / "emerge_read_tax.tsv")
        query = WORK / f"temporal_{package}.faa"
        subset_fasta(TEMPORAL / package / "combined.faa", set(expected), query)
        observed = graft(query, package, WORK / f"A_{package}")
        summary["A_temporal"][package] = compare("A", package, expected, observed, mismatches)
    for package in ["amoA_pmoA", "amoB_pmoB", "amoC_pmoC", "dsrA", "dsrB", "narG_nxrA", "narH_nxrB",
                    "norZ_norB", "mmoX_sMMO"]:
        expected = read_tax(SYNCOMM / "graftm_results" / f"{package}_graftm" / "all_genes" / "all_genes_read_tax.tsv")
        observed = graft(SYNCOMM / "all_genes.faa", package, WORK / f"B_{package}")
        summary["B_syncomm"][package] = compare("B", package, expected, observed, mismatches)
    versions = subprocess.run(["hmmsearch", "-h"], capture_output=True, text=True).stdout.splitlines()[1]
    summary["hmmer"] = versions.strip("# ").strip()
    summary["pplacer"] = subprocess.run(["pplacer", "--version"], capture_output=True, text=True).stdout.strip()
    summary["environment"] = sys.prefix
    (HERE / f"graftm_hmmer_check{LABEL}.json").write_text(json.dumps(summary, indent=2) + "\n")
    with open(HERE / f"graftm_hmmer_mismatches{LABEL}.tsv", "w") as handle:
        handle.write("check\tpackage\tsequence\thmmer_3.2.1\thmmer_3.4\n")
        for row in mismatches:
            handle.write("\t".join(row) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
