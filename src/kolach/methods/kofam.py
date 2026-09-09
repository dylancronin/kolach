"""KOfam search and annotation following anvi'o v9's default KO decision rules.

Reference: merenlab/anvio, tag v9, anvio/metabolism/annotate.py.
This implementation operates on protein FASTA files, without a contigs database.
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen

REFERENCE = "https://github.com/merenlab/anvio/blob/v9/anvio/metabolism/annotate.py"
BASE_URL = "https://www.genome.jp/ftp/db/kofam/"
FIELDS = ["gene_id", "ko", "assignment", "score_type", "threshold", "bit_score",
          "e_value", "domain_bit_score", "domain_e_value", "definition"]


def read_profiles(path):
    profiles = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["threshold"] == "-":
                continue  # anvi'o's default excludes KOs without thresholds
            threshold = float(row["threshold"])
            if not math.isfinite(threshold) or row["score_type"] not in {"full", "domain"}:
                raise ValueError(f"Invalid threshold/score type for {row['knum']}")
            profiles[row["knum"]] = {**row, "threshold": threshold}
    if not profiles:
        raise ValueError("No KOfam profiles with numeric thresholds found")
    return profiles


def download_kofam(directory):
    """Download, validate and combine threshold-bearing profiles; retain checksums."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        tmp = Path(tmp)
        manifest = {"method_reference": REFERENCE, "downloads": {}}
        for name in ("ko_list.gz", "profiles.tar.gz"):
            digest = hashlib.sha256()
            with urlopen(BASE_URL + name, timeout=120) as response, open(tmp / name, "wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
                    digest.update(block)
            manifest["downloads"][name] = {"url": BASE_URL + name, "sha256": digest.hexdigest()}
        with gzip.open(tmp / "ko_list.gz", "rb") as source, open(tmp / "ko_list", "wb") as out:
            shutil.copyfileobj(source, out)
        profiles = read_profiles(tmp / "ko_list")
        found = set()
        with tarfile.open(tmp / "profiles.tar.gz", "r:gz") as archive, open(tmp / "profiles.hmm", "wb") as out:
            # Stream regular files only; do not extract archive paths or links.
            for member in archive:
                ko = Path(member.name).stem
                if member.isfile() and member.name.endswith(".hmm") and ko in profiles:
                    if ko in found:
                        raise ValueError(f"Duplicate profile: {ko}")
                    with archive.extractfile(member) as source:
                        shutil.copyfileobj(source, out)
                    found.add(ko)
        if found != set(profiles):
            raise ValueError(f"Missing {len(set(profiles) - found)} KOfam profiles")
        (tmp / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        for name in ("ko_list", "profiles.hmm", "manifest.json"):
            (tmp / name).replace(directory / name)


def count_proteins(path):
    ids = set()
    has_sequence = False
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if ids and not has_sequence:
                    raise ValueError("FASTA contains an empty protein sequence")
                header = line[1:].split()
                if not header or header[0] in ids:
                    raise ValueError("FASTA headers must have unique, nonempty IDs")
                ids.add(header[0])
                has_sequence = False
            else:
                if not ids:
                    raise ValueError("FASTA sequence appears before its header")
                has_sequence = True
    if not ids or not has_sequence:
        raise ValueError("FASTA contains no proteins or an empty final sequence")
    return len(ids)


def parse_hits(path):
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split(maxsplit=18)
            if len(cols) < 18:
                raise ValueError("Malformed hmmsearch --tblout row")
            yield {"gene_id": cols[0], "ko": cols[2],
                   "e_value": float(cols[4]), "bit_score": float(cols[5]),
                   "domain_e_value": float(cols[7]), "domain_bit_score": float(cols[8])}


def assign_hits(hits, profiles, *, fraction=0.75, evalue=1e-5, rescue=True):
    """Apply initial >= thresholds, then unique-KO rescue with strict > scores."""
    if not math.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError("Heuristic bit-score fraction must be between 0 and 1")
    if not math.isfinite(evalue) or evalue < 0:
        raise ValueError("Heuristic E-value must be finite and nonnegative")
    accepted, candidates, annotated = [], {}, set()
    for hit in hits:
        profile = profiles[hit["ko"]]
        domain = profile["score_type"] == "domain"
        score = hit["domain_bit_score" if domain else "bit_score"]
        ev = hit["domain_e_value" if domain else "e_value"]
        row = {**hit, "score_type": profile["score_type"],
               "threshold": profile["threshold"], "definition": profile.get("definition", "")}
        if score >= profile["threshold"]:
            accepted.append({**row, "assignment": "threshold"})
            annotated.add(hit["gene_id"])
        elif rescue and ev <= evalue and score > fraction * profile["threshold"]:
            gene = hit["gene_id"]
            if gene not in candidates:
                candidates[gene] = [set(), float("inf"), None]
            candidate = candidates[gene]
            candidate[0].add(hit["ko"])
            if ev <= candidate[1]:
                candidate[1:] = [ev, row]
    for gene, (kos, _, row) in candidates.items():
        if gene not in annotated and len(kos) == 1:
            accepted.append({**row, "assignment": "rescued"})
    return sorted(accepted, key=lambda row: (row["gene_id"], row["ko"]))


def annotate(args, output):
    db = Path(args.database_dir).expanduser().resolve() / "kofam"
    profiles = read_profiles(db / "ko_list")
    fasta = Path(args.protein_fasta).expanduser().resolve()
    count = count_proteins(fasta)
    threads = getattr(args, "threads", 1)
    if threads < 1:
        raise ValueError("Threads must be at least 1")
    with tempfile.TemporaryDirectory(prefix="kolach-kofam-") as tmp:
        table = Path(tmp) / "hits.tbl"
        cmd = ["hmmsearch", "--cpu", str(threads), "--noali", "-o", "/dev/null",
               "--tblout", str(table), "-Z", str(count), "--domZ", str(count)]
        if getattr(args, "no_hmmer_prefiltering", False):
            cmd.extend(["-T", "-20", "--domT", "-20"])
        cmd.extend([str(db / "profiles.hmm"), str(fasta)])
        subprocess.run(cmd, check=True)
        rows = assign_hits(parse_hits(table), profiles,
                           fraction=getattr(args, "heuristic_bitscore_fraction", 0.75),
                           evalue=getattr(args, "heuristic_e_value", 1e-5),
                           rescue=not getattr(args, "skip_bitscore_heuristic", False))
        writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
