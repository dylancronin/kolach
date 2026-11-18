"""Refine shared KOs with lineage-aware GraftM placement.

Several functions share one KO (pmoA/amoA, dsrA/rdsrA, narG/nxrA, ...). The
EMERGE GraftM packages place a protein in a curated tree; the clade it lands in
names its function. This module follows the temporal-paper workflow
(`GraftM-analysis/update_annotations.R`) rather than AnnoGuild's
`process_graftm.py`:

- Only proteins already carrying the package's KO are placed (accepted or
  alternative KO), so a package's HMM never classifies unrelated enzymes.
- A placement is interpreted by an exact (package, placement) lookup covering
  every clade in every package, never by substring matching. Substring matching
  turned every rdsrA/rdsrB into dsrA/dsrB, and every mmoX outgroup homologue into
  mmoX.
- Placements without a functional clade (the root, unknown clades, outgroups)
  are reported as unresolved and refine nothing.
"""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile

from kolach.download import download_atomic

# Package name -> KO whose proteins it is allowed to classify.
PACKAGES = {
    "dsrA": "K11180",
    "dsrB": "K11181",
    "narG_nxrA": "K00370",
    "narH_nxrB": "K00371",
    "amoA_pmoA": "K10944",
    "amoB_pmoB": "K10945",
    "amoC_pmoC": "K10946",
    "norZ_norB": "K04561",
    "mmoX_sMMO": "K16157",
}
# Refined feature IDs used by the EMERGE rules. Shared KO + 0/1 encodes the
# function; mmoX keeps its own, already specific KO.
MARKER_FEATURES = {
    "dsrA": "K111800", "rdsrA": "K111801", "dsrB": "K111810", "rdsrB": "K111811",
    "narG": "K003700", "nxrA": "K003701", "narH": "K003710", "nxrB": "K003711",
    "pmoA": "K109440", "amoA": "K109441", "pmoB": "K109450", "amoB": "K109451",
    "pmoC": "K109460", "amoC": "K109461", "norB": "K045610", "norZ": "K045611",
    "mmoX": "K16157",
}
BUNDLE_VERSION = "1"
BUNDLE_NAME = f"kolach-graftm-packages-v{BUNDLE_VERSION}.tar.gz"
BUNDLE_URL = f"https://github.com/dylancronin/kolach/releases/download/graftm-packages-v{BUNDLE_VERSION}/{BUNDLE_NAME}"
BUNDLE_SHA256 = "40c4136d6e91cf6c673638deae466aa54b3b4cb68605021610ee438c3c27cfec"
PLACEMENTS_FILE = Path(__file__).resolve().parents[1] / "metabolism" / "data" / "graftm_placements.tsv"
MARKER_COLUMNS = [
    "gene_id", "package", "gate_ko", "gate_source", "placement",
    "marker", "feature", "status", "unresolved_reason",
]
EXECUTABLES = ("graftM", "hmmsearch", "hmmalign", "pplacer", "guppy")


def check_dependencies():
    """Fail before launching the workflow when GraftM or its tools are missing."""
    for executable in EXECUTABLES:
        if not shutil.which(executable):
            raise FileNotFoundError(
                f"The '{executable}' executable was not found in PATH. "
                "Please activate the kolach conda environment with GraftM installed."
            )


def load_placements(path=PLACEMENTS_FILE) -> dict:
    """Return {(package, placement): (marker, unresolved_reason)}."""
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    table = {}
    for row in rows:
        key = (row["package"], row["placement"])
        marker, reason = row["marker"], row["unresolved_reason"]
        if key in table or row["package"] not in PACKAGES or bool(marker) == bool(reason):
            raise ValueError(f"Invalid GraftM placement row: {row}")
        if marker and marker not in MARKER_FEATURES:
            raise ValueError(f"Unknown marker {marker!r} in GraftM placement table")
        table[key] = (marker, reason)
    return table


def package_lineages(package_dir) -> set[str]:
    """Every placement string GraftM can report for a package.

    GraftM prints the tax_name chain from the root, joined by '; '.
    """
    files = list(Path(package_dir).glob("*.refpkg/*taxonomy*.csv"))
    if len(files) != 1:
        raise ValueError(f"Expected one taxonomy file in {package_dir}")
    with open(files[0], newline="", encoding="utf-8") as handle:
        rows = {row["tax_id"]: row for row in csv.DictReader(handle)}
    lineages = set()
    for row in rows.values():
        chain = []
        while True:
            chain.append(row["tax_name"])
            if row["tax_id"] == row["parent_id"]:
                break
            row = rows[row["parent_id"]]
        lineages.add("; ".join(reversed(chain)))
    return lineages


def validate_database(directory, placements=None):
    """Check every package is present and every clade it can report is mapped."""
    directory = Path(directory)
    placements = placements or load_placements()
    for package in PACKAGES:
        package_dir = directory / f"{package}.gpkg"
        if not (package_dir / "CONTENTS.json").is_file():
            raise FileNotFoundError(f"Missing GraftM package: {package_dir}")
        expected = {placement for name, placement in placements if name == package}
        found = package_lineages(package_dir)
        if found != expected:
            raise ValueError(
                f"GraftM package {package} clades do not match the placement table: "
                f"unmapped={sorted(found - expected)}, absent={sorted(expected - found)}")


def _safe_members(archive: tarfile.TarFile):
    for member in archive.getmembers():
        path = Path(member.name)
        if path.is_absolute() or ".." in path.parts or not (member.isfile() or member.isdir()):
            raise ValueError(f"Unsafe entry in GraftM package bundle: {member.name}")
        yield member


def install_bundle(archive_path, directory, url=None):
    """Extract a verified bundle into `directory`, publishing only valid packages."""
    directory = Path(directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    (directory / ".download_complete").unlink(missing_ok=True)
    digest = hashlib.sha256(Path(archive_path).read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory(prefix="install-", dir=directory) as tmp:
        stage = Path(tmp)
        with tarfile.open(archive_path, "r:gz") as archive:
            members = list(_safe_members(archive))
            # The "data" filter (Python 3.11.4+) also rejects unsafe paths and modes.
            options = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}
            archive.extractall(stage, members=members, **options)
        validate_database(stage)
        manifest = {"bundle_version": BUNDLE_VERSION, "url": url, "sha256": digest,
                    "packages": sorted(f"{name}.gpkg" for name in PACKAGES)}
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for name in PACKAGES:
            target = directory / f"{name}.gpkg"
            if target.exists():
                shutil.rmtree(target)
            (stage / f"{name}.gpkg").replace(target)
        (stage / "manifest.json").replace(directory / "manifest.json")
    (directory / ".download_complete").touch()


def download_graftm(directory, url=BUNDLE_URL, sha256=BUNDLE_SHA256):
    """Download the pinned package bundle, check its hash, then install it."""
    directory = Path(directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="download-", dir=directory) as tmp:
        archive = Path(tmp) / BUNDLE_NAME
        digest = download_atomic(url, archive)
        if digest != sha256:
            raise ValueError(f"GraftM bundle checksum mismatch: expected {sha256}, got {digest}")
        install_bundle(archive, directory, url=url)


def build_bundle(package_dirs: dict, output) -> str:
    """Write a reproducible bundle (fixed order, owners and times); return its SHA256."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.PAX_FORMAT) as archive:
        for name in sorted(package_dirs):
            root = Path(package_dirs[name])
            for path in sorted([root, *root.rglob("*")]):
                if path.is_symlink():
                    raise ValueError(f"Symlinks are not bundled: {path}")
                info = archive.gettarinfo(path, arcname=str(Path(f"{name}.gpkg") / path.relative_to(root)))
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                info.mtime = 0
                info.mode = 0o755 if path.is_dir() else 0o644
                if path.is_file():
                    with open(path, "rb") as handle:
                        archive.addfile(info, handle)
                else:
                    archive.addfile(info)
    data = gzip.compress(buffer.getvalue(), mtime=0)
    Path(output).write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def _kos(value):
    return {token for token in value.replace(",", ";").split(";") if token.strip() not in ("", "-")}


def gate_genes(annotation_table) -> dict:
    """Return {package: [(gene_id, gate_source), ...]} in table order."""
    with open(annotation_table, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"gene_id", "accepted_ko", "alternative_kos"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Annotation table lacks columns needed for GraftM gating: {sorted(missing)}")
        rows = list(reader)
    gated = {package: [] for package in PACKAGES}
    for row in rows:
        accepted, alternatives = _kos(row["accepted_ko"]), _kos(row["alternative_kos"])
        for package, ko in PACKAGES.items():
            if ko in accepted:
                gated[package].append((row["gene_id"], "accepted"))
            elif ko in alternatives:
                gated[package].append((row["gene_id"], "alternative"))
    return gated


def read_fasta(path, wanted) -> dict:
    """Sequences for `wanted` IDs (first header token, as in integration)."""
    from kolach.integrate import _open_text

    sequences, name, chunks = {}, None, []
    with _open_text(Path(path)) as handle:
        for line in handle:
            if line.startswith(">"):
                if name in wanted:
                    sequences[name] = "".join(chunks)
                parts = line[1:].split()
                name, chunks = (parts[0] if parts else None), []
            elif name in wanted:
                chunks.append(line.strip())
    if name in wanted:
        sequences[name] = "".join(chunks)
    return sequences


def read_placements(path) -> dict:
    """Parse GraftM's <name>_read_tax.tsv: sequence ID and placement string."""
    placements = {}
    if not Path(path).is_file():
        return placements
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                sequence_id, placement = line.rstrip("\n").split("\t", 1)
                placements[sequence_id] = placement.strip()
    return placements


def interpret(package, placement, table) -> dict:
    """Map one placement to a marker row (status refined or unresolved)."""
    if (package, placement) not in table:
        return {"placement": placement, "marker": "", "feature": "", "status": "unresolved",
                "unresolved_reason": "placement_not_in_lookup"}
    marker, reason = table[package, placement]
    if marker:
        return {"placement": placement, "marker": marker, "feature": MARKER_FEATURES[marker],
                "status": "refined", "unresolved_reason": ""}
    return {"placement": placement, "marker": "", "feature": "", "status": "unresolved",
            "unresolved_reason": reason}


def annotate(protein_fasta, annotation_table, package_dir, output, threads=1, workdir=None):
    """Place gated proteins with every package and write the marker table.

    Every gated protein gets one row per package: refined, unresolved, or
    not_placed (GraftM's own search did not keep it).
    """
    table = load_placements()
    package_dir = Path(package_dir)
    validate_database(package_dir, table)
    gated = gate_genes(annotation_table)
    sequences = read_fasta(protein_fasta, {gene for genes in gated.values() for gene, _ in genes})
    output = Path(output)
    workdir = Path(workdir or output.parent / "graftm")
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    rows = []
    for package, genes in gated.items():
        if not genes:
            continue
        # Short positional IDs: GraftM and pplacer never see the original names.
        ids = {f"q{index}": gene for index, (gene, _) in enumerate(genes)}
        query = workdir / f"{package}.faa"
        with open(query, "w", encoding="utf-8") as handle:
            for safe, gene in ids.items():
                if gene not in sequences:
                    raise ValueError(f"Gene {gene!r} is in the annotation table but not the FASTA")
                handle.write(f">{safe}\n{sequences[gene]}\n")
        out_dir = workdir / package
        log = workdir / f"{package}.log"
        with open(log, "w") as handle:
            subprocess.run(["graftM", "graft", "--forward", str(query),
                            "--graftm_package", str(package_dir / f"{package}.gpkg"),
                            "--output_directory", str(out_dir), "--threads", str(threads), "--force"],
                           stdout=handle, stderr=subprocess.STDOUT, check=True)
        placed = read_placements(out_dir / package / f"{package}_read_tax.tsv")
        unknown = set(placed) - set(ids)
        if unknown:
            raise ValueError(f"GraftM reported unexpected sequence IDs for {package}: {sorted(unknown)[:5]}")
        for safe, (gene, source) in zip(ids, genes):
            row = {"gene_id": gene, "package": package, "gate_ko": PACKAGES[package], "gate_source": source}
            if safe in placed:
                row.update(interpret(package, placed[safe], table))
            else:
                row.update({"placement": "", "marker": "", "feature": "", "status": "not_placed",
                            "unresolved_reason": "no_graftm_hit"})
            rows.append(row)
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MARKER_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    shutil.rmtree(workdir)
    return output
