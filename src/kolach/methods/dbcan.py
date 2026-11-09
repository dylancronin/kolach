"""Run pinned dbCAN searches and retain domain-level, non-KO evidence.

Reference: DRAM2 4536cd93, modules/nf-core/rundbcan/cazymeannotation.
Only upstream-filtered hits are imported; recommendations are kept separately.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

VERSION = "5.2.9"
RELEASE = "db_v5-2-9_5-5-2026"
BASE_URL = f"https://dbcan.s3.us-west-2.amazonaws.com/{RELEASE}"
ASSETS = {
    "CAZy.dmnd": "CAZy.dmnd", "dbCAN.hmm": "dbCAN.hmm",
    "dbCAN-sub.hmm": "dbCAN_sub.hmm",
    "fam-substrate-mapping.tsv": "fam-substrate-mapping.tsv",
}
SUMMARY_COLUMNS = [
    "dbcan_id", "dbcan_sub_id", "dbcan_sub_composition", "dbcan_sub_ec",
    "dbcan_sub_substrate", "dbcan_family_ec", "dbcan_cazy_id",
    "dbcan_methods", "dbcan_tools_count", "dbcan_recommended",
]
EVIDENCE_COLUMNS = [
    "gene_id", "method", "feature_id", "query_start", "query_end",
    "e_value", "coverage", "bit_score", "ec", "substrate", "source_row",
]
FAMILY = re.compile(r"(?:GH|GT|CBM|AA|CE|PL)\d+(?:_\d+)?")
EC = re.compile(r"\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)")
OWNER_FILE = ".kolach_dbcan_native"
OWNER_TEXT = "Kolach dbCAN native outputs v1\n"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _run(command, log):
    with open(log, "w") as handle:
        subprocess.run([str(x) for x in command], stdout=handle,
                       stderr=subprocess.STDOUT, check=True)
    # Some upstream stages catch errors and still return exit status zero.
    if re.search(r"(?:^|\s)(?:ERROR|CRITICAL)(?::|\s)", Path(log).read_text(), re.M):
        raise RuntimeError(f"dbCAN reported an error; see {log}")


def tool_version():
    result = subprocess.run(["run_dbcan", "version"], check=True,
                            capture_output=True, text=True)
    if result.stdout.strip() != f"dbCAN version: {VERSION}":
        raise RuntimeError(f"Expected dbCAN {VERSION}, got {result.stdout.strip()!r}")
    return VERSION


def _read_table(path, required):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"Missing required columns in {path}: {required}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Duplicate columns in {path}")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"Malformed TSV row in {path}")
        return rows


def write_table(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_database(directory):
    directory = Path(directory)
    for name in ASSETS:
        path = directory / name
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(f"Missing or empty dbCAN database asset: {path}")
    for name in ("dbCAN.hmm", "dbCAN-sub.hmm"):
        subprocess.run(["hmmstat", str(directory / name)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    subprocess.run(["diamond", "dbinfo", "--db", str(directory / "CAZy.dmnd")],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    with open(directory / "fam-substrate-mapping.tsv") as handle:
        header = [column.strip() for column in handle.readline().rstrip().split("\t")]
        first = handle.readline().rstrip().split("\t")
        if not {"Family", "EC_Number", "Substrate_curated"} <= set(header) or len(first) != len(header):
            raise ValueError("Invalid dbCAN family/substrate mapping table")


def download_dbcan(directory):
    """Use the upstream fixed-release downloader; publish only validated assets."""
    directory = Path(directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    version = tool_version()
    # Invalidate any old success marker before attempting an update.
    (directory / ".download_complete").unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix="download-", dir=directory) as tmp:
        stage = Path(tmp)
        command = ["run_dbcan", "database", "--db_dir", str(stage), "--aws_s3", "--no-cgc"]
        _run(command, directory / "download.log")
        validate_database(stage)
        manifest = {"tool_version": version, "release": RELEASE, "command": command,
                    "files": {name: {"url": f"{BASE_URL}/{remote}",
                                     "sha256": sha256(stage / name)}
                              for name, remote in ASSETS.items()}}
        for name in ASSETS:
            (stage / name).replace(directory / name)
        (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (directory / ".download_complete").touch()


def _parts(value, separator=";"):
    return {part.strip() for part in value.split(separator) if part.strip() not in ("", "-")}


def _number(value, label, minimum=0, maximum=None, integer=False):
    try:
        number = float(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"Invalid {label}: {value!r}") from error
    if (not math.isfinite(number) or number < minimum or
            maximum is not None and number > maximum or integer and not number.is_integer()):
        raise ValueError(f"Invalid {label}: {value!r}")
    return number


def parse_outputs(native, gene_ids, id_map=None):
    """Return stable summaries and all native records, without dropping domains."""
    native = Path(native)
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError("Duplicate gene IDs in dbCAN protein inventory")
    id_map = id_map or {}
    genes = {gene: {col: set() for col in SUMMARY_COLUMNS} for gene in gene_ids}
    evidence, methods = [], {gene: set() for gene in gene_ids}
    coordinates = ["Target From", "Target To", "i-Evalue", "Coverage"]
    sources = [
        ("hmm", "dbCAN_hmm_results.tsv", "HMM Name", "dbcan_id"),
        ("dbCANsub", "dbCANsub_hmm_results.tsv", "Subfam Name", "dbcan_sub_id"),
        ("diamond", "diamond.out", "CAZy ID", "dbcan_cazy_id"),
    ]
    for method, filename, feature_col, summary_col in sources:
        required = (["Gene ID", feature_col, "Gene Start", "Gene End", "E Value", "Bit Score"]
                    if method == "diamond" else ["Target Name", feature_col, *coordinates])
        if method == "dbCANsub":
            required += ["Subfam Composition", "Subfam EC", "Substrate"]
        for row in _read_table(native / filename, required):
            source_gene = row["Gene ID" if method == "diamond" else "Target Name"]
            gene = id_map.get(source_gene, source_gene)
            if gene not in genes:
                raise ValueError(f"dbCAN returned unknown gene ID: {gene!r}")
            feature = row[feature_col].removesuffix(".hmm")
            if not feature or feature == "-":
                raise ValueError(f"Missing dbCAN feature for {gene!r}")
            start, end = ((row["Gene Start"], row["Gene End"]) if method == "diamond"
                          else (row["Target From"], row["Target To"]))
            _number(start, "query start", minimum=1, integer=True)
            _number(end, "query end", minimum=float(start), integer=True)
            evalue = row["E Value" if method == "diamond" else "i-Evalue"]
            _number(evalue, "E-value")
            coverage, score = row.get("Coverage", "-"), row.get("Bit Score", "-")
            if method != "diamond":
                _number(coverage, "coverage", maximum=1)
            else:
                _number(score, "bit score")
            summary = genes[gene]
            methods[gene].add(method)
            if method == "diamond":
                summary[summary_col].update(part for part in feature.split("|") if FAMILY.fullmatch(part))
            else:
                summary[summary_col].update(_parts(feature))
            if method == "dbCANsub":
                for source, dest in [("Subfam Composition", "dbcan_sub_composition"),
                                     ("Subfam EC", "dbcan_sub_ec"), ("Substrate", "dbcan_sub_substrate")]:
                    summary[dest].update(_parts(row[source]))
                # The family is encoded in the domain's subfamily name; its
                # composition field lists source-sequence counts, not families.
                families = {match.split("_")[0] for match in FAMILY.findall(feature)}
                ecs = {part.split(":")[0] for part in _parts(row["Subfam EC"])}
                if any(not EC.fullmatch(ec) for ec in ecs):
                    raise ValueError(f"Malformed subfamily EC for {gene}: {row['Subfam EC']}")
                # Ambiguous multi-family domains do not establish correlated pairs.
                if len(families) == 1:
                    summary["dbcan_family_ec"].update(f"{next(iter(families))};EC:{ec}" for ec in ecs)
            evidence.append(dict(zip(EVIDENCE_COLUMNS, [
                gene, method, feature, start, end, evalue, coverage, score,
                row.get("Subfam EC", "-"), row.get("Substrate", "-"),
                json.dumps(row, sort_keys=True, separators=(",", ":")),
            ])))
    overview = _read_table(native / "overview.tsv", ["Gene ID", "#ofTools", "Recommend Results"])
    seen = set()
    for row in overview:
        gene = id_map.get(row["Gene ID"], row["Gene ID"])
        if gene not in genes or gene in seen:
            raise ValueError(f"Unknown or duplicate overview gene ID: {gene!r}")
        seen.add(gene)
        count = _number(row["#ofTools"], "tool count", maximum=3, integer=True)
        if count != len(methods[gene]):
            raise ValueError(f"Inconsistent overview method count for {gene}")
        genes[gene]["dbcan_recommended"] = row["Recommend Results"] or "-"
    if any(methods[gene] and gene not in seen for gene in genes):
        raise ValueError("dbCAN overview omits genes with search hits")
    summaries = []
    for gene, values in genes.items():
        row = {"gene_id": gene}
        for col in SUMMARY_COLUMNS:
            value = values[col]
            row[col] = value if isinstance(value, str) else (
                ("|" if col == "dbcan_family_ec" else ";").join(sorted(value)) or "-")
        row["dbcan_recommended"] = values["dbcan_recommended"] if gene in seen else "-"
        row["dbcan_methods"] = ";".join(sorted(methods[gene])) or "-"
        row["dbcan_tools_count"] = str(len(methods[gene]))
        summaries.append(row)
    evidence.sort(key=lambda row: (row["gene_id"], row["method"], int(row["query_start"]),
                                   row["feature_id"], row["source_row"]))
    return summaries, evidence


def prepare_input(fasta, directory, gene_ids):
    """Protect IDs from upstream pandas type inference without changing sequences."""
    from kolach.integrate import _open_text

    id_map = {f"kolach_dbcan_{i:08d}": gene for i, gene in enumerate(gene_ids)}
    forward = {gene: safe for safe, gene in id_map.items()}
    directory = Path(directory)
    with _open_text(Path(fasta)) as source, open(directory / "input.faa", "w") as target:
        for line in source:
            if line.startswith(">"):
                target.write(f">{forward[line[1:].split()[0]]}\n")
            else:
                target.write(line)
    write_table(directory / "gene_ids.tsv", ["search_id", "gene_id"],
                [{"search_id": safe, "gene_id": gene} for safe, gene in id_map.items()])
    return id_map


def read_summary(path):
    rows = _read_table(path, ["gene_id", *SUMMARY_COLUMNS])
    seen = set()
    for row in rows:
        if not row["gene_id"] or row["gene_id"] in seen:
            raise ValueError(f"Empty or duplicate dbCAN gene ID in {path}")
        seen.add(row["gene_id"])
        _number(row["dbcan_tools_count"], "tool count", maximum=3, integer=True)
    return [{key: row[key] for key in ["gene_id", *SUMMARY_COLUMNS]} for row in rows]


def _check_output_paths(database, fasta, output, published):
    """Refuse collisions before running searches or replacing any output."""
    if published.is_symlink():
        raise ValueError(f"dbCAN native output directory must not be a symlink: {published}")
    destination = published.resolve()
    if destination.is_relative_to(database) or database.is_relative_to(destination):
        raise ValueError("dbCAN database and native output directories must not overlap")
    products = [output, output.parent / "dbcan_evidence.tsv", output.parent / "dbcan_metadata.json"]
    if fasta.is_relative_to(destination) or any(fasta == path.resolve() for path in products):
        raise ValueError("dbCAN output paths must not contain or overwrite the input FASTA")
    protected = [fasta, *(database / name for name in ASSETS), database / "manifest.json"]
    for path in products:
        if path.resolve().is_relative_to(database) or any(
                path.exists() and source.exists() and path.samefile(source) for source in protected):
            raise ValueError("dbCAN output files must not alias the database or input FASTA")
    if published.exists():
        if not published.is_dir():
            raise ValueError(f"dbCAN native output path is not a directory: {published}")
        marker = published / OWNER_FILE
        if any(published.iterdir()) and (not marker.is_file() or marker.read_text() != OWNER_TEXT):
            raise ValueError(f"Refusing to overwrite a nonempty directory not owned by Kolach: {published}")


def annotate(args, output):
    from kolach.integrate import read_fasta_ids

    fasta = Path(args.protein_fasta).expanduser().resolve()
    database = (Path(args.database_dir).expanduser().resolve() / "dbcan").resolve()
    output = Path(output).expanduser().resolve()
    published = output.parent / "dbcan"
    _check_output_paths(database, fasta, output, published)
    gene_ids = read_fasta_ids(fasta)
    if not gene_ids or len(gene_ids) != len(set(gene_ids)):
        raise ValueError("dbCAN requires a nonempty FASTA with unique gene IDs")
    threads = int(args.threads)
    if threads < 1:
        raise ValueError("dbCAN threads must be positive")
    validate_database(database)
    version = tool_version()
    output.parent.mkdir(parents=True, exist_ok=True)
    native = Path(tempfile.mkdtemp(prefix="dbcan-run-", dir=output.parent))
    id_map = prepare_input(fasta, native, gene_ids)
    command = ["run_dbcan", "CAZyme_annotation", "--mode", "protein",
               "--input_raw_data", str(native / "input.faa"), "--db_dir", str(database),
               "--output_dir", str(native), "--threads", str(threads),
               "--methods", "diamond,hmm,dbCANsub"]
    _run(command, native / "run.log")
    summaries, evidence = parse_outputs(native, gene_ids, id_map)
    metadata = {
        "tool_version": version, "command": command, "threads": threads,
        "diamond_version": subprocess.run(["diamond", "version"], check=True,
                                          capture_output=True, text=True).stdout.strip(),
        "policy": "all upstream passing hits; recommendations retained separately; no KO voting",
        "thresholds": {"diamond_evalue": 1e-102, "hmm_evalue": 1e-15,
                       "subfamily_evalue": 1e-15, "hmm_coverage": 0.35, "subfamily_coverage": 0.35},
        "protein_fasta": {"path": str(fasta), "sha256": sha256(fasta)},
        "database": {name: {"path": str(database / name), "sha256": sha256(database / name)} for name in ASSETS},
        "native_directory": str(output.parent / "dbcan"),
        "search_id_map": str(output.parent / "dbcan" / "gene_ids.tsv"),
    }
    # Keep failed runs for diagnosis; publish validated native data and summaries.
    _check_output_paths(database, fasta, output, published)
    files = list(native.iterdir())
    if any(not path.is_file() or (published / path.name).is_dir() for path in files):
        raise ValueError("Unexpected directory in dbCAN native output files")
    published.mkdir(parents=True, exist_ok=True)
    # Replace generated files individually; preserve any unrelated files added
    # to a previously owned native directory. Never remove the directory tree.
    for path in files:
        path.replace(published / path.name)
    (published / OWNER_FILE).write_text(OWNER_TEXT)
    native.rmdir()
    write_table(output, ["gene_id", *SUMMARY_COLUMNS], summaries)
    write_table(output.parent / "dbcan_evidence.tsv", EVIDENCE_COLUMNS, evidence)
    (output.parent / "dbcan_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
