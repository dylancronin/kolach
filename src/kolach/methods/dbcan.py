"""Run pinned dbCAN searches and retain domain-level, non-KO evidence.

Reference: DRAM2 4536cd93, modules/nf-core/rundbcan/cazymeannotation.
Only upstream-filtered hits are imported; recommendations are kept separately.
"""
from collections import namedtuple
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from kolach.download import download_atomic

VERSION = "5.2.9"
RELEASE = "db_v5-2-9_5-5-2026"
BASE_URL = f"https://dbcan.s3.us-west-2.amazonaws.com/{RELEASE}"
ASSETS = {
    "CAZy.dmnd": "CAZy.dmnd", "dbCAN.hmm": "dbCAN.hmm",
    "dbCAN-sub.hmm": "dbCAN_sub.hmm",
    "fam-substrate-mapping.tsv": "fam-substrate-mapping.tsv",
}
# Summary columns that collect the distinct values seen across a gene's hits.
SET_COLUMNS = [
    "dbcan_id", "dbcan_sub_id", "dbcan_sub_composition", "dbcan_sub_ec",
    "dbcan_sub_substrate", "dbcan_family_ec", "dbcan_cazy_id", "dbcan_methods",
]
SUMMARY_COLUMNS = [*SET_COLUMNS, "dbcan_tools_count", "dbcan_recommended"]
EVIDENCE_COLUMNS = [
    "gene_id", "method", "feature_id", "query_start", "query_end",
    "e_value", "coverage", "bit_score", "ec", "substrate", "source_row",
]
# Names used when domain evidence shares a table with KO evidence.
COMBINED_EVIDENCE_NAMES = {
    column: f"dbcan_{column}" for column in (
        "feature_id", "query_start", "query_end", "coverage", "ec", "substrate", "source_row",
    )
}
# Captures the family of a (sub)family name, e.g. GH13 from GH13_1.
FAMILY = re.compile(r"((?:GH|GT|CBM|AA|CE|PL)\d+)(?:_\d+)?")
EC = re.compile(r"\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:n?\d+|-)")
OWNER_FILE = ".kolach_dbcan_native"
OWNER_TEXT = "Kolach dbCAN native outputs v1\n"
RAW_SIDECARS = ("dbcan_evidence.tsv", "dbcan_metadata.json")


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


def check_dependencies():
    """Fail before launching the workflow when dbCAN executables are missing."""
    for executable in ("run_dbcan", "hmmstat", "diamond"):
        if not shutil.which(executable):
            raise FileNotFoundError(
                f"The '{executable}' executable was not found in PATH. "
                "Please activate the kolach conda environment with dbCAN installed."
            )


def tool_version():
    result = subprocess.run(["run_dbcan", "version"], check=True,
                            capture_output=True, text=True)
    if result.stdout.strip() != f"dbCAN version: {VERSION}":
        raise RuntimeError(f"Expected dbCAN {VERSION}, got {result.stdout.strip()!r}")
    return VERSION


def _read_table(path, required, aliases=None):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and aliases:
            reader.fieldnames = [aliases.get(column, column) for column in reader.fieldnames]
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
    """Download the pinned upstream bundle; publish only validated assets."""
    directory = Path(directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    # Invalidate any old success marker before attempting an update.
    (directory / ".download_complete").unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix="download-", dir=directory) as tmp:
        stage = Path(tmp)
        manifest = {"compatible_tool_version": VERSION, "release": RELEASE, "files": {}}
        for name, remote in ASSETS.items():
            url = f"{BASE_URL}/{remote}"
            digest = download_atomic(url, stage / name)
            manifest["files"][name] = {"url": url, "sha256": digest}
        validate_database(stage)
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        for name in [*ASSETS, "manifest.json"]:
            (stage / name).replace(directory / name)
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


# Native column names for one search; `score` is coverage (HMMs) or bit score (DIAMOND).
Source = namedtuple("Source", "file gene feature start end e_value score")
SOURCES = {
    "hmm": Source("dbCAN_hmm_results.tsv", "Target Name", "HMM Name",
                  "Target From", "Target To", "i-Evalue", "Coverage"),
    "dbCANsub": Source("dbCANsub_hmm_results.tsv", "Target Name", "Subfam Name",
                       "Target From", "Target To", "i-Evalue", "Coverage"),
    "diamond": Source("diamond.out", "Gene ID", "CAZy ID",
                      "Gene Start", "Gene End", "E Value", "Bit Score"),
}
# Extra columns required in the subfamily table.
SUBFAMILY_COLUMNS = ("Subfam Composition", "Subfam EC", "Substrate")


def _validate_record(record):
    """Feature and numeric checks shared by fresh searches and re-read evidence."""
    if not record["feature_id"] or record["feature_id"] == "-":
        raise ValueError(f"Missing dbCAN feature for {record['gene_id']!r}")
    _number(record["query_start"], "query start", minimum=1, integer=True)
    _number(record["query_end"], "query end", minimum=float(record["query_start"]), integer=True)
    _number(record["e_value"], "E-value")
    if record["method"] == "diamond":
        _number(record["bit_score"], "bit score")
    else:
        _number(record["coverage"], "coverage", maximum=1)


def _add_features(summary, method, gene, feature, row):
    """Add the summary values contributed by one native hit to a gene's sets."""
    summary["dbcan_methods"].add(method)
    if method == "hmm":
        summary["dbcan_id"] |= _parts(feature)
    elif method == "diamond":
        # Labels after the CAZy accession, minus EC numbers.
        labels = _parts(feature.partition("|")[2], "|")
        summary["dbcan_cazy_id"] |= {label for label in labels if not EC.fullmatch(label)}
    elif method == "dbCANsub":
        sub_ec = _parts(row["Subfam EC"])
        ecs = {part.split(":")[0] for part in sub_ec}
        if any(not EC.fullmatch(ec) for ec in ecs):
            raise ValueError(f"Malformed subfamily EC for {gene}: {row['Subfam EC']}")
        summary["dbcan_sub_id"] |= _parts(feature)
        summary["dbcan_sub_ec"] |= sub_ec
        summary["dbcan_sub_composition"] |= _parts(row["Subfam Composition"])
        summary["dbcan_sub_substrate"] |= _parts(row["Substrate"])
        # The family is encoded in the domain's subfamily name; its composition
        # field lists source-sequence counts, not families. Ambiguous
        # multi-family domains do not establish correlated family/EC pairs.
        families = set(FAMILY.findall(feature))
        if len(families) == 1:
            (family,) = families
            summary["dbcan_family_ec"] |= {f"{family};EC:{ec}" for ec in ecs}
    else:
        raise ValueError(f"Unknown dbCAN method: {method!r}")


def _restore(to_gene, source_gene):
    # With a search-ID map, every upstream ID must be one Kolach issued.
    if source_gene not in to_gene:
        raise ValueError(f"dbCAN returned unknown gene ID: {source_gene!r}")
    return to_gene[source_gene]


def parse_outputs(native, gene_ids, id_map=None):
    """Return stable summaries and all native records, without dropping domains."""
    native = Path(native)
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError("Duplicate gene IDs in dbCAN protein inventory")
    genes = {gene: {column: set() for column in SET_COLUMNS} for gene in gene_ids}
    to_gene = {search: gene for search, gene in (id_map or {g: g for g in gene_ids}).items() if gene in genes}
    evidence = []
    for method, source in SOURCES.items():
        required = [source.gene, source.feature, source.start, source.end, source.e_value,
                    source.score, *(SUBFAMILY_COLUMNS if method == "dbCANsub" else ())]
        for row in _read_table(native / source.file, required):
            gene = _restore(to_gene, row[source.gene])
            feature = row[source.feature].removesuffix(".hmm")
            record = {
                "gene_id": gene, "method": method, "feature_id": feature,
                "query_start": row[source.start], "query_end": row[source.end],
                "e_value": row[source.e_value], "coverage": row.get("Coverage", "-"),
                "bit_score": row.get("Bit Score", "-"), "ec": row.get("Subfam EC", "-"),
                "substrate": row.get("Substrate", "-"),
                "source_row": json.dumps(row, sort_keys=True, separators=(",", ":")),
            }
            _validate_record(record)
            _add_features(genes[gene], method, gene, feature, row)
            evidence.append(record)
    recommended = {}
    for row in _read_table(native / "overview.tsv", ["Gene ID", "#ofTools", "Recommend Results"]):
        gene = _restore(to_gene, row["Gene ID"])
        if gene in recommended:
            raise ValueError(f"Duplicate overview gene ID: {gene!r}")
        count = _number(row["#ofTools"], "tool count", maximum=3, integer=True)
        if count != len(genes[gene]["dbcan_methods"]):
            raise ValueError(f"Inconsistent overview method count for {gene}")
        recommended[gene] = row["Recommend Results"] or "-"
    if any(found["dbcan_methods"] and gene not in recommended for gene, found in genes.items()):
        raise ValueError("dbCAN overview omits genes with search hits")
    summaries = [{
        "gene_id": gene,
        **{column: ("|" if column == "dbcan_family_ec" else ";").join(sorted(found[column])) or "-"
           for column in SET_COLUMNS},
        "dbcan_tools_count": str(len(found["dbcan_methods"])),
        "dbcan_recommended": recommended.get(gene, "-"),
    } for gene, found in genes.items()]
    evidence.sort(key=lambda row: (row["gene_id"], row["method"], float(row["query_start"]),
                                   row["feature_id"], row["source_row"]))
    return summaries, evidence


def prepare_input(fasta, directory, gene_ids):
    """Protect IDs from upstream pandas type inference without changing sequences."""
    from kolach.integrate import _open_text

    id_map = {f"kolach_dbcan_{i:08d}": gene for i, gene in enumerate(gene_ids)}
    forward = {gene: safe for safe, gene in id_map.items()}
    directory = Path(directory)
    with _open_text(Path(fasta)) as source, open(directory / "input.faa", "w") as target:
        seen_header = False
        for number, line in enumerate(source, 1):
            if line.startswith(">"):
                parts = line[1:].split()
                if not parts:
                    raise ValueError(f"Empty FASTA header at line {number} of {fasta}")
                seen_header = True
                target.write(f">{forward[parts[0]]}\n")
            elif not seen_header and line.strip():
                raise ValueError(f"Sequence data before the first FASTA header in {fasta}")
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
        count = _number(row["dbcan_tools_count"], "tool count", maximum=3, integer=True)
        if count != len(_parts(row["dbcan_methods"])):
            raise ValueError(f"dbCAN tool count disagrees with methods for {row['gene_id']!r}")
    return [{key: row[key] for key in ["gene_id", *SUMMARY_COLUMNS]} for row in rows]


def read_evidence(path, gene_ids):
    """Read domain evidence with the same inventory and numeric checks as searches."""
    rows = _read_table(path, EVIDENCE_COLUMNS,
                       aliases={combined: native for native, combined in COMBINED_EVIDENCE_NAMES.items()})
    # A standard Kolach evidence table can be reused without retaining raw files.
    if rows and "dbcan_method" in rows[0]:
        rows = [{**row, "method": row["dbcan_method"]} for row in rows if row["method"] == "dbcan"]
    inventory = set(gene_ids)
    for row in rows:
        if row["gene_id"] not in inventory:
            raise ValueError(f"Unknown dbCAN evidence gene ID: {row['gene_id']!r}")
        if row["method"] not in {"hmm", "dbCANsub", "diamond"}:
            raise ValueError(f"Unknown dbCAN evidence method: {row['method']!r}")
        _validate_record(row)
        if not isinstance(json.loads(row["source_row"]), dict):
            raise ValueError("dbCAN source row must be a JSON object")
    return [{column: row[column] for column in EVIDENCE_COLUMNS} for row in rows]


def _check_output_paths(database, fasta, output, published):
    """Refuse collisions before running searches or replacing any output."""
    if published.is_symlink():
        raise ValueError(f"dbCAN native output directory must not be a symlink: {published}")
    destination = published.resolve()
    if destination.is_relative_to(database) or database.is_relative_to(destination):
        raise ValueError("dbCAN database and native output directories must not overlap")
    products = [output, output.parent / "dbcan_evidence.tsv", output.parent / ".dbcan_evidence.tsv",
                output.parent / "dbcan_metadata.json"]
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
        if any(published.iterdir()) and (not marker.is_file() or not marker.read_text().startswith(OWNER_TEXT)):
            raise ValueError(f"Refusing to overwrite a nonempty directory not owned by Kolach: {published}")


def _clear_raw_outputs(published, sidecars=()):
    """Remove files a previous Kolach run published; leave any user-added files."""
    marker = published / OWNER_FILE
    if marker.is_file():
        # The marker lists the generated files after its header line.
        for name in marker.read_text().splitlines()[1:]:
            path = published / name
            if name and path.parent == published and (path.is_file() or path.is_symlink()):
                path.unlink()
        if any(path != marker for path in published.iterdir()):
            marker.write_text(OWNER_TEXT)
        else:
            marker.unlink()
            published.rmdir()
    for path in sidecars:
        if path.is_file() or path.is_symlink():
            path.unlink()


def annotate(args, output):
    from kolach.integrate import read_fasta_ids

    fasta = Path(args.protein_fasta).expanduser().resolve()
    database = (Path(args.database_dir).expanduser().resolve() / "dbcan").resolve()
    output = Path(output).expanduser().resolve()
    published = output.parent / "dbcan"
    keep_raw = bool(getattr(args, "dbcan_keep_raw", False))
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
    _check_output_paths(database, fasta, output, published)
    write_table(output, ["gene_id", *SUMMARY_COLUMNS], summaries)
    write_table(output.parent / ".dbcan_evidence.tsv", EVIDENCE_COLUMNS, evidence)
    # Raw files from an earlier run would no longer describe this overview; in
    # particular, an old search-ID map would assign hits to the wrong genes.
    _clear_raw_outputs(published, [output.parent / name for name in RAW_SIDECARS])
    if not keep_raw:
        shutil.rmtree(native)
        return
    metadata = {
        "tool_version": version, "command": command, "threads": threads,
        "diamond_version": subprocess.run(["diamond", "version"], check=True,
                                          capture_output=True, text=True).stdout.strip(),
        "policy": "all upstream passing hits; recommendations retained separately; no KO voting",
        "thresholds": {"diamond_evalue": 1e-102, "hmm_evalue": 1e-15,
                       "subfamily_evalue": 1e-15, "hmm_coverage": 0.35, "subfamily_coverage": 0.35},
        "protein_fasta": {"path": str(fasta), "sha256": sha256(fasta)},
        "database": {name: {"path": str(database / name), "sha256": sha256(database / name)} for name in ASSETS},
        "native_directory": str(published),
        "search_id_map": str(published / "gene_ids.tsv"),
    }
    # Keep failed runs for diagnosis; publish validated native data and summaries.
    files = list(native.iterdir())
    if any(not path.is_file() or (published / path.name).is_dir() for path in files):
        raise ValueError("Unexpected directory in dbCAN native output files")
    published.mkdir(parents=True, exist_ok=True)
    # Replace generated files individually; preserve any unrelated files added
    # to a previously owned native directory. Never remove the directory tree.
    for path in files:
        path.replace(published / path.name)
    (published / OWNER_FILE).write_text(OWNER_TEXT + "".join(f"{path.name}\n" for path in files))
    native.rmdir()
    shutil.copyfile(output.parent / ".dbcan_evidence.tsv", output.parent / "dbcan_evidence.tsv")
    (output.parent / "dbcan_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
