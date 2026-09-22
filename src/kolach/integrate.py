"""Consensus and integration of annotation tables for kolach.

This module integrates results from KOfam, DeepKOALA, and eggNOG into a
unified, high-confidence annotation table using pandas, preserving all
method-specific metrics (bit scores, E-values, probabilities, and thresholds)
and evidence provenance. The gene summary table's KOfam metrics (kofam_bit_score,
kofam_evalue) follow each profile's score_type (full vs domain) as indicated by
kofam_score_type, while the long-form evidence table preserves all original full
and domain metrics.
"""
import gzip
from pathlib import Path
import re
from typing import Any, Optional, Union

import numpy as np
import pandas as pd

KO_REGEX = re.compile(r"\bK\d{5}\b")

EVIDENCE_COLUMNS = [
    "gene_id",
    "ko",
    "method",
    "bit_score",
    "e_value",
    "domain_bit_score",
    "domain_e_value",
    "score_type",
    "bit_score_threshold",
    "e_value_threshold",
    "deepkoala_probability",
    "deepkoala_threshold",
    "eggnog_shared_seed_hit",
    "original_status",
    "cross_method_status",
    "consensus_level",
    "supporting_methods",
]


def _validate_gene_ids(values: pd.Series, path: Path, method: str) -> None:
    """Reject missing identifiers before pandas string conversion can create 'nan'."""
    invalid = values.isna() | values.astype("string").str.strip().eq("")
    if invalid.any():
        rows = ", ".join(str(i + 2) for i in values.index[invalid][:5])
        raise ValueError(f"{method} file '{path}' contains missing or blank gene IDs at row(s): {rows}")


def parse_kos(val: Any) -> set[str]:
    """Extract valid KEGG Orthology identifiers (matching K\d{5}) from a string or value.

    Splits strings on common delimiters (commas, semicolons, whitespace, pipes, pluses)
    and uses regular expression matching to isolate standard 6-character KO numbers.

    Args:
        val: Input value (string, float, NaN, or None) containing raw KO annotations.

    Returns:
        set[str]: Set of valid KO identifiers (e.g. {'K00001', 'K00002'}), or an empty set
        if the input is missing, empty, or contains no valid KO tokens.
    """
    # Guard against missing values, nulls, and standard missing-value strings
    if pd.isna(val) or val is None or str(val).strip() in ("", "-", "None", "nan"):
        return set()

    # Split on diverse delimiter formats used across tools (e.g., commas in eggNOG, semicolons/spaces)
    tokens = re.split(r"[,+;\s|]+", str(val).strip())
    kos = set()
    for t in tokens:
        m = KO_REGEX.search(t)
        if m:
            kos.add(m.group(0))
    return kos


def format_kos(kos: set[str]) -> str:
    """Format a set of KO identifiers into a deterministic, sorted, comma-separated string.

    Args:
        kos (set[str]): Set of KO identifier strings to format.

    Returns:
        str: Alphabetically sorted, comma-separated KOs (e.g. 'K00001,K00002'),
        or '-' if the set is empty.
    """
    # Alphabetical sorting guarantees reproducible, byte-identical output across runs
    return ",".join(sorted(kos)) if kos else "-"


def format_metric_value(value: Any) -> str:
    """Format a numeric value compactly for a KO-qualified summary field."""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "-"
    return f"{numeric:.6g}" if np.isfinite(numeric) else "-"


def read_fasta_ids(fasta_path: Union[str, Path]) -> list[str]:
    """Read all protein/gene sequence IDs from a FASTA file in original sequence order.

    Determines whether the input is gzip-compressed based on file extension and extracts
    the primary identifier (the first non-whitespace token after '>') from each header line.

    Args:
        fasta_path (Union[str, Path]): Path to the protein FASTA file (.fasta, .faa, or .gz).

    Returns:
        list[str]: Sequence identifiers in the exact order they appear in the FASTA file.
    """
    fasta_path = Path(fasta_path).expanduser().resolve()
    gene_ids = []

    # Check for gzip compression by inspecting the suffix
    is_gz = fasta_path.suffix.lower() in (".gz", ".gzip") or str(fasta_path).endswith((".fasta.gz", ".faa.gz", ".fa.gz"))
    open_fn = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8", errors="replace"))
        if is_gz
        else (lambda p: open(p, "r", encoding="utf-8", errors="replace"))
    )

    with open_fn(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # Header line: strip '>' and take the first whitespace-delimited word as the ID
                parts = line[1:].split()
                if parts:
                    gene_ids.append(parts[0])
    return gene_ids


def load_ko_definitions(ko_list_path: Optional[Union[str, Path]]) -> dict[str, str]:
    """Load authoritative KO functional definitions from KOfam's master ko_list reference file.

    Parses the tab-delimited ko_list file where column 1 is the KO number (Kxxxxx)
    and column 2 (or the final column) contains the functional definition text.

    Args:
        ko_list_path (Optional[Union[str, Path]]): Path to the ko_list file, or None if omitted.

    Returns:
        dict[str, str]: Mapping from KO ID (str) to definition text (str). Returns empty dict
        if path is None or the file does not exist.
    """
    definitions = {}
    if not ko_list_path:
        return definitions
    path = Path(ko_list_path).expanduser().resolve()
    if not path.is_file():
        return definitions

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 2:
                    ko = parts[0].strip()
                    # Validate that the first column matches standard Kxxxxx format
                    if KO_REGEX.fullmatch(ko):
                        definition = parts[-1].strip() if len(parts) > 1 else ""
                        definitions[ko] = definition
    except Exception:
        # Fall back gracefully to empty definitions on read or parsing errors
        pass
    return definitions


def select_kofam_scores(
    score_type: Optional[str],
    bit_score: Any,
    e_value: Any,
    domain_bit_score: Any,
    domain_e_value: Any,
) -> tuple[float, float, str]:
    """Select the applicable KOfam bit score and E-value based on profile score_type.

    KOfam models use either full-sequence ('full') or domain-level ('domain') thresholds.
    This function routes the corresponding score and E-value according to the profile definition.
    Crucially, it never substitutes a full-sequence score when a required domain score is missing.

    Args:
        score_type (Optional[str]): Profile score model type, typically 'full' or 'domain'.
        bit_score (Any): Full-sequence alignment bit score.
        e_value (Any): Full-sequence alignment E-value.
        domain_bit_score (Any): Domain-level alignment bit score.
        domain_e_value (Any): Domain-level alignment E-value.

    Returns:
        tuple[float, float, str]: A 3-tuple containing:
            - selected_bit_score (float): Applicable bit score (or np.nan if missing/nonfinite).
            - selected_e_value (float): Applicable E-value (or np.nan if missing/nonfinite).
            - normalized_score_type (str): 'full', 'domain', or '-' for unrecognized score types.
    """
    st_raw = str(score_type).strip().lower() if pd.notna(score_type) else ""
    if st_raw == "full":
        # Full-sequence model: select full bit score and full E-value
        bs = float(bit_score) if (pd.notna(bit_score) and np.isfinite(bit_score)) else np.nan
        ev = float(e_value) if (pd.notna(e_value) and np.isfinite(e_value)) else np.nan
        return bs, ev, "full"
    elif st_raw == "domain":
        # Domain model: select domain bit score and domain E-value (never substitute full scores)
        dbs = float(domain_bit_score) if (pd.notna(domain_bit_score) and np.isfinite(domain_bit_score)) else np.nan
        dev = float(domain_e_value) if (pd.notna(domain_e_value) and np.isfinite(domain_e_value)) else np.nan
        return dbs, dev, "domain"
    else:
        # Unrecognized or missing score type
        return np.nan, np.nan, "-"


def make_evidence_record(
    gene_id: str,
    ko: str,
    method: str,
    original_status: str,
    bit_score: float = np.nan,
    e_value: float = np.nan,
    domain_bit_score: float = np.nan,
    domain_e_value: float = np.nan,
    score_type: str = "-",
    threshold: float = np.nan,
    bit_score_threshold: float = np.nan,
    e_value_threshold: float = np.nan,
    deepkoala_probability: float = np.nan,
    deepkoala_threshold: float = np.nan,
    eggnog_shared_seed_hit: bool = False,
    shared_seed_hit: bool = False,
    definition: str = "-",
    **tool_kwargs: Any,
) -> dict[str, Any]:
    """Construct a normalized evidence record with standard field defaults.

    Provides uniform dictionary structure across all three annotation extractors
    while allowing tool-specific metadata (such as KOfam's 'assignment' or eggNOG's
    'eggnog_description') via keyword arguments.

    Args:
        gene_id (str): Gene/protein sequence identifier.
        ko (str): Target KEGG Orthology identifier (Kxxxxx).
        method (str): Prediction tool name ('kofam', 'deepkoala', or 'eggnog').
        original_status (str): Initial threshold evaluation status ('threshold_passing',
            'heuristic_rescued', 'below_threshold', or 'unknown').
        bit_score (float, optional): Alignment bit score. Defaults to np.nan.
        e_value (float, optional): Alignment E-value. Defaults to np.nan.
        domain_bit_score (float, optional): Domain alignment bit score. Defaults to np.nan.
        domain_e_value (float, optional): Domain alignment E-value. Defaults to np.nan.
        score_type (str, optional): Metric model type ('full', 'domain', 'seed_hit', etc.). Defaults to '-'.
        threshold (float, optional): Profile or method score threshold. Defaults to np.nan.
        bit_score_threshold (float, optional): Applicable bit score threshold. Defaults to np.nan.
        e_value_threshold (float, optional): Applicable E-value cutoff. Defaults to np.nan.
        deepkoala_probability (float, optional): DeepKOALA model prediction score. Defaults to np.nan.
        deepkoala_threshold (float, optional): DeepKOALA decision cutoff. Defaults to np.nan.
        eggnog_shared_seed_hit (bool, optional): True if hit originated from a multi-KO seed hit. Defaults to False.
        shared_seed_hit (bool, optional): Alias for eggnog_shared_seed_hit. Defaults to False.
        definition (str, optional): Functional description text. Defaults to '-'.
        **tool_kwargs (Any): Tool-specific attributes (e.g. assignment='threshold', eggnog_description='...').

    Returns:
        dict[str, Any]: Standardized dictionary representing one record in the evidence table.
    """
    rec = {
        "gene_id": gene_id,
        "ko": ko,
        "method": method,
        "original_status": original_status,
        "bit_score": bit_score,
        "e_value": e_value,
        "domain_bit_score": domain_bit_score,
        "domain_e_value": domain_e_value,
        "score_type": score_type if score_type else "-",
        "threshold": threshold,
        "bit_score_threshold": bit_score_threshold,
        "e_value_threshold": e_value_threshold,
        "deepkoala_probability": deepkoala_probability,
        "deepkoala_threshold": deepkoala_threshold,
        "eggnog_shared_seed_hit": eggnog_shared_seed_hit,
        "shared_seed_hit": shared_seed_hit,
        "definition": definition,
    }
    # Append tool-specific attributes without leaking default keys into other tools' records
    rec.update(tool_kwargs)
    return rec


def extract_kofam_records(
    tsv_path: Union[str, Path],
    heuristic_bitscore_fraction: float = 0.75,
    heuristic_e_value: float = 1e-5,
) -> list[dict]:
    """Extract per-hit KO records from KOfam annotations TSV.

    Duplicate hits for a (gene_id, ko) pair prefer threshold-passing evidence,
    then heuristic-rescued, below-threshold, and unknown evidence. Applicable bit
    score, E-value, and input order break ties within a status category.

    original_status is strictly determined:
    - 'threshold' or '*' -> threshold_passing
    - 'rescued' -> heuristic_rescued
    - 'below_threshold', 'below', or 'none' -> below_threshold
    - missing / unrecognized: evaluated using score_type and threshold.
      If valid numeric applicable score and threshold exist:
        threshold_passing if score >= threshold else below_threshold
      Else (including missing required domain score for domain profiles):
        unknown
    Preserves original full and domain metrics in the underlying evidence records.

    For heuristic_rescued hits the reported bit_score_threshold is the actual
    relaxed cutoff used upstream (heuristic_bitscore_fraction * threshold) and
    e_value_threshold is the actual heuristic E-value cutoff, so the evidence
    table's thresholds match the configured KOfam rescue parameters.

    Args:
        tsv_path (Union[str, Path]): Path to the KOfam annotations TSV file.
        heuristic_bitscore_fraction (float, optional): Fraction of the profile
            threshold used for rescued hits. Defaults to 0.75.
        heuristic_e_value (float, optional): Maximum E-value for rescue
            consideration. Defaults to 1e-5.

    Returns:
        list[dict]: List of normalized evidence records conforming to EVIDENCE_COLUMNS.
    """
    path = Path(tsv_path).expanduser().resolve()
    # Return empty if file does not exist or has zero bytes
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"gene_id", "ko"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"KOfam file '{path}' missing required column(s): {', '.join(sorted(missing))}")
    if df.empty:
        return []
    _validate_gene_ids(df["gene_id"], path, "KOfam")

    # Clean identifiers and coerce score columns to numeric floats (NaN on conversion failure)
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["ko"] = df["ko"].astype(str).str.strip()
    df["bit_score"] = pd.to_numeric(df.get("bit_score"), errors="coerce")
    df["e_value"] = pd.to_numeric(df.get("e_value"), errors="coerce")
    df["domain_bit_score"] = pd.to_numeric(df.get("domain_bit_score"), errors="coerce")
    df["domain_e_value"] = pd.to_numeric(df.get("domain_e_value"), errors="coerce")
    df["threshold"] = pd.to_numeric(df.get("threshold"), errors="coerce")
    if "score_type" not in df.columns:
        df["score_type"] = "-"
    if "definition" not in df.columns:
        df["definition"] = "-"
    if "assignment" not in df.columns:
        df["assignment"] = "-"

    # Filter to valid KEGG format identifiers (matching Kxxxxx)
    df = df[df["ko"].str.fullmatch(KO_REGEX, na=False)].copy()
    if df.empty:
        return []

    # Compute selected score and E-value for deterministic deduplication and inferred status
    sel_scores = [
        select_kofam_scores(st, bs, ev, dbs, dev)
        for st, bs, ev, dbs, dev in zip(
            df["score_type"], df["bit_score"], df["e_value"], df["domain_bit_score"], df["domain_e_value"]
        )
    ]
    df["_sel_bit_score"] = [s[0] for s in sel_scores]
    df["_sel_e_value"] = [s[1] for s in sel_scores]
    df["_orig_idx"] = np.arange(len(df))

    def status_for_row(row: pd.Series) -> str:
        assignment = str(row["assignment"]).strip().lower() if pd.notna(row["assignment"]) else ""
        if assignment in ("threshold", "*"):
            return "threshold_passing"
        if assignment == "rescued":
            return "heuristic_rescued"
        if assignment in ("below_threshold", "below", "none"):
            return "below_threshold"
        if (
            str(row["score_type"]).strip().lower() in ("full", "domain")
            and pd.notna(row["_sel_bit_score"])
            and np.isfinite(row["_sel_bit_score"])
            and pd.notna(row["threshold"])
            and np.isfinite(row["threshold"])
        ):
            return "threshold_passing" if row["_sel_bit_score"] >= row["threshold"] else "below_threshold"
        return "unknown"

    df["_original_status"] = df.apply(status_for_row, axis=1)
    df["_status_rank"] = df["_original_status"].map(
        {"threshold_passing": 3, "heuristic_rescued": 2, "below_threshold": 1, "unknown": 0}
    )

    # Preserve stronger status before ranking metrics within each category.
    df = df.sort_values(
        by=["gene_id", "ko", "_status_rank", "_sel_bit_score", "_sel_e_value", "_orig_idx"],
        ascending=[True, True, False, False, True, True],
        na_position="last",
    )
    df = df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in df.iterrows():
        assign_raw = str(row["assignment"]).strip() if pd.notna(row["assignment"]) else ""
        score_type = str(row["score_type"]).strip() if pd.notna(row["score_type"]) else "-"
        thresh = row["threshold"]
        bs = row["bit_score"]
        ev = row["e_value"]
        dbs = row["domain_bit_score"]
        dev = row["domain_e_value"]
        orig_status = row["_original_status"]

        definition = str(row["definition"]).strip() if pd.notna(row["definition"]) else "-"
        if definition in ("", "-"):
            definition = "-"

        # Determine bit_score_threshold and e_value_threshold for evidence reporting
        # For heuristic rescued hits, report the actual relaxed cutoff used upstream
        if orig_status == "heuristic_rescued":
            bs_thresh = heuristic_bitscore_fraction * thresh if (pd.notna(thresh) and np.isfinite(thresh)) else np.nan
            ev_thresh = float(heuristic_e_value)
        else:
            bs_thresh = thresh
            ev_thresh = np.nan

        records.append(
            make_evidence_record(
                gene_id=row["gene_id"],
                ko=row["ko"],
                method="kofam",
                original_status=orig_status,
                bit_score=bs,
                e_value=ev,
                domain_bit_score=dbs,
                domain_e_value=dev,
                score_type=score_type,
                threshold=thresh,
                bit_score_threshold=bs_thresh,
                e_value_threshold=ev_thresh,
                definition=definition,
                assignment=assign_raw if assign_raw else "-",
            )
        )
    return records


def extract_deepkoala_records(tsv_path: Union[str, Path]) -> list[dict]:
    """Extract per-hit KO records from DeepKOALA annotations TSV.

    Duplicate hits for a (gene_id, ko) pair prefer threshold-passing evidence,
    then below-threshold evidence, then unknown evidence; probability and input
    order break ties within a status category.

    Status determination:
    - If annotate column is present and contains '*': threshold_passing
    - If annotate column is empty/non-passing or absent:
      If valid numeric score and threshold exist:
        threshold_passing if score >= threshold else below_threshold
      Else:
        unknown
    - Documented simple format (name, predict_label exactly):
      emitted KOs have passed upstream filtering -> threshold_passing with
      score_type = 'upstream_filtered'.
    - Tables without documented gene/KO headers are rejected.

    Args:
        tsv_path (Union[str, Path]): Path to the DeepKOALA annotations TSV file.

    Returns:
        list[dict]: List of normalized evidence records conforming to EVIDENCE_COLUMNS.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = pd.read_csv(path, sep="\t", dtype=str)
    # Clean header column names for case-insensitive and hash-stripped matching
    cleaned_raw_cols = [c.lstrip("#").strip() for c in df.columns]
    cols_lower = [c.lower() for c in cleaned_raw_cols]

    # Check if this matches documented DeepKOALA simple-output format (name + predict_label only)
    is_documented_simple = (
        len(df.columns) == 2
        and "name" in cols_lower
        and "predict_label" in cols_lower
        and "annotate" not in cols_lower
        and "probability" not in cols_lower
        and "score" not in cols_lower
    )

    # Build column rename mapping for flexible DeepKOALA headers
    col_map = {}
    if "name" in cols_lower:
        col_map[df.columns[cols_lower.index("name")]] = "gene_id"
    elif "gene_id" in cols_lower:
        col_map[df.columns[cols_lower.index("gene_id")]] = "gene_id"
    else:
        raise ValueError(f"DeepKOALA file '{path}' missing required 'name' or 'gene_id' column")

    if "predict_label" in cols_lower:
        col_map[df.columns[cols_lower.index("predict_label")]] = "predict_label"
    elif "ko" in cols_lower:
        col_map[df.columns[cols_lower.index("ko")]] = "predict_label"
    else:
        raise ValueError(f"DeepKOALA file '{path}' missing required 'predict_label' or 'ko' column")

    if "probability" in cols_lower:
        col_map[df.columns[cols_lower.index("probability")]] = "deepkoala_score"
    elif "score" in cols_lower:
        col_map[df.columns[cols_lower.index("score")]] = "deepkoala_score"

    if "threshold" in cols_lower:
        col_map[df.columns[cols_lower.index("threshold")]] = "deepkoala_threshold"

    if "annotate" in cols_lower:
        col_map[df.columns[cols_lower.index("annotate")]] = "deepkoala_annotate"

    df = df.rename(columns=col_map)
    if df.empty:
        return []
    _validate_gene_ids(df["gene_id"], path, "DeepKOALA")
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["deepkoala_score"] = pd.to_numeric(df.get("deepkoala_score"), errors="coerce")
    df["deepkoala_threshold"] = pd.to_numeric(df.get("deepkoala_threshold"), errors="coerce")

    exploded_rows = []
    for _, row in df.iterrows():
        raw_val = row.get("predict_label", "-")
        # Explode multi-KO predict labels so each (gene_id, ko) pair is an independent record
        kos = parse_kos(raw_val)
        score = row.get("deepkoala_score", np.nan)
        thresh = row.get("deepkoala_threshold", np.nan)
        ann = str(row.get("deepkoala_annotate", "")).strip() if "deepkoala_annotate" in df.columns and pd.notna(row.get("deepkoala_annotate")) else ""

        # Determine original status per row based on DeepKOALA output format
        if is_documented_simple:
            # Documented 2-col output contains KOs that already passed upstream model threshold
            orig_status = "threshold_passing"
            score_type = "upstream_filtered"
        elif "deepkoala_annotate" in df.columns:
            if "*" in ann:
                orig_status = "threshold_passing"
                score_type = "probability"
            elif pd.notna(score) and np.isfinite(score) and pd.notna(thresh) and np.isfinite(thresh):
                orig_status = "threshold_passing" if score >= thresh else "below_threshold"
                score_type = "probability"
            else:
                orig_status = "unknown"
                score_type = "-"
        elif pd.notna(score) and np.isfinite(score) and pd.notna(thresh) and np.isfinite(thresh):
            orig_status = "threshold_passing" if score >= thresh else "below_threshold"
            score_type = "probability"
        else:
            orig_status = "unknown"
            score_type = "-"

        for ko in kos:
            exploded_rows.append({
                "gene_id": row["gene_id"],
                "ko": ko,
                "deepkoala_score": score,
                "deepkoala_threshold": thresh,
                "score_type": score_type,
                "original_status": orig_status,
            })

    if not exploded_rows:
        return []

    exp_df = pd.DataFrame(exploded_rows)
    # Preserve passing evidence before ranking within each status category.
    exp_df["_status_rank"] = exp_df["original_status"].map(
        {"threshold_passing": 2, "below_threshold": 1, "unknown": 0}
    )
    exp_df["_orig_idx"] = np.arange(len(exp_df))
    exp_df = exp_df.sort_values(
        by=["gene_id", "ko", "_status_rank", "deepkoala_score", "_orig_idx"],
        ascending=[True, True, False, False, True],
        na_position="last",
    )
    exp_df = exp_df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in exp_df.iterrows():
        records.append(
            make_evidence_record(
                gene_id=row["gene_id"],
                ko=row["ko"],
                method="deepkoala",
                original_status=row["original_status"],
                score_type=row["score_type"],
                deepkoala_probability=row["deepkoala_score"],
                deepkoala_threshold=row["deepkoala_threshold"],
            )
        )
    return records


def read_eggnog_tsv(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Read eggNOG annotations table, properly handling native and normalized headers.

    Supports:
    - Native .emapper.annotations files containing '##' metadata and a '#query' header.
    - Kolach-normalized TSVs containing an unprefixed 'query' header.
    - Preserves column aliases ('query'/'gene_id', 'kegg_ko'/'ko', 'score', 'evalue', 'description').
    - Raises ValueError when required identifying or KO columns are missing.
    - Valid header-only files return an empty DataFrame with the parsed columns.

    Args:
        tsv_path (Union[str, Path]): Path to the eggNOG annotations file (plain text or gzip).

    Returns:
        pd.DataFrame: Parsed DataFrame of eggNOG records with clean, un-prefixed column names.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the file is empty (0 bytes) or lacks required identifying/KO headers.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"eggNOG file not found: '{path}'")
    if path.stat().st_size == 0:
        raise ValueError(f"eggNOG annotations file is empty (0 bytes): '{path}'")

    # Detect gzip compression from extension
    is_gz = path.suffix.lower() in (".gz", ".gzip") or str(path).endswith(".annotations.gz")
    open_fn = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8", errors="replace"))
        if is_gz
        else (lambda p: open(p, "r", encoding="utf-8", errors="replace"))
    )

    header_line = None
    skip_indices: list[int] = []
    with open_fn(path) as f:
        for idx, line in enumerate(f):
            stripped = line.strip()
            if not stripped:
                skip_indices.append(idx)
                continue
            # Skip emapper comment metadata lines
            if stripped.startswith("##"):
                skip_indices.append(idx)
                continue
            if header_line is None:
                # Skip comment lines before the header except the recognized #query header
                if stripped.startswith("#"):
                    first_cell = re.sub(r"\s+", "", stripped.split("\t")[0].lower())
                    if first_cell not in ("#query", "#gene_id"):
                        skip_indices.append(idx)
                        continue
                header_line = line.rstrip("\r\n")
                header_idx = idx
            else:
                if stripped.startswith("#"):
                    skip_indices.append(idx)

    if header_line is None:
        raise ValueError(f"eggNOG annotations file contains no valid header: '{path}'")

    raw_cols = header_line.split("\t")
    clean_cols = [c.lstrip("#").strip() for c in raw_cols]

    cols_lower = [c.lower() for c in clean_cols]
    has_id = any(c in cols_lower for c in ("query", "gene_id"))
    has_ko = any(c in cols_lower for c in ("kegg_ko", "ko"))
    if not has_id or not has_ko:
        raise ValueError(
            f"eggNOG annotations file '{path}' missing required identifying ('query'/'gene_id') "
            f"or KO ('kegg_ko'/'ko') column (found columns: {clean_cols})"
        )

    # Re-open and parse natively from the stream, skipping metadata/comment rows in
    # a single pass (avoids buffering the whole file into memory via StringIO).
    skip_rows = sorted(set(skip_indices) | {header_idx})
    with open_fn(path) as f:
        df = pd.read_csv(f, sep="\t", header=None, names=clean_cols, dtype=str, skiprows=skip_rows)

    # Valid header-only files return an empty DataFrame with the parsed columns
    if df.empty:
        return pd.DataFrame(columns=clean_cols)
    return df


def extract_eggnog_records(
    tsv_path: Union[str, Path],
    min_bitscore: float = 60.0,
    max_evalue: float = 1e-5,
) -> list[dict]:
    """Extract per-hit KO records from eggNOG annotations TSV.

    Shared seed-hit metrics for multi-KO assignments are explicitly tagged.
    Duplicate hits for a (gene_id, ko) pair prefer threshold-passing evidence,
    then below-threshold evidence, then unknown evidence; bit score, E-value,
    and input order break ties within a status category. All source seed-group
    memberships remain available internally for group-scoped disambiguation.

    Status determination:
    - Missing, malformed, or nonfinite bit score or E-value -> unknown
    - Valid finite metrics -> threshold_passing if bit_score >= min_bitscore and evalue <= max_evalue,
      else below_threshold.
    Preserves legitimate E-values of zero.
    eggNOG descriptions are generic seed-hit descriptions and are preserved as
    metadata, NOT authoritative KO-specific definitions.

    Args:
        tsv_path (Union[str, Path]): Path to eggNOG annotations TSV file.
        min_bitscore (float, optional): Minimum alignment bit score for confident hits. Defaults to 60.0.
        max_evalue (float, optional): Maximum alignment E-value cutoff for confident hits. Defaults to 1e-5.

    Returns:
        list[dict]: List of normalized evidence records conforming to EVIDENCE_COLUMNS.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    # Read eggNOG TSV handling native headers, gzip, and comment lines
    df = read_eggnog_tsv(path)
    if df.empty:
        return []

    # Build case-insensitive column map for flexible eggNOG outputs
    cols_lower = [c.lower() for c in df.columns]
    col_map = {}
    if "query" in cols_lower:
        col_map[df.columns[cols_lower.index("query")]] = "gene_id"
    elif "gene_id" in cols_lower:
        col_map[df.columns[cols_lower.index("gene_id")]] = "gene_id"

    if "kegg_ko" in cols_lower:
        col_map[df.columns[cols_lower.index("kegg_ko")]] = "eggnog_raw_ko"
    elif "ko" in cols_lower:
        col_map[df.columns[cols_lower.index("ko")]] = "eggnog_raw_ko"

    if "score" in cols_lower:
        col_map[df.columns[cols_lower.index("score")]] = "eggnog_bit_score"
    if "evalue" in cols_lower:
        col_map[df.columns[cols_lower.index("evalue")]] = "eggnog_evalue"
    if "description" in cols_lower:
        col_map[df.columns[cols_lower.index("description")]] = "eggnog_description"

    df = df.rename(columns=col_map)
    _validate_gene_ids(df["gene_id"], path, "eggNOG")
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["eggnog_bit_score"] = pd.to_numeric(df.get("eggnog_bit_score"), errors="coerce")
    df["eggnog_evalue"] = pd.to_numeric(df.get("eggnog_evalue"), errors="coerce")
    if "eggnog_description" not in df.columns:
        df["eggnog_description"] = "-"
    if "eggnog_raw_ko" not in df.columns:
        df["eggnog_raw_ko"] = "-"

    exploded_rows = []
    for row_idx, (_, row) in enumerate(df.iterrows()):
        raw_val = row["eggnog_raw_ko"]
        # Explode multi-KO annotations into individual per-KO records
        kos = parse_kos(raw_val)
        if not kos:
            continue
        # Flag if this seed hit called multiple KOs (triggers disambiguation)
        is_multi = len(kos) > 1
        bs = row["eggnog_bit_score"]
        ev = row["eggnog_evalue"]
        # Stable identifier for the originating eggNOG row; all KOs exploded from
        # the same seed hit share this group id so disambiguation can be scoped to
        # the KOs that genuinely co-occurred in a multi-KO hit.
        group_key = f"en_row{row_idx}"

        # Validate finite numeric metrics (legitimate 0.0 E-values evaluate to True)
        is_finite_bs = pd.notna(bs) and np.isfinite(bs)
        is_finite_ev = pd.notna(ev) and np.isfinite(ev)

        if not (is_finite_bs and is_finite_ev):
            orig_status = "unknown"
        elif bs >= min_bitscore and ev <= max_evalue:
            orig_status = "threshold_passing"
        else:
            orig_status = "below_threshold"

        for ko in kos:
            exploded_rows.append({
                "gene_id": row["gene_id"],
                "ko": ko,
                "bit_score": bs,
                "e_value": ev,
                "eggnog_shared_seed_hit": is_multi,
                "shared_seed_hit": is_multi,
                "eggnog_seed_group": group_key,
                "source_is_multi": is_multi,
                "original_status": orig_status,
                "eggnog_description": str(row["eggnog_description"]).strip() if pd.notna(row["eggnog_description"]) else "-",
            })

    if not exploded_rows:
        return []

    exp_df = pd.DataFrame(exploded_rows)
    source_map: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for source in exploded_rows:
        source_map.setdefault((source["gene_id"], source["ko"]), []).append(
            {
                "group": source["eggnog_seed_group"],
                "is_multi": source["source_is_multi"],
                "status": source["original_status"],
            }
        )

    # Preserve passing evidence before ranking metrics within each status category.
    exp_df["_status_rank"] = exp_df["original_status"].map(
        {"threshold_passing": 2, "below_threshold": 1, "unknown": 0}
    )
    exp_df["_orig_idx"] = np.arange(len(exp_df))
    exp_df = exp_df.sort_values(
        by=["gene_id", "ko", "_status_rank", "bit_score", "e_value", "_orig_idx"],
        ascending=[True, True, False, False, True, True],
        na_position="last",
    )
    exp_df = exp_df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in exp_df.iterrows():
        sources = source_map[(row["gene_id"], row["ko"])]
        is_shared = any(source["is_multi"] for source in sources)
        records.append(
            make_evidence_record(
                gene_id=row["gene_id"],
                ko=row["ko"],
                method="eggnog",
                original_status=row["original_status"],
                bit_score=row["bit_score"],
                e_value=row["e_value"],
                score_type="seed_hit",
                threshold=min_bitscore,
                bit_score_threshold=min_bitscore,
                e_value_threshold=max_evalue,
                eggnog_shared_seed_hit=is_shared,
                shared_seed_hit=is_shared,
                eggnog_seed_group=row["eggnog_seed_group"],
                eggnog_sources=sources,
                definition="-",  # Generic seed description is NOT an authoritative KO definition
                eggnog_description=row["eggnog_description"],
            )
        )
    return records


def get_single_ko_definition(
    ko: str,
    ko_definitions: dict[str, str],
    evidence_definitions: Optional[dict[str, str]] = None,
) -> str:
    """Resolve functional description for a single KO identifier following authority hierarchy.

    Strict priority:
    1. Master ko_list definition.
    2. Matching KO-specific definition from evidence records (e.g. KOfam profile definition).
    3. Return '-' if no matching definition is available.

    Args:
        ko (str): Target KO identifier (Kxxxxx).
        ko_definitions (dict[str, str]): Master definitions dictionary from ko_list.
        evidence_definitions (Optional[dict[str, str]], optional): Fallback definitions from KOfam profiles.

    Returns:
        str: Functional definition text, or '-' if not found.
    """
    # 1. Master ko_list definition is highest authority
    if ko in ko_definitions and ko_definitions[ko] and ko_definitions[ko] != "-":
        return str(ko_definitions[ko]).strip()
    # 2. Fall back to profile definition embedded in KOfam HMM output
    if evidence_definitions and ko in evidence_definitions:
        edef = evidence_definitions[ko]
        if edef and edef != "-":
            return str(edef).strip()
    # 3. Default placeholder
    return "-"


def resolve_definition_for_kos(
    kos: set[str],
    ko_definitions: dict[str, str],
    evidence_definitions: Optional[dict[str, str]] = None,
) -> str:
    """Resolve functional definition text for a set of KOs.

    Preserves unambiguous KO-definition associations:
    - Single KO: returns definition string directly (or '-')
    - Multiple KOs: returns 'Kxxxxx: def1; Kyyyyy: def2' for KOs with definitions,
      or '-' if none found.

    Args:
        kos (set[str]): Set of KO identifiers to annotate.
        ko_definitions (dict[str, str]): Master definitions dictionary from ko_list.
        evidence_definitions (Optional[dict[str, str]], optional): Fallback definitions from KOfam profiles.

    Returns:
        str: Formatted definition string, or '-' if the set is empty or lacks definitions.
    """
    if not kos:
        return "-"
    sorted_kos = sorted(kos)
    # Single KO case: return clean definition string without KO prefix
    if len(sorted_kos) == 1:
        return get_single_ko_definition(sorted_kos[0], ko_definitions, evidence_definitions)

    # Multi-KO case: format each known KO as 'Kxxxxx: def' and join with semicolons
    parts = []
    for k in sorted_kos:
        d = get_single_ko_definition(k, ko_definitions, evidence_definitions)
        if d != "-":
            parts.append(f"{k}: {d}")
    return "; ".join(parts) if parts else "-"


def disambiguate_eggnog(
    en_agreed: set[str],
    tool_recs: dict[str, list[dict]],
    active_tools: list[str],
    confident_kos: dict[str, set[str]],
    candidate_kos: dict[str, set[str]],
    rescued_kos: dict[str, set[str]],
    eggnog_filter_multi: str = "disambiguate",
) -> tuple[set[str], set[str], set[Any], set[Any]]:
    """Disambiguate eggNOG multi-KO assignments against other methods' predictions.

    When an eggNOG assignment represents a multi-KO seed hit, checks for corroboration
    first against confident calls from other tools, then against candidate/rescued calls.
    Narrowing is scoped to the KOs that genuinely shared the multi-KO seed hit:
    confident singleton calls from other eggNOG rows for the same gene are preserved
    unconditionally.
    Dropped KOs are returned in en_dropped and preserved in alternative_kos.

    Args:
        en_agreed (set[str]): Confident KOs called by eggNOG for this gene.
        tool_recs (dict[str, list[dict]]): Categorized evidence records per tool for this gene.
        active_tools (list[str]): List of tools participating in integration.
        confident_kos (dict[str, set[str]]): Confident KOs called by each tool.
        candidate_kos (dict[str, set[str]]): Sub-threshold candidate KOs from each tool.
        rescued_kos (dict[str, set[str]]): Heuristic rescued KOs from KOfam.
        eggnog_filter_multi (str, optional): Strategy ('disambiguate' or 'none'). Defaults to 'disambiguate'.

    Returns:
        tuple[set[str], set[str], set[Any], set[Any]]: Retained KOs, globally
        dropped KOs, narrowed source-group IDs, and unresolved source-group IDs.
    """
    en_dropped = set()
    narrowed_groups: set[Any] = set()
    unresolved_groups: set[Any] = set()
    if "eggnog" in active_tools and en_agreed:
        groups: dict[Any, set[str]] = {}
        singleton_support: set[str] = set()
        for r in tool_recs.get("eggnog", []):
            sources = r.get("eggnog_sources") or [{
                "group": r.get("eggnog_seed_group", id(r)),
                "is_multi": r.get("eggnog_shared_seed_hit", False),
                "status": r["original_status"],
            }]
            for source in sources:
                if source["status"] == "threshold_passing" and source["is_multi"]:
                    groups.setdefault(source["group"], set()).add(r["ko"])
                elif source["status"] == "threshold_passing":
                    singleton_support.add(r["ko"])

        grouped_kos = set().union(*groups.values()) if groups else set()
        retained = (en_agreed - grouped_kos) | singleton_support
        trusted_other = set().union(
            *[confident_kos.get(t, set()) for t in active_tools if t != "eggnog"]
        )
        candidate_other = set().union(
            *[
                rescued_kos.get(t, set()) | candidate_kos.get(t, set())
                for t in active_tools
                if t != "eggnog"
            ]
        )

        for group_id, group_kos in groups.items():
            kept = set(group_kos)
            if eggnog_filter_multi == "disambiguate" and len(group_kos) > 1:
                overlap = group_kos & trusted_other
                if not overlap:
                    overlap = group_kos & candidate_other
                if overlap:
                    kept = overlap
                    narrowed_groups.add(group_id)
                else:
                    unresolved_groups.add(group_id)
            retained |= kept

        en_agreed = en_agreed & retained
        # A KO is globally dropped only when no independent confident source survives.
        en_dropped = set(confident_kos.get("eggnog", set())) - en_agreed

    return en_agreed, en_dropped, narrowed_groups, unresolved_groups


def adjudicate_gene(
    active_tools: list[str],
    confident_kos: dict[str, set[str]],
    candidate_kos: dict[str, set[str]],
    rescued_kos: dict[str, set[str]],
    en_dropped: set[str],
    conflict_strategy: str = "multiple",
    priority_order: Optional[list[str]] = None,
) -> tuple[str, set[str], str, str]:
    """Adjudicate consensus KO, alternative KOs, consensus level, and evidence string for a gene.

    Evaluates confident (threshold-passing), candidate (below-threshold), and
    heuristic-rescued calls across active tools according to specified hierarchical
    resolution rules.

    Alternative KOs include confident minority calls, eggNOG KOs dropped during
    multi-KO disambiguation, and corroborated candidate KOs when multiple candidates
    remain unresolved. Uncorroborated below-threshold candidates are never injected.

    Args:
        active_tools: List of active annotation tools (e.g. ['kofam', 'deepkoala', 'eggnog']).
        confident_kos: Mapping of tool name to set of confident (threshold-passing) KO IDs.
        candidate_kos: Mapping of tool name to set of candidate (below-threshold) KO IDs.
        rescued_kos: Mapping of tool name to set of heuristic-rescued KO IDs (specifically KOfam).
        en_dropped: Set of KO IDs dropped from eggNOG multi-KO predictions during disambiguation.
        conflict_strategy: Resolution strategy for disjoint conflicts ('multiple' or 'priority').
        priority_order: Ordered list of tools for priority-based conflict resolution.

    Returns:
        tuple[str, set[str], str, str]: (accepted_ko, alt_set, consensus_level, evidence_str)
            - accepted_ko: Primary consensus KO ID, or '-' if none/conflicting/ambiguous.
            - alt_set: Set of alternative, conflicting, or dropped KO IDs.
            - consensus_level: Classification label ('unanimous', 'majority', 'single_tool',
              'single_tool_with_candidate', 'dual_candidate', 'conflict',
              'conflict_priority', 'unannotated').
            - evidence_str: Comma-separated string describing calling or supporting tools.
    """
    # Filter to tools that made confident (threshold-passing) calls for this gene
    tool_calls = {t: confident_kos[t] for t in active_tools if t in confident_kos}
    num_calling_tools = len(tool_calls)

    accepted_ko = "-"
    alt_set = set()
    consensus_level = "unannotated"
    evidence_str = "-"

    # -------------------------------------------------------------------------
    # CASE 1: No tool produced a confident call (num_calling_tools == 0)
    # Check if candidate or heuristic-rescued hits can establish orthogonal consensus.
    # -------------------------------------------------------------------------
    if num_calling_tools == 0:
        dk_cand = candidate_kos.get("deepkoala", set())
        en_cand = candidate_kos.get("eggnog", set())
        kf_resc = rescued_kos.get("kofam", set())

        # Zero-calling-tools decision hierarchy:
        # 1. KOfam rescue in any form first:
        #    a. KOfam rescue agreed with both DeepKOALA and eggNOG candidates (3-way)
        #    b. KOfam rescue agreed with eggNOG candidate
        #    c. KOfam rescue agreed with DeepKOALA candidate
        #    d. KOfam rescue by itself (single_tool)
        # 2. If no KOfam rescue, test DeepKOALA candidate agreeing with eggNOG candidate
        # 3. Otherwise: remain unannotated.
        cand_agree = set()
        agreeing_cands = []

        if "kofam" in active_tools and kf_resc:
            # Check for candidate agreement with KOfam rescue across tools
            if "deepkoala" in active_tools and "eggnog" in active_tools and (kf_resc & en_cand & dk_cand):
                # 3-way candidate consensus (KOfam rescued + DeepKOALA candidate + eggNOG candidate)
                cand_agree = kf_resc & en_cand & dk_cand
                agreeing_cands = ["deepkoala(candidate)", "eggnog(candidate)", "kofam(rescued)"]
            elif "eggnog" in active_tools and (kf_resc & en_cand):
                # 2-way candidate consensus (KOfam rescued + eggNOG candidate)
                cand_agree = kf_resc & en_cand
                agreeing_cands = ["eggnog(candidate)", "kofam(rescued)"]
            elif "deepkoala" in active_tools and (kf_resc & dk_cand):
                # 2-way candidate consensus (KOfam rescued + DeepKOALA candidate)
                cand_agree = kf_resc & dk_cand
                agreeing_cands = ["deepkoala(candidate)", "kofam(rescued)"]

            if cand_agree:
                consensus_level = "dual_candidate"
                evidence_str = ",".join(sorted(agreeing_cands))
                if len(cand_agree) == 1:
                    # Exactly one unanimous candidate KO accepted
                    accepted_ko = next(iter(cand_agree))
                    alt_set = en_dropped - {accepted_ko}
                else:
                    # Multiple candidate KOs agreed; remain conservative and assign to alternatives
                    accepted_ko = "-"
                    alt_set = cand_agree | en_dropped
            else:
                # KOfam rescue accepted by itself without candidate corroboration
                consensus_level = "single_tool"
                evidence_str = "kofam(rescued)"
                if len(kf_resc) == 1:
                    accepted_ko = next(iter(kf_resc))
                    alt_set = en_dropped - {accepted_ko}
                else:
                    accepted_ko = "-"
                    alt_set = kf_resc | en_dropped

        elif "deepkoala" in active_tools and "eggnog" in active_tools and (dk_cand & en_cand):
            # Orthogonal agreement between below-threshold DeepKOALA and eggNOG candidates
            cand_agree = dk_cand & en_cand
            agreeing_cands = ["deepkoala(candidate)", "eggnog(candidate)"]
            consensus_level = "dual_candidate"
            evidence_str = ",".join(sorted(agreeing_cands))
            if len(cand_agree) == 1:
                accepted_ko = next(iter(cand_agree))
                alt_set = en_dropped - {accepted_ko}
            else:
                accepted_ko = "-"
                alt_set = cand_agree | en_dropped

        else:
            # No confident or corroborated candidate calls
            accepted_ko = "-"
            alt_set = set()
            consensus_level = "unannotated"
            evidence_str = "-"

    # -------------------------------------------------------------------------
    # CASE 2: Exactly one tool produced a confident call (num_calling_tools == 1)
    # Check if other non-calling tools have sub-threshold candidates supporting it.
    # -------------------------------------------------------------------------
    elif num_calling_tools == 1:
        tool_name = next(iter(tool_calls))
        called_kos = tool_calls[tool_name]

        dk_cand = candidate_kos.get("deepkoala", set())
        en_cand = candidate_kos.get("eggnog", set())
        kf_resc = rescued_kos.get("kofam", set())

        # Check for sub-threshold candidate corroboration from other active tools
        rescued_by = []
        if tool_name != "deepkoala" and "deepkoala" in active_tools and (called_kos & dk_cand):
            rescued_by.append("deepkoala(candidate)")
        if tool_name != "eggnog" and "eggnog" in active_tools and (called_kos & en_cand):
            rescued_by.append("eggnog(candidate)")
        if tool_name != "kofam" and "kofam" in active_tools and (called_kos & kf_resc):
            rescued_by.append("kofam(rescued)")

        # ... keep 'en_dropped' as the single source of disambiguated eggNOG alternatives
        if len(called_kos) == 1:
            # Single unambiguous confident call
            accepted_ko = next(iter(called_kos))
            alt_set = en_dropped - {accepted_ko}
        else:
            # Multiple confident KOs from the same tool; keep accepted_ko as '-' to avoid ambiguity
            accepted_ko = "-"
            alt_set = called_kos | en_dropped

        if rescued_by:
            # Confident call supported by sub-threshold candidate from another tool
            consensus_level = "single_tool_with_candidate"
            evidence_str = ",".join([tool_name] + rescued_by)
        else:
            # Uncorroborated single-tool confident call
            consensus_level = "single_tool"
            evidence_str = tool_name

    # -------------------------------------------------------------------------
    # CASE 3: Two or more tools produced confident calls (num_calling_tools >= 2)
    # Evaluate intersection (unanimous or majority agreement) or resolve conflicts.
    # -------------------------------------------------------------------------
    else:
        # Extract sets of called KOs across all calling tools
        all_sets = list(tool_calls.values())
        union_all_calls = set.union(*all_sets)
        common_all = set.intersection(*all_sets)

        if common_all:
            # Complete agreement across all calling tools
            is_unanimous = (num_calling_tools == len(active_tools))
            consensus_level = "unanimous" if is_unanimous else "majority"
            evidence_str = ",".join(sorted(tool_calls.keys()))
            unselected_minority = union_all_calls - common_all
            extra_alts = unselected_minority | en_dropped

            if len(common_all) == 1:
                accepted_ko = next(iter(common_all))
                alt_set = extra_alts - {accepted_ko}
            else:
                accepted_ko = "-"
                alt_set = common_all | extra_alts
        else:
            # Check for pairwise majority agreement between pairs of calling tools (e.g. 2 out of 3)
            pairwise_agreed = set()
            agreeing_tools = set()
            tools_list = list(tool_calls.keys())
            for i in range(len(tools_list)):
                for j in range(i + 1, len(tools_list)):
                    inter = tool_calls[tools_list[i]] & tool_calls[tools_list[j]]
                    if inter:
                        pairwise_agreed.update(inter)
                        agreeing_tools.add(tools_list[i])
                        agreeing_tools.add(tools_list[j])

            if pairwise_agreed:
                # Majority consensus reached among a subset of calling tools
                consensus_level = "majority"
                evidence_str = ",".join(sorted(agreeing_tools))
                unselected_minority = union_all_calls - pairwise_agreed
                extra_alts = unselected_minority | en_dropped

                if len(pairwise_agreed) == 1:
                    accepted_ko = next(iter(pairwise_agreed))
                    alt_set = extra_alts - {accepted_ko}
                else:
                    accepted_ko = "-"
                    alt_set = pairwise_agreed | extra_alts
            else:
                # Disjoint conflict: calling tools produced completely non-overlapping KO predictions
                all_alts = union_all_calls | en_dropped
                if conflict_strategy == "priority":
                    # Priority resolution: pick top-ranking tool according to predefined tool hierarchy
                    if priority_order is None:
                        priority_order = [t for t in ["kofam", "eggnog", "deepkoala"] if t in active_tools]
                    top_tool = next((t for t in priority_order if t in tool_calls), tools_list[0])
                    top_kos = tool_calls[top_tool]
                    if len(top_kos) == 1:
                        accepted_ko = next(iter(top_kos))
                        alt_set = ((union_all_calls - top_kos) | en_dropped) - {accepted_ko}
                    else:
                        accepted_ko = "-"
                        alt_set = all_alts
                    consensus_level = "conflict_priority"
                    evidence_str = top_tool
                else:  # "multiple" (default conservative strategy)
                    # Keep accepted_ko as '-' and place all conflicting calls into alternative_kos
                    # to avoid spurious multifunctional enzyme assertions downstream
                    accepted_ko = "-"
                    alt_set = all_alts
                    consensus_level = "conflict"
                    evidence_str = ",".join(sorted(tool_calls.keys()))

    return accepted_ko, alt_set, consensus_level, evidence_str


def build_gene_summary_row(
    gid: str,
    accepted_ko: str,
    alt_set: set[str],
    consensus_level: str,
    evidence_str: str,
    tool_recs: dict[str, list[dict]],
    en_agreed: set[str],
    master_ko_defs: dict[str, str],
    evidence_defs: dict[str, str],
) -> dict[str, Any]:
    """Build summary metrics and annotations dictionary for a single gene.

    Formats consensus calls, functional definitions, and method-specific scores
    (KOfam, DeepKOALA, eggNOG) into a single row dictionary matching the final
    gene summary table schema.

    Args:
        gid: Gene identifier.
        accepted_ko: Primary accepted consensus KO ID, or '-' if none/ambiguous.
        alt_set: Set of alternative, conflicting, or dropped KO IDs.
        consensus_level: Gene consensus classification label.
        evidence_str: Comma-separated string describing calling or supporting tools.
        tool_recs: Mapping of tool name to raw evidence records for this gene.
        en_agreed: Retained eggNOG KO IDs after multi-KO disambiguation.
        master_ko_defs: Master dictionary mapping KO IDs to definitions from ko_list.
        evidence_defs: Fallback dictionary mapping KO IDs to definitions from tool outputs.

    Returns:
        dict[str, Any]: Formatted row dictionary ready for DataFrame assembly.
    """
    alternative_kos = format_kos(alt_set)
    accepted_set = {accepted_ko} if accepted_ko != "-" else set()
    # Resolve definitions for consensus accepted KO and alternative KOs
    definition = resolve_definition_for_kos(accepted_set, master_ko_defs, evidence_defs)
    alternative_definition = resolve_definition_for_kos(alt_set, master_ko_defs, evidence_defs)

    # -------------------------------------------------------------------------
    # Format KOfam metrics
    # If single record: format scalar values.
    # If multiple records: format as comma-separated 'KO:value' strings.
    # -------------------------------------------------------------------------
    # External tables can contain weak KOfam records; keep those in evidence but do
    # not present them as confident KOfam summary assignments.
    kf_recs = [
        r for r in tool_recs.get("kofam", [])
        if r["original_status"] in ("threshold_passing", "heuristic_rescued")
    ]
    if kf_recs:
        kf_kos = sorted({r["ko"] for r in kf_recs})
        kofam_ko_val = format_kos(set(kf_kos))
        if len(kf_recs) == 1:
            r0 = kf_recs[0]
            sel_bs, sel_ev, norm_st = select_kofam_scores(
                r0.get("score_type"),
                r0.get("bit_score"),
                r0.get("e_value"),
                r0.get("domain_bit_score"),
                r0.get("domain_e_value"),
            )
            kf_st_val = norm_st
            kf_assign_val = str(r0["assignment"])
            kf_bs_val = sel_bs
            kf_ev_val = sel_ev
            kf_th_val = r0["threshold"]
        else:
            # Multi-hit KOfam output: format score types ordered by formatted KO list
            kos_order = kofam_ko_val.split(",")
            ko_to_rec = {r["ko"]: r for r in kf_recs}
            st_list = []
            for k in kos_order:
                if k in ko_to_rec:
                    r = ko_to_rec[k]
                    _, _, norm_st = select_kofam_scores(
                        r.get("score_type"),
                        r.get("bit_score"),
                        r.get("e_value"),
                        r.get("domain_bit_score"),
                        r.get("domain_e_value"),
                    )
                    st_list.append(norm_st)
            kf_st_val = ",".join(st_list)
            # Multi-hit metrics: prefix with KO ID to keep associations unambiguous
            kf_assign_val = ",".join(f"{r['ko']}:{r['assignment']}" for r in kf_recs)
            kf_bs_val = ",".join(
                f"{r['ko']}:{format_metric_value(select_kofam_scores(r.get('score_type'), r.get('bit_score'), r.get('e_value'), r.get('domain_bit_score'), r.get('domain_e_value'))[0])}"
                for r in kf_recs
            )
            kf_ev_val = ",".join(
                f"{r['ko']}:{format_metric_value(select_kofam_scores(r.get('score_type'), r.get('bit_score'), r.get('e_value'), r.get('domain_bit_score'), r.get('domain_e_value'))[1])}"
                for r in kf_recs
            )
            kf_th_val = ",".join(f"{r['ko']}:{format_metric_value(r['threshold'])}" for r in kf_recs)
    else:
        kofam_ko_val = "-"
        kf_st_val = "-"
        kf_assign_val = "-"
        kf_bs_val = np.nan
        kf_ev_val = np.nan
        kf_th_val = np.nan

    # -------------------------------------------------------------------------
    # Format DeepKOALA metrics
    # Distinguish confident (threshold_passing) from candidate (below_threshold).
    # Preserve KO-to-metric associations when more than one record exists.
    # -------------------------------------------------------------------------
    dk_recs = tool_recs.get("deepkoala", [])
    if dk_recs:
        dk_thresh_kos = {r["ko"] for r in dk_recs if r["original_status"] == "threshold_passing"}
        dk_cands_kos = {r["ko"] for r in dk_recs if r["original_status"] in ("threshold_passing", "below_threshold")}
        dk_ko_val = format_kos(dk_thresh_kos)
        dk_cand_val = format_kos(dk_cands_kos)
        ordered_dk = sorted(
            (r for r in dk_recs if r["original_status"] in ("threshold_passing", "below_threshold")),
            key=lambda r: r["ko"],
        )
        if not ordered_dk:
            dk_score_val = np.nan
            dk_th_val = np.nan
        elif len(ordered_dk) == 1:
            dk_score_val = ordered_dk[0]["deepkoala_probability"]
            dk_th_val = ordered_dk[0]["deepkoala_threshold"]
        else:
            dk_score_val = ",".join(
                f"{r['ko']}:{format_metric_value(r['deepkoala_probability'])}" for r in ordered_dk
            )
            dk_th_val = ",".join(
                f"{r['ko']}:{format_metric_value(r['deepkoala_threshold'])}" for r in ordered_dk
            )
    else:
        dk_ko_val = "-"
        dk_cand_val = "-"
        dk_score_val = np.nan
        dk_th_val = np.nan

    # -------------------------------------------------------------------------
    # Format eggNOG metrics
    # Report disambiguated KOs in eggnog_ko and all candidates in eggnog_candidate_ko.
    # Preserve KO-to-metric associations when more than one record exists.
    # -------------------------------------------------------------------------
    en_recs = tool_recs.get("eggnog", [])
    if en_recs:
        en_cands_kos = {r["ko"] for r in en_recs if r["original_status"] in ("threshold_passing", "below_threshold")}
        en_ko_val = format_kos(en_agreed)
        en_cand_val = format_kos(en_cands_kos)
        ordered_en = sorted(
            (r for r in en_recs if r["original_status"] in ("threshold_passing", "below_threshold")),
            key=lambda r: r["ko"],
        )
        if not ordered_en:
            en_bs_val = np.nan
            en_ev_val = np.nan
        elif len(ordered_en) == 1:
            en_bs_val = ordered_en[0]["bit_score"]
            en_ev_val = ordered_en[0]["e_value"]
        else:
            en_bs_val = ",".join(f"{r['ko']}:{format_metric_value(r['bit_score'])}" for r in ordered_en)
            en_ev_val = ",".join(f"{r['ko']}:{format_metric_value(r['e_value'])}" for r in ordered_en)
    else:
        en_ko_val = "-"
        en_cand_val = "-"
        en_bs_val = np.nan
        en_ev_val = np.nan

    return {
        "gene_id": gid,
        "accepted_ko": accepted_ko,
        "alternative_kos": alternative_kos,
        "definition": definition,
        "alternative_definition": alternative_definition,
        "consensus_level": consensus_level,
        "evidence": evidence_str,
        "kofam_ko": kofam_ko_val,
        "kofam_score_type": kf_st_val,
        "kofam_bit_score": kf_bs_val,
        "kofam_evalue": kf_ev_val,
        "kofam_assignment": kf_assign_val,
        "kofam_threshold": kf_th_val,
        "deepkoala_ko": dk_ko_val,
        "deepkoala_candidate_ko": dk_cand_val,
        "deepkoala_score": dk_score_val,
        "deepkoala_threshold": dk_th_val,
        "eggnog_ko": en_ko_val,
        "eggnog_candidate_ko": en_cand_val,
        "eggnog_bit_score": en_bs_val,
        "eggnog_evalue": en_ev_val,
    }


def annotate_evidence_records(
    gid: str,
    recs: list[dict],
    tool_recs: dict[str, list[dict]],
    active_tools: list[str],
    accepted_set: set[str],
    alt_set: set[str],
    accepted_ko: str,
    consensus_level: str,
    en_dropped: set[str],
    narrowed_eggnog_groups: set[Any],
    unresolved_eggnog_groups: set[Any],
) -> list[dict]:
    """Annotate raw evidence records for a gene with cross-method consensus and support.

    Computes supporting methods and final cross-method status for each evidence
    record based on gene-level adjudication outcomes.

    Args:
        gid: Gene identifier.
        recs: Raw evidence records for this gene across all methods.
        tool_recs: Mapping of tool name to evidence records for this gene.
        active_tools: List of active annotation tools.
        accepted_set: Set containing the accepted KO ID (or empty if none).
        alt_set: Set of alternative/conflicting KO IDs.
        accepted_ko: Primary accepted KO ID string ('-' if none).
        consensus_level: Gene consensus classification label.
        en_dropped: Set of eggNOG KO IDs dropped during multi-KO disambiguation.
        narrowed_eggnog_groups: eggNOG source groups narrowed during disambiguation.
        unresolved_eggnog_groups: eggNOG multi-KO groups lacking corroboration.

    Returns:
        list[dict]: Adjudicated evidence records populated with 'cross_method_status'
            and 'supporting_methods'.
    """
    adjudicated = []
    # eggNOG seed-hit groups that actually lost a KO during multi-KO disambiguation.
    # A retained eggNOG KO is only labeled 'disambiguated_retained' when its own
    # multi-KO seed group was genuinely narrowed (i.e. it shares a group with a
    # dropped counterpart). With eggnog_filter_multi='none' no group ever drops a KO,
    # so the label cannot occur.
    for r in recs:
        ko = r["ko"]
        method = r["method"]
        orig_status = r["original_status"]

        # ---------------------------------------------------------------------
        # Determine supporting methods:
        # Check other active tools for matching KO predictions.
        # Record all independent matching evidence. Retention is checked separately
        # before a below-threshold record can be called rescued.
        # ---------------------------------------------------------------------
        supporting = []
        for other_t in active_tools:
            if other_t != method:
                other_t_recs = tool_recs.get(other_t, [])
                other_matching = [
                    o for o in other_t_recs
                    if o["ko"] == ko and o["original_status"] in ("threshold_passing", "heuristic_rescued", "below_threshold")
                ]
                if other_matching:
                    supporting.append(other_t)

        supp_str = ",".join(sorted(supporting)) if supporting else "-"

        # ---------------------------------------------------------------------
        # Determine cross_method_status state machine:
        # - unknown: missing/uninterpretable status
        # - disambiguated_dropped: eggNOG multi-KO dropped during disambiguation
        # - accepted / disambiguated_retained: hit matches accepted consensus KO
        # - heuristic_rescued: KOfam hit rescued by score threshold margin
        # - rescued: below-threshold candidate elevated by orthogonal consensus
        # - conflict: threshold-passing hit that conflicted under disjoint calls
        # - unresolved_multi_ko: eggNOG multi-KO hit that could not be disambiguated
        # - alternative: threshold-passing call not chosen as primary consensus
        # - below_threshold_unrescued: candidate hit without sufficient support
        # - unselected: threshold-passing or heuristic-rescued hit not selected in consensus
        # ---------------------------------------------------------------------
        if orig_status == "unknown":
            cross_status = "unknown"
        elif method == "eggnog" and ko in en_dropped:
            cross_status = "disambiguated_dropped"
        elif ko in accepted_set:
            if orig_status == "threshold_passing":
                if (
                    method == "eggnog"
                    and any(
                        source["group"] in narrowed_eggnog_groups and source["status"] == "threshold_passing"
                        for source in r.get("eggnog_sources", [])
                    )
                ):
                    cross_status = "disambiguated_retained"
                else:
                    cross_status = "accepted"
            elif orig_status == "heuristic_rescued":
                cross_status = "heuristic_rescued"
            elif orig_status == "below_threshold":
                cross_status = "rescued"
            else:
                cross_status = "unknown"
        elif ko in alt_set:
            if orig_status == "threshold_passing":
                if consensus_level.startswith("conflict"):
                    cross_status = "conflict"
                elif method == "eggnog" and any(
                    source["group"] in unresolved_eggnog_groups
                    and source["status"] == "threshold_passing"
                    for source in r.get("eggnog_sources", [])
                ):
                    cross_status = "unresolved_multi_ko"
                else:
                    cross_status = "alternative"
            elif orig_status == "heuristic_rescued":
                cross_status = "heuristic_rescued"
            elif orig_status == "below_threshold":
                cross_status = "rescued" if supporting else "below_threshold_unrescued"
            else:
                cross_status = "unknown"
        else:
            if orig_status == "below_threshold":
                cross_status = "below_threshold_unrescued"
            elif orig_status == "threshold_passing":
                cross_status = "unselected"
            elif orig_status == "heuristic_rescued":
                cross_status = "unselected"  # rescue lost to another tool's call; rescue preserved in original_status
            else:
                cross_status = "unselected"

        adjudicated.append({
            "gene_id": gid,
            "ko": ko,
            "method": method,
            "consensus_level": consensus_level,
            "original_status": orig_status,
            "bit_score": r.get("bit_score", np.nan),
            "e_value": r.get("e_value", np.nan),
            "domain_bit_score": r.get("domain_bit_score", np.nan),
            "domain_e_value": r.get("domain_e_value", np.nan),
            "score_type": r.get("score_type", "-"),
            "bit_score_threshold": r.get("bit_score_threshold", r.get("threshold", np.nan)),
            "e_value_threshold": r.get("e_value_threshold", np.nan),
            "deepkoala_probability": r.get("deepkoala_probability", np.nan),
            "deepkoala_threshold": r.get("deepkoala_threshold", np.nan),
            "eggnog_shared_seed_hit": bool(r.get("eggnog_shared_seed_hit", r.get("shared_seed_hit", False))),
            "cross_method_status": cross_status,
            "supporting_methods": supp_str,
        })
    return adjudicated


def export_tables(
    result_df: pd.DataFrame,
    evidence_df: pd.DataFrame,
    output_tsv: Optional[Union[str, Path]] = None,
    evidence_tsv: Optional[Union[str, Path]] = None,
) -> None:
    """Export gene summary and evidence audit tables to TSV format.

    Ensures parent directories exist, formats float NaNs as '-', booleans as
    'True'/'False', and writes clean tab-separated files.

    Args:
        result_df: Gene-level integrated consensus DataFrame.
        evidence_df: Long-form evidence audit DataFrame.
        output_tsv: Destination file path for gene summary table.
        evidence_tsv: Destination file path for evidence table. If None and
            output_tsv is specified, defaults to 'kolach_evidence.tsv' in the
            same directory as output_tsv.

    Returns:
        None
    """
    if output_tsv:
        out_path = Path(output_tsv).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Create a formatted copy so internal DataFrame dtypes remain unaffected
        tsv_export = result_df.copy()
        for col in tsv_export.columns:
            # Replace NaN in float columns with '-' string while preserving numeric string representation
            if pd.api.types.is_float_dtype(tsv_export[col]):
                tsv_export[col] = tsv_export[col].apply(lambda x: "-" if pd.isna(x) else str(x))
            else:
                tsv_export[col] = tsv_export[col].fillna("-")
        tsv_export.to_csv(out_path, sep="\t", index=False)

        # Automatically export evidence TSV alongside annotations if not separately specified
        if evidence_tsv is None:
            evidence_tsv = out_path.parent / "kolach_evidence.tsv"

    if evidence_tsv:
        ev_path = Path(evidence_tsv).expanduser().resolve()
        ev_path.parent.mkdir(parents=True, exist_ok=True)
        # Create a formatted copy for evidence export
        ev_export = evidence_df.copy()
        for col in ev_export.columns:
            # Handle float NaNs, booleans, and general missing values
            if pd.api.types.is_float_dtype(ev_export[col]):
                ev_export[col] = ev_export[col].apply(lambda x: "-" if pd.isna(x) else str(x))
            elif pd.api.types.is_bool_dtype(ev_export[col]):
                ev_export[col] = ev_export[col].apply(lambda x: "True" if x else "False")
            else:
                ev_export[col] = ev_export[col].fillna("-")
        ev_export.to_csv(ev_path, sep="\t", index=False)


def integrate_annotations(
    protein_fasta: Optional[Union[str, Path]] = None,
    kofam_tsv: Optional[Union[str, Path]] = None,
    deepkoala_tsv: Optional[Union[str, Path]] = None,
    eggnog_tsv: Optional[Union[str, Path]] = None,
    output_tsv: Optional[Union[str, Path]] = None,
    evidence_tsv: Optional[Union[str, Path]] = None,
    database_dir: Optional[Union[str, Path]] = None,
    eggnog_min_bitscore: float = 60.0,
    eggnog_max_evalue: float = 1e-5,
    eggnog_filter_multi: str = "disambiguate",
    conflict_strategy: str = "multiple",
    heuristic_bitscore_fraction: float = 0.75,
    heuristic_e_value: float = 1e-5,
) -> pd.DataFrame:
    """Main evidence-first integration pipeline for kolach annotation outputs.

    Parses each tool output once into normalized per-gene, per-KO evidence records.
    Directly adjudicates filtering, eggNOG disambiguation, rescue, and consensus
    from those records, ensuring 1-to-1 consistency between the gene summary table
    and the long-form evidence table.

    Args:
        protein_fasta: Path to query protein FASTA. When provided, defines the
            complete universe and ordering of gene IDs (including unannotated genes).
        kofam_tsv: Path to KOfam tabular output file.
        deepkoala_tsv: Path to DeepKOALA tabular output file.
        eggnog_tsv: Path to eggNOG-mapper tabular output file.
        output_tsv: Optional file path to export integrated gene summary TSV.
        evidence_tsv: Optional file path to export long-form evidence TSV.
        database_dir: Directory containing KEGG/KOfam database (e.g. ko_list).
        eggnog_min_bitscore: Minimum bit score threshold for eggNOG hits (default: 60.0).
        eggnog_max_evalue: Maximum e-value threshold for eggNOG hits (default: 1e-5).
        eggnog_filter_multi: Strategy for disambiguating multi-KO eggNOG seed hits
            ('disambiguate' or 'none').
        conflict_strategy: Resolution strategy for disjoint calling tools
            ('multiple' or 'priority').
        heuristic_bitscore_fraction: Fraction of the KOfam profile threshold used
            to rescue below-threshold hits (default: 0.75). Reported in the
            evidence table for heuristic_rescued hits.
        heuristic_e_value: Maximum E-value for KOfam rescue consideration
            (default: 1e-5). Reported in the evidence table for heuristic_rescued hits.

    Returns:
        pd.DataFrame: Integrated gene-level summary table with the long-form evidence
            audit DataFrame stored in result_df.attrs['evidence'].
    """
    if conflict_strategy not in {"multiple", "priority"}:
        raise ValueError("conflict_strategy must be 'multiple' or 'priority'")
    if eggnog_filter_multi not in {"disambiguate", "none"}:
        raise ValueError("eggnog_filter_multi must be 'disambiguate' or 'none'")

    numeric_parameters = {
        "eggnog_min_bitscore": eggnog_min_bitscore,
        "eggnog_max_evalue": eggnog_max_evalue,
        "heuristic_bitscore_fraction": heuristic_bitscore_fraction,
        "heuristic_e_value": heuristic_e_value,
    }
    for name, value in numeric_parameters.items():
        try:
            finite = np.isfinite(value)
        except TypeError as exc:
            raise ValueError(f"{name} must be a finite number") from exc
        if not finite:
            raise ValueError(f"{name} must be a finite number")
    if eggnog_min_bitscore < 0:
        raise ValueError("eggnog_min_bitscore must be nonnegative")
    if eggnog_max_evalue < 0 or heuristic_e_value < 0:
        raise ValueError("E-value cutoffs must be nonnegative")
    if not 0 <= heuristic_bitscore_fraction <= 1:
        raise ValueError("heuristic_bitscore_fraction must be within [0, 1]")

    active_tools = []
    all_evidence_records: list[dict] = []
    evidence_defs: dict[str, str] = {}  # KO-specific definitions harvested from KOfam

    # Validate that explicitly supplied input files exist before processing
    explicit_inputs = [
        ("protein_fasta", protein_fasta),
        ("kofam_tsv", kofam_tsv),
        ("deepkoala_tsv", deepkoala_tsv),
        ("eggnog_tsv", eggnog_tsv),
    ]
    for arg_name, p in explicit_inputs:
        if p is not None:
            norm_p = Path(p).expanduser().resolve()
            if not norm_p.is_file():
                raise FileNotFoundError(
                    f"Supplied {arg_name} file does not exist or is not a regular file: '{p}' (resolved: '{norm_p}')"
                )

    # 1. Ingest KOfam: parse hits, select scores, and harvest KO definitions
    if kofam_tsv is not None:
        kf_records = extract_kofam_records(
            kofam_tsv,
            heuristic_bitscore_fraction=heuristic_bitscore_fraction,
            heuristic_e_value=heuristic_e_value,
        )
        all_evidence_records.extend(kf_records)
        for r in kf_records:
            if r.get("definition") and r["definition"] != "-":
                evidence_defs[r["ko"]] = r["definition"]
        active_tools.append("kofam")

    # 2. Ingest DeepKOALA: parse predictions, handle 2-col vs 3-col formats
    if deepkoala_tsv is not None:
        dk_records = extract_deepkoala_records(deepkoala_tsv)
        all_evidence_records.extend(dk_records)
        active_tools.append("deepkoala")

    # 3. Ingest eggNOG: parse hits, filter on bit score / e-value, tag multi-KO hits
    if eggnog_tsv is not None:
        en_records = extract_eggnog_records(
            eggnog_tsv,
            min_bitscore=eggnog_min_bitscore,
            max_evalue=eggnog_max_evalue,
        )
        all_evidence_records.extend(en_records)
        active_tools.append("eggnog")

    # Determine complete universe and order of gene IDs (FASTA order takes precedence)
    if protein_fasta is not None:
        all_gene_ids = read_fasta_ids(protein_fasta)
    else:
        # Preserve first-observed order across evidence records
        seen = set()
        ordered_ids = []
        for r in all_evidence_records:
            gid = r["gene_id"]
            if gid not in seen:
                seen.add(gid)
                ordered_ids.append(gid)
        all_gene_ids = ordered_ids

    # Load master KO definitions from ko_list if available in database directory
    ko_list_file = None
    if database_dir:
        db_path = Path(database_dir).expanduser().resolve()
        if (db_path / "kofam" / "ko_list").is_file():
            ko_list_file = db_path / "kofam" / "ko_list"
        elif (db_path / "ko_list").is_file():
            ko_list_file = db_path / "ko_list"
    master_ko_defs = load_ko_definitions(ko_list_file)

    # Group all ingested evidence records by gene_id
    records_by_gene: dict[str, list[dict]] = {gid: [] for gid in all_gene_ids}
    for r in all_evidence_records:
        gid = r["gene_id"]
        if gid not in records_by_gene:
            # Include any gene IDs found in evidence tables even if omitted from FASTA
            records_by_gene[gid] = []
            all_gene_ids.append(gid)
        records_by_gene[gid].append(r)

    # Adjudicate consensus and cross-method evidence per gene directly from records
    gene_summary_rows = []
    adjudicated_records = []

    priority_order = [t for t in ["kofam", "eggnog", "deepkoala"] if t in active_tools]

    for gid in all_gene_ids:
        recs = records_by_gene.get(gid, [])

        # Categorize evidence records by tool
        tool_recs: dict[str, list[dict]] = {t: [] for t in active_tools}
        for r in recs:
            t = r["method"]
            if t in tool_recs:
                tool_recs[t].append(r)

        # Categorize calls into confident, candidate, and heuristic-rescued sets
        confident_kos: dict[str, set[str]] = {}
        candidate_kos: dict[str, set[str]] = {}
        rescued_kos: dict[str, set[str]] = {}

        for t in active_tools:
            t_list = tool_recs[t]
            thresh = {r["ko"] for r in t_list if r["original_status"] == "threshold_passing"}
            cands = {r["ko"] for r in t_list if r["original_status"] == "below_threshold"}
            resc = {r["ko"] for r in t_list if r["original_status"] == "heuristic_rescued"}

            if thresh:
                confident_kos[t] = thresh
            if cands:
                candidate_kos[t] = cands
            if resc:
                rescued_kos[t] = resc

        # eggNOG multi-KO disambiguation against other tools
        en_agreed = set(confident_kos.get("eggnog", set()))
        en_agreed, en_dropped, narrowed_eggnog_groups, unresolved_eggnog_groups = disambiguate_eggnog(
            en_agreed=en_agreed,
            tool_recs=tool_recs,
            active_tools=active_tools,
            confident_kos=confident_kos,
            candidate_kos=candidate_kos,
            rescued_kos=rescued_kos,
            eggnog_filter_multi=eggnog_filter_multi,
        )
        if en_agreed:
            confident_kos["eggnog"] = en_agreed

        # Adjudicate gene-level consensus
        accepted_ko, alt_set, consensus_level, evidence_str = adjudicate_gene(
            active_tools=active_tools,
            confident_kos=confident_kos,
            candidate_kos=candidate_kos,
            rescued_kos=rescued_kos,
            en_dropped=en_dropped,
            conflict_strategy=conflict_strategy,
            priority_order=priority_order,
        )

        # Build gene-level summary row
        summary_row = build_gene_summary_row(
            gid=gid,
            accepted_ko=accepted_ko,
            alt_set=alt_set,
            consensus_level=consensus_level,
            evidence_str=evidence_str,
            tool_recs=tool_recs,
            en_agreed=en_agreed,
            master_ko_defs=master_ko_defs,
            evidence_defs=evidence_defs,
        )
        gene_summary_rows.append(summary_row)

        # Annotate evidence records with cross-method status and supporting tools
        accepted_set = {accepted_ko} if accepted_ko != "-" else set()
        gene_adjudicated = annotate_evidence_records(
            gid=gid,
            recs=recs,
            tool_recs=tool_recs,
            active_tools=active_tools,
            accepted_set=accepted_set,
            alt_set=alt_set,
            accepted_ko=accepted_ko,
            consensus_level=consensus_level,
            en_dropped=en_dropped,
            narrowed_eggnog_groups=narrowed_eggnog_groups,
            unresolved_eggnog_groups=unresolved_eggnog_groups,
        )
        adjudicated_records.extend(gene_adjudicated)

    # Build long-form evidence DataFrame sorted by gene_id, ko, and method
    if adjudicated_records:
        evidence_df = pd.DataFrame(adjudicated_records)[EVIDENCE_COLUMNS]
        evidence_df = evidence_df.sort_values(by=["gene_id", "ko", "method"]).reset_index(drop=True)
    else:
        evidence_df = pd.DataFrame(columns=EVIDENCE_COLUMNS)

    # Define final ordered gene-level schema
    final_cols = [
        "gene_id",
        "accepted_ko",
        "alternative_kos",
        "definition",
        "alternative_definition",
        "consensus_level",
        "evidence",
        "kofam_ko",
        "kofam_score_type",
        "kofam_bit_score",
        "kofam_evalue",
        "kofam_assignment",
        "kofam_threshold",
        "deepkoala_ko",
        "deepkoala_candidate_ko",
        "deepkoala_score",
        "deepkoala_threshold",
        "eggnog_ko",
        "eggnog_candidate_ko",
        "eggnog_bit_score",
        "eggnog_evalue",
    ]

    # Build gene summary DataFrame matching exact schema
    if gene_summary_rows:
        merged_summary = pd.DataFrame(gene_summary_rows)
        output_cols = [c for c in final_cols if c in merged_summary.columns]
        result_df = merged_summary[output_cols].copy()
    else:
        result_df = pd.DataFrame(columns=final_cols)

    # Attach long-form evidence table to attrs for programmatic access
    result_df.attrs["evidence"] = evidence_df

    # Export tables to disk if requested
    export_tables(
        result_df=result_df,
        evidence_df=evidence_df,
        output_tsv=output_tsv,
        evidence_tsv=evidence_tsv,
    )

    return result_df
