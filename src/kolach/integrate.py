"""Consensus and integration of annotation tables for kolach.

This module integrates results from KOfam, DeepKOALA, and eggNOG into a
unified, high-confidence annotation table using pandas, preserving all
method-specific metrics (bit scores, E-values, probabilities, and thresholds)
and evidence provenance. The gene summary table's KOfam metrics (kofam_bit_score,
kofam_evalue) follow each profile's score_type (full vs domain) as indicated by
kofam_score_type, while the long-form evidence table preserves all original full
and domain metrics.
"""
import argparse
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
    "original_status",
    "bit_score",
    "e_value",
    "domain_bit_score",
    "domain_e_value",
    "score_type",
    "threshold",
    "deepkoala_probability",
    "deepkoala_threshold",
    "shared_seed_hit",
    "cross_method_status",
    "supporting_methods",
]

VALID_ORIGINAL_STATUSES = {
    "threshold_passing",
    "heuristic_rescued",
    "below_threshold",
    "unknown",
}


def parse_kos(val) -> set[str]:
    """Extract valid KEGG Orthology identifiers (Kxxxxx) from a string."""
    if pd.isna(val) or val is None or str(val).strip() in ("", "-", "None", "nan"):
        return set()
    tokens = re.split(r"[,+;\s|]+", str(val).strip())
    kos = set()
    for t in tokens:
        m = KO_REGEX.search(t)
        if m:
            kos.add(m.group(0))
    return kos


def format_kos(kos: set[str]) -> str:
    """Format a set of KOs as a sorted comma-separated string."""
    return ",".join(sorted(kos)) if kos else "-"


def read_fasta_ids(fasta_path: Union[str, Path]) -> list[str]:
    """Read all protein IDs from FASTA in original sequence order."""
    fasta_path = Path(fasta_path).expanduser().resolve()
    gene_ids = []
    is_gz = fasta_path.suffix.lower() in (".gz", ".gzip") or str(fasta_path).endswith((".fasta.gz", ".faa.gz", ".fa.gz"))
    open_fn = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8", errors="replace"))
        if is_gz
        else (lambda p: open(p, "r", encoding="utf-8", errors="replace"))
    )

    with open_fn(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split()
                if parts:
                    gene_ids.append(parts[0])
    return gene_ids


def load_ko_definitions(ko_list_path: Optional[Union[str, Path]]) -> dict[str, str]:
    """Load KO definitions from KOfam's ko_list file if available."""
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
                    if KO_REGEX.fullmatch(ko):
                        definition = parts[-1].strip() if len(parts) > 1 else ""
                        definitions[ko] = definition
    except Exception:
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

    - score_type == 'full': selects bit_score and e_value.
    - score_type == 'domain': selects domain_bit_score and domain_e_value.
    Never substitutes a full-sequence score when a required domain score is missing.
    For unrecognized or missing score types, returns (nan, nan, '-').

    Returns (selected_bit_score, selected_e_value, normalized_score_type).
    """
    st_raw = str(score_type).strip().lower() if pd.notna(score_type) else ""
    if st_raw == "full":
        bs = float(bit_score) if (pd.notna(bit_score) and np.isfinite(bit_score)) else np.nan
        ev = float(e_value) if (pd.notna(e_value) and np.isfinite(e_value)) else np.nan
        return bs, ev, "full"
    elif st_raw == "domain":
        dbs = float(domain_bit_score) if (pd.notna(domain_bit_score) and np.isfinite(domain_bit_score)) else np.nan
        dev = float(domain_e_value) if (pd.notna(domain_e_value) and np.isfinite(domain_e_value)) else np.nan
        return dbs, dev, "domain"
    else:
        return np.nan, np.nan, "-"


def extract_kofam_records(tsv_path: Union[str, Path]) -> list[dict]:
    """Extract per-hit KO records from KOfam annotations TSV.

    Duplicate hits for a (gene_id, ko) pair are deterministically resolved by
    selecting the hit with the highest applicable bit score (and lowest E-value)
    according to each profile's score_type (full vs domain).

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
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return []

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

    # Filter to valid K-numbers
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

    # Deterministic deduplication: highest selected bit_score, then lowest selected e_value, then original index
    df = df.sort_values(
        by=["gene_id", "ko", "_sel_bit_score", "_sel_e_value", "_orig_idx"],
        ascending=[True, True, False, True, True],
        na_position="last",
    )
    df = df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in df.iterrows():
        assign_raw = str(row["assignment"]).strip() if pd.notna(row["assignment"]) else ""
        assign_lower = assign_raw.lower()
        score_type = str(row["score_type"]).strip() if pd.notna(row["score_type"]) else "-"
        thresh = row["threshold"]
        bs = row["bit_score"]
        ev = row["e_value"]
        dbs = row["domain_bit_score"]
        dev = row["domain_e_value"]
        sel_bs = row["_sel_bit_score"]
        sel_ev = row["_sel_e_value"]

        # Evaluate status
        if assign_lower in ("threshold", "*"):
            orig_status = "threshold_passing"
        elif assign_lower == "rescued":
            orig_status = "heuristic_rescued"
        elif assign_lower in ("below_threshold", "below", "none"):
            orig_status = "below_threshold"
        else:
            # Missing or unrecognized assignment -> evaluate using score_type and threshold
            st_norm = str(score_type).strip().lower()
            if (
                st_norm in ("full", "domain")
                and pd.notna(sel_bs)
                and np.isfinite(sel_bs)
                and pd.notna(thresh)
                and np.isfinite(thresh)
            ):
                orig_status = "threshold_passing" if sel_bs >= thresh else "below_threshold"
            else:
                orig_status = "unknown"

        definition = str(row["definition"]).strip() if pd.notna(row["definition"]) else "-"
        if definition in ("", "-"):
            definition = "-"

        records.append({
            "gene_id": row["gene_id"],
            "ko": row["ko"],
            "method": "kofam",
            "original_status": orig_status,
            "bit_score": bs,
            "e_value": ev,
            "domain_bit_score": dbs,
            "domain_e_value": dev,
            "score_type": score_type if score_type else "-",
            "threshold": thresh,
            "deepkoala_probability": np.nan,
            "deepkoala_threshold": np.nan,
            "shared_seed_hit": False,
            "definition": definition,
            "assignment": assign_raw if assign_raw else "-",
        })
    return records


def load_kofam(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load and aggregate KOfam annotations table into a gene-level summary.

    Expected columns in input: gene_id, ko, assignment, score_type, threshold,
                              bit_score, e_value, domain_bit_score, domain_e_value, definition.
    Per-KO statuses and scores are preserved without attaching one KO's status
    or score to other KOs. Summary scores and E-values reflect the profile's score_type.
    """
    records = extract_kofam_records(tsv_path)
    keep_cols = [
        "gene_id", "kofam_ko", "kofam_score_type", "kofam_bit_score",
        "kofam_evalue", "kofam_assignment", "kofam_threshold", "kofam_definition"
    ]
    if not records:
        return pd.DataFrame(columns=keep_cols)

    rec_df = pd.DataFrame(records)
    gene_rows = []
    for gene_id, group in rec_df.groupby("gene_id", sort=False):
        sel_tuples = [
            select_kofam_scores(r["score_type"], r["bit_score"], r["e_value"], r["domain_bit_score"], r["domain_e_value"])
            for _, r in group.iterrows()
        ]
        group = group.copy()
        group["_sel_bs"] = [t[0] for t in sel_tuples]
        group["_sel_ev"] = [t[1] for t in sel_tuples]
        group["_sel_st"] = [t[2] for t in sel_tuples]
        group = group.sort_values(by=["_sel_bs"], ascending=False, na_position="last")
        kos = sorted(set(group["ko"]))
        ko_str = format_kos(set(kos))

        if len(group) == 1:
            row = group.iloc[0]
            gene_rows.append({
                "gene_id": gene_id,
                "kofam_ko": ko_str,
                "kofam_score_type": str(row["_sel_st"]),
                "kofam_assignment": str(row["assignment"]),
                "kofam_bit_score": row["_sel_bs"],
                "kofam_evalue": row["_sel_ev"],
                "kofam_threshold": row["threshold"],
                "kofam_definition": str(row["definition"]),
            })
        else:
            # Multi-KO gene: map each KO unambiguously to its own assignment and scores.
            # Plain score types ordered matching kofam_ko
            kos_order = ko_str.split(",")
            ko_to_st = {r["ko"]: str(r["_sel_st"]) for _, r in group.iterrows()}
            score_type_str = ",".join(ko_to_st[k] for k in kos_order if k in ko_to_st)

            assign_str = ",".join(f"{r['ko']}:{r['assignment']}" for _, r in group.iterrows())
            score_str = ",".join(f"{r['ko']}:{r['_sel_bs']}" for _, r in group.iterrows())
            evalue_str = ",".join(f"{r['ko']}:{r['_sel_ev']}" for _, r in group.iterrows())
            thresh_str = ",".join(f"{r['ko']}:{r['threshold']}" for _, r in group.iterrows())
            defs = [f"{r['ko']}: {r['definition']}" for _, r in group.iterrows() if str(r["definition"]) not in ("", "-")]
            def_str = "; ".join(defs) if defs else "-"

            gene_rows.append({
                "gene_id": gene_id,
                "kofam_ko": ko_str,
                "kofam_score_type": score_type_str,
                "kofam_assignment": assign_str,
                "kofam_bit_score": score_str,
                "kofam_evalue": evalue_str,
                "kofam_threshold": thresh_str,
                "kofam_definition": def_str,
            })

    return pd.DataFrame(gene_rows)[keep_cols]


def extract_deepkoala_records(tsv_path: Union[str, Path]) -> list[dict]:
    """Extract per-hit KO records from DeepKOALA annotations TSV.

    Duplicate hits for a (gene_id, ko) pair are deterministically resolved by
    selecting the hit with the highest probability.

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
    - Arbitrary two-column tables without documented headers: unknown.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return []

    cleaned_raw_cols = [c.lstrip("#").strip() for c in df.columns]
    cols_lower = [c.lower() for c in cleaned_raw_cols]

    # Check if this matches documented DeepKOALA simple-output format:
    is_documented_simple = (
        len(df.columns) == 2
        and "name" in cols_lower
        and "predict_label" in cols_lower
        and "annotate" not in cols_lower
        and "probability" not in cols_lower
        and "score" not in cols_lower
    )

    col_map = {}
    if "name" in cols_lower:
        col_map[df.columns[cols_lower.index("name")]] = "gene_id"
    elif "gene_id" in cols_lower:
        col_map[df.columns[cols_lower.index("gene_id")]] = "gene_id"
    else:
        col_map[df.columns[0]] = "gene_id"

    if "predict_label" in cols_lower:
        col_map[df.columns[cols_lower.index("predict_label")]] = "predict_label"
    elif "ko" in cols_lower:
        col_map[df.columns[cols_lower.index("ko")]] = "predict_label"
    elif len(df.columns) > 1:
        col_map[df.columns[1]] = "predict_label"

    if "probability" in cols_lower:
        col_map[df.columns[cols_lower.index("probability")]] = "deepkoala_score"
    elif "score" in cols_lower:
        col_map[df.columns[cols_lower.index("score")]] = "deepkoala_score"

    if "threshold" in cols_lower:
        col_map[df.columns[cols_lower.index("threshold")]] = "deepkoala_threshold"

    if "annotate" in cols_lower:
        col_map[df.columns[cols_lower.index("annotate")]] = "deepkoala_annotate"

    df = df.rename(columns=col_map)
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["deepkoala_score"] = pd.to_numeric(df.get("deepkoala_score"), errors="coerce")
    df["deepkoala_threshold"] = pd.to_numeric(df.get("deepkoala_threshold"), errors="coerce")

    exploded_rows = []
    for _, row in df.iterrows():
        raw_val = row.get("predict_label", "-")
        kos = parse_kos(raw_val)
        score = row.get("deepkoala_score", np.nan)
        thresh = row.get("deepkoala_threshold", np.nan)
        ann = str(row.get("deepkoala_annotate", "")).strip() if "deepkoala_annotate" in df.columns and pd.notna(row.get("deepkoala_annotate")) else ""

        # Determine original status per row
        if is_documented_simple:
            orig_status = "threshold_passing"
            score_type = "upstream_filtered"
        elif "deepkoala_annotate" in df.columns:
            if "*" in ann:
                orig_status = "threshold_passing"
                score_type = "probability"
            elif pd.notna(score) and pd.notna(thresh):
                orig_status = "threshold_passing" if score >= thresh else "below_threshold"
                score_type = "probability"
            else:
                orig_status = "unknown"
                score_type = "-"
        elif pd.notna(score) and pd.notna(thresh):
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
    # Deduplicate by highest score (unknown or NaN scores sort to end)
    exp_df = exp_df.sort_values(by=["gene_id", "ko", "deepkoala_score"], ascending=[True, True, False])
    exp_df = exp_df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in exp_df.iterrows():
        records.append({
            "gene_id": row["gene_id"],
            "ko": row["ko"],
            "method": "deepkoala",
            "original_status": row["original_status"],
            "bit_score": np.nan,
            "e_value": np.nan,
            "domain_bit_score": np.nan,
            "domain_e_value": np.nan,
            "score_type": row["score_type"],
            "threshold": np.nan,
            "deepkoala_probability": row["deepkoala_score"],
            "deepkoala_threshold": row["deepkoala_threshold"],
            "shared_seed_hit": False,
            "definition": "-",
        })
    return records


def load_deepkoala(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load and format DeepKOALA annotations table."""
    keep_cols = ["gene_id", "deepkoala_ko", "deepkoala_candidate_ko", "deepkoala_score", "deepkoala_threshold"]
    records = extract_deepkoala_records(tsv_path)
    if not records:
        return pd.DataFrame(columns=keep_cols)

    rec_df = pd.DataFrame(records)
    gene_rows = []
    for gene_id, group in rec_df.groupby("gene_id", sort=False):
        # Candidates are valid below_threshold or threshold_passing hits
        cands = sorted(set(group.loc[group["original_status"].isin(["threshold_passing", "below_threshold"]), "ko"]))
        thresh_kos = sorted(set(group.loc[group["original_status"] == "threshold_passing", "ko"]))

        best_row = group.sort_values(by="deepkoala_probability", ascending=False).iloc[0]
        gene_rows.append({
            "gene_id": gene_id,
            "deepkoala_ko": format_kos(set(thresh_kos)),
            "deepkoala_candidate_ko": format_kos(set(cands)),
            "deepkoala_score": best_row["deepkoala_probability"],
            "deepkoala_threshold": best_row["deepkoala_threshold"],
        })
    return pd.DataFrame(gene_rows)[keep_cols]


def read_eggnog_tsv(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Read eggNOG annotations table, properly handling native and normalized headers.

    Supports:
    - Native .emapper.annotations files containing '##' metadata and a '#query' header.
    - Kolach-normalized TSVs containing an unprefixed 'query' header.
    - Preserves column aliases ('query'/'gene_id', 'kegg_ko'/'ko', 'score', 'evalue', 'description').
    - Raises ValueError when required identifying or KO columns are missing.
    - Valid header-only files return an empty DataFrame with the parsed columns.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"eggNOG file not found: '{path}'")
    if path.stat().st_size == 0:
        raise ValueError(f"eggNOG annotations file is empty (0 bytes): '{path}'")

    is_gz = path.suffix.lower() in (".gz", ".gzip") or str(path).endswith(".annotations.gz")
    open_fn = (
        (lambda p: gzip.open(p, "rt", encoding="utf-8", errors="replace"))
        if is_gz
        else (lambda p: open(p, "r", encoding="utf-8", errors="replace"))
    )

    header_line = None
    data_lines = []
    with open_fn(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("##"):
                continue
            if header_line is None:
                # Skip comment lines before the header except the recognized #query header
                if stripped.startswith("#"):
                    first_cell = re.sub(r"\s+", "", stripped.split("\t")[0].lower())
                    if first_cell not in ("#query", "#gene_id"):
                        continue
                header_line = line.rstrip("\r\n")
            else:
                if stripped.startswith("#"):
                    continue
                data_lines.append(line)

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

    if not data_lines:
        return pd.DataFrame(columns=clean_cols)

    import io
    content = "\t".join(clean_cols) + "\n" + "".join(data_lines)
    return pd.read_csv(io.StringIO(content), sep="\t", dtype=str)


def extract_eggnog_records(
    tsv_path: Union[str, Path],
    min_bitscore: float = 60.0,
    max_evalue: float = 1e-5,
) -> list[dict]:
    """Extract per-hit KO records from eggNOG annotations TSV.

    Shared seed-hit metrics for multi-KO assignments are explicitly tagged.
    Duplicate hits for a (gene_id, ko) pair are deterministically resolved by
    selecting the hit with the highest bit_score (and lowest e_value).

    Status determination:
    - Missing, malformed, or nonfinite bit score or E-value -> unknown
    - Valid finite metrics -> threshold_passing if bit_score >= min_bitscore and evalue <= max_evalue,
      else below_threshold.
    Preserves legitimate E-values of zero.
    eggNOG descriptions are generic seed-hit descriptions and are preserved as
    metadata, NOT authoritative KO-specific definitions.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = read_eggnog_tsv(path)
    if df.empty:
        return []

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
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["eggnog_bit_score"] = pd.to_numeric(df.get("eggnog_bit_score"), errors="coerce")
    df["eggnog_evalue"] = pd.to_numeric(df.get("eggnog_evalue"), errors="coerce")
    if "eggnog_description" not in df.columns:
        df["eggnog_description"] = "-"
    if "eggnog_raw_ko" not in df.columns:
        df["eggnog_raw_ko"] = "-"

    exploded_rows = []
    for _, row in df.iterrows():
        raw_val = row["eggnog_raw_ko"]
        kos = parse_kos(raw_val)
        if not kos:
            continue
        is_multi = len(kos) > 1
        bs = row["eggnog_bit_score"]
        ev = row["eggnog_evalue"]

        # Status evaluation: require finite bit score and E-value
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
                "shared_seed_hit": is_multi,
                "original_status": orig_status,
                "eggnog_description": str(row["eggnog_description"]).strip() if pd.notna(row["eggnog_description"]) else "-",
            })

    if not exploded_rows:
        return []

    exp_df = pd.DataFrame(exploded_rows)
    # Deduplicate: prefer records with valid finite metrics, then highest bit_score, lowest e_value
    exp_df["_valid_metrics"] = exp_df["original_status"].isin(["threshold_passing", "below_threshold"])
    exp_df["_orig_idx"] = np.arange(len(exp_df))
    exp_df = exp_df.sort_values(
        by=["gene_id", "ko", "_valid_metrics", "bit_score", "e_value", "_orig_idx"],
        ascending=[True, True, False, False, True, True],
        na_position="last",
    )
    exp_df = exp_df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in exp_df.iterrows():
        records.append({
            "gene_id": row["gene_id"],
            "ko": row["ko"],
            "method": "eggnog",
            "original_status": row["original_status"],
            "bit_score": row["bit_score"],
            "e_value": row["e_value"],
            "domain_bit_score": np.nan,
            "domain_e_value": np.nan,
            "score_type": "seed_hit",
            "threshold": min_bitscore,
            "deepkoala_probability": np.nan,
            "deepkoala_threshold": np.nan,
            "shared_seed_hit": bool(row["shared_seed_hit"]),
            "definition": "-",  # generic seed description is NOT authoritative KO definition
            "eggnog_description": row["eggnog_description"],
        })
    return records


def load_eggnog(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load raw eggNOG-mapper annotations table."""
    path = Path(tsv_path).expanduser().resolve()
    keep_cols = ["gene_id", "eggnog_raw_ko", "eggnog_candidate_ko", "eggnog_bit_score", "eggnog_evalue", "eggnog_description"]
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=keep_cols)

    df = read_eggnog_tsv(path)
    if df.empty:
        return pd.DataFrame(columns=keep_cols)

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
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    if "eggnog_raw_ko" not in df.columns:
        df["eggnog_raw_ko"] = "-"
    df["eggnog_bit_score"] = pd.to_numeric(df.get("eggnog_bit_score"), errors="coerce")
    df["eggnog_evalue"] = pd.to_numeric(df.get("eggnog_evalue"), errors="coerce")
    if "eggnog_description" not in df.columns:
        df["eggnog_description"] = "-"

    # Only valid finite metrics can produce candidates
    df["eggnog_candidate_ko"] = df.apply(
        lambda r: format_kos(parse_kos(r["eggnog_raw_ko"]))
        if (
            pd.notna(r["eggnog_bit_score"])
            and pd.notna(r["eggnog_evalue"])
            and np.isfinite(r["eggnog_bit_score"])
            and np.isfinite(r["eggnog_evalue"])
        )
        else "-",
        axis=1,
    )
    df["_valid_metrics"] = (
        pd.notna(df["eggnog_bit_score"])
        & pd.notna(df["eggnog_evalue"])
        & np.isfinite(df["eggnog_bit_score"])
        & np.isfinite(df["eggnog_evalue"])
    )
    df["_orig_idx"] = np.arange(len(df))
    df = df.sort_values(
        by=["gene_id", "_valid_metrics", "eggnog_bit_score", "eggnog_evalue", "_orig_idx"],
        ascending=[True, False, False, True, True],
        na_position="last",
    )
    return df[keep_cols].drop_duplicates(subset=["gene_id"])


def extract_kofam_status_sets(kofam_ko: str, kofam_assignment: str) -> tuple[set[str], set[str]]:
    """Separate confident threshold-passing KOs from heuristic-rescued KOs.

    Missing or unrecognized assignments do NOT default to threshold_passing.
    """
    kos = parse_kos(kofam_ko)
    if not kos:
        return set(), set()
    assign_str = str(kofam_assignment).strip()
    if not assign_str or assign_str == "-":
        return set(), set()

    if ":" in assign_str:
        thresh_kos = set()
        rescued_kos = set()
        for part in assign_str.split(","):
            if ":" in part:
                k, a = part.split(":", 1)
                k = k.strip()
                a = a.strip().lower()
                if a in ("threshold", "*"):
                    thresh_kos.add(k)
                elif a == "rescued":
                    rescued_kos.add(k)
        return (thresh_kos & kos), (rescued_kos & kos)

    if assign_str.lower() in ("threshold", "*"):
        return kos, set()
    if assign_str.lower() == "rescued":
        return set(), kos
    return set(), set()


def get_single_ko_definition(
    ko: str,
    ko_definitions: dict[str, str],
    evidence_definitions: Optional[dict[str, str]] = None,
) -> str:
    """Resolve functional definition for a single KO.

    Strict priority:
    1. Master ko_list definition.
    2. Matching KO-specific definition from evidence records (e.g. KOfam profile definition).
    3. Return '-' if no matching definition is available.
    """
    if ko in ko_definitions and ko_definitions[ko] and ko_definitions[ko] != "-":
        return str(ko_definitions[ko]).strip()
    if evidence_definitions and ko in evidence_definitions:
        edef = evidence_definitions[ko]
        if edef and edef != "-":
            return str(edef).strip()
    return "-"


def resolve_definition_for_kos(
    kos: set[str],
    ko_definitions: dict[str, str],
    evidence_definitions: Optional[dict[str, str]] = None,
) -> str:
    """Resolve functional definition text for a set of KOs.

    Preserves unambiguous KO-definition associations:
    - Single KO: returns definition string (or '-')
    - Multiple KOs: returns 'Kxxxxx: def1; Kyyyyy: def2' for KOs with definitions,
      or '-' if none found.
    """
    if not kos:
        return "-"
    sorted_kos = sorted(kos)
    if len(sorted_kos) == 1:
        return get_single_ko_definition(sorted_kos[0], ko_definitions, evidence_definitions)

    parts = []
    for k in sorted_kos:
        d = get_single_ko_definition(k, ko_definitions, evidence_definitions)
        if d != "-":
            parts.append(f"{k}: {d}")
    return "; ".join(parts) if parts else "-"


def _resolve_definition(
    kos: set[str],
    row: Union[pd.Series, dict],
    ko_definitions: dict[str, str],
) -> str:
    """Find the best definition text strictly for the specified KO(s).

    Never leaks an unrelated KO's definition (e.g. losing KOfam KO) to another KO.
    """
    if not kos:
        return "-"

    # Build KO-specific evidence definitions from row if available
    evidence_defs: dict[str, str] = {}
    if isinstance(row, pd.Series) or isinstance(row, dict):
        kf_ko = parse_kos(row.get("kofam_ko", "-"))
        kf_def = str(row.get("kofam_definition", "-")).strip()
        if kf_def and kf_def != "-" and len(kf_ko) == 1:
            evidence_defs[next(iter(kf_ko))] = kf_def
        elif kf_def and kf_def != "-" and ":" in kf_def:
            # Format: 'K00001: def1; K00002: def2'
            for part in kf_def.split(";"):
                if ":" in part:
                    k, d = part.split(":", 1)
                    evidence_defs[k.strip()] = d.strip()

    return resolve_definition_for_kos(kos, ko_definitions, evidence_defs)


def _assign_accepted_or_alt(
    candidate_kos: set[str],
    row: pd.Series,
    ko_definitions: dict[str, str],
    additional_alternatives: Optional[set[str]] = None,
) -> tuple[str, str, str, str]:
    """Return (accepted_ko, alternative_kos, definition, alternative_definition).

    Enforces a strict single-hit policy for accepted_ko: accepted_ko contains
    single hits only (len == 1). Multiple hits (len > 1) are placed into
    alternative_kos to avoid downstream pathway reconstruction tools misinterpreting
    comma-separated multi-hits as multifunctional enzymes.
    """
    alts = set(additional_alternatives) if additional_alternatives else set()

    if len(candidate_kos) == 1:
        accepted = format_kos(candidate_kos)
        alt_set = alts - candidate_kos
        alt = format_kos(alt_set) if alt_set else "-"
        defn = _resolve_definition(candidate_kos, row, ko_definitions)
        alt_defn = _resolve_definition(alt_set, row, ko_definitions) if alt_set else "-"
    elif len(candidate_kos) > 1:
        accepted = "-"
        alt_set = candidate_kos | alts
        alt = format_kos(alt_set)
        defn = "-"
        alt_defn = _resolve_definition(alt_set, row, ko_definitions)
    else:
        accepted = "-"
        alt_set = alts
        alt = format_kos(alt_set) if alt_set else "-"
        defn = "-"
        alt_defn = _resolve_definition(alt_set, row, ko_definitions) if alt_set else "-"

    return accepted, alt, defn, alt_defn


def filter_and_disambiguate_eggnog(
    eggnog_df: pd.DataFrame,
    merged_df: pd.DataFrame,
    min_bitscore: float = 60.0,
    max_evalue: float = 1e-5,
    filter_multi_mode: str = "disambiguate",
) -> pd.DataFrame:
    """Filter eggNOG hits by bitscore and disambiguate multi-KO assignments.

    Below-threshold eggNOG hits are NEVER promoted to confident calls in eggnog_ko.
    """
    if eggnog_df.empty:
        merged_df["eggnog_ko"] = "-"
        merged_df["eggnog_candidate_ko"] = "-"
        merged_df["eggnog_bit_score"] = np.nan
        merged_df["eggnog_evalue"] = np.nan
        merged_df["eggnog_description"] = "-"
        return merged_df

    merged_df = merged_df.merge(eggnog_df, on="gene_id", how="left")
    merged_df["eggnog_raw_ko"] = merged_df.get("eggnog_raw_ko", pd.Series("-", index=merged_df.index)).fillna("-")
    merged_df["eggnog_candidate_ko"] = merged_df.get("eggnog_candidate_ko", pd.Series("-", index=merged_df.index)).fillna("-")

    filtered_kos = []
    for _, row in merged_df.iterrows():
        raw_val = row.get("eggnog_raw_ko", "-")
        bitscore = row.get("eggnog_bit_score", np.nan)
        evalue = row.get("eggnog_evalue", np.nan)
        en_kos = parse_kos(raw_val)

        kf_thresh, kf_rescued = extract_kofam_status_sets(row.get("kofam_ko", "-"), row.get("kofam_assignment", "-"))
        dk_kos = parse_kos(row.get("deepkoala_ko", "-"))
        dk_cand = parse_kos(row.get("deepkoala_candidate_ko", "-")) - dk_kos
        trusted_other = kf_thresh | dk_kos
        cand_other = dk_cand | kf_rescued

        # 1. Below threshold, missing, or nonfinite metrics: never promote to confident eggnog_ko!
        has_valid_metrics = (
            pd.notna(bitscore) and pd.notna(evalue)
            and np.isfinite(bitscore) and np.isfinite(evalue)
            and bitscore >= min_bitscore and evalue <= max_evalue
        )
        if not en_kos or not has_valid_metrics:
            filtered_kos.append("-")
            continue

        # 2. Single-KO meeting threshold
        if len(en_kos) == 1 or filter_multi_mode == "none":
            filtered_kos.append(format_kos(en_kos))
            continue

        # 3. Multi-KO meeting threshold
        overlap = trusted_other & en_kos
        if overlap:
            filtered_kos.append(format_kos(overlap))
        elif cand_other & en_kos:
            filtered_kos.append(format_kos(cand_other & en_kos))
        elif trusted_other:
            if filter_multi_mode == "strict":
                filtered_kos.append("-")
            else:
                filtered_kos.append(format_kos(en_kos))
        else:
            if filter_multi_mode == "strict":
                filtered_kos.append("-")
            else:
                filtered_kos.append(format_kos(en_kos))

    merged_df["eggnog_ko"] = filtered_kos
    merged_df.drop(columns=["eggnog_raw_ko"], inplace=True, errors="ignore")
    return merged_df


def adjudicate_consensus(
    merged_df: pd.DataFrame,
    active_tools: list[str],
    conflict_strategy: str = "multiple",
    ko_definitions: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    """Evaluate consensus KO, consensus level, evidence, and functional definition.

    - Computes confident agreement using ONLY calls that independently passed
      thresholds (KOfam heuristic rescues, sub-threshold candidates, and unknown
      evidence do NOT count as confident votes).
    - Preserves unselected minority threshold-passing calls from all methods in
      alternative_kos without breaking majority consensus.
    - Resolves definitions strictly by matching KO.
    """
    if ko_definitions is None:
        ko_definitions = {}

    accepted_kos = []
    alternative_kos_list = []
    consensus_levels = []
    evidence_list = []
    definitions = []
    alt_definitions = []

    priority_order = [t for t in ["kofam", "eggnog", "deepkoala"] if t in active_tools]

    for _, row in merged_df.iterrows():
        # Identify confident threshold-passing calls per tool
        tool_calls: dict[str, set[str]] = {}

        kf_thresh, kf_rescued = extract_kofam_status_sets(row.get("kofam_ko", "-"), row.get("kofam_assignment", "-"))
        if "kofam" in active_tools and kf_thresh:
            tool_calls["kofam"] = kf_thresh

        dk_thresh = parse_kos(row.get("deepkoala_ko", "-"))
        if "deepkoala" in active_tools and dk_thresh:
            tool_calls["deepkoala"] = dk_thresh

        en_thresh = parse_kos(row.get("eggnog_ko", "-"))
        if "eggnog" in active_tools and en_thresh:
            tool_calls["eggnog"] = en_thresh

        # Sub-threshold candidates (strictly excluding unknown)
        dk_cand = parse_kos(row.get("deepkoala_candidate_ko", "-")) - dk_thresh
        en_cand = parse_kos(row.get("eggnog_candidate_ko", "-")) - en_thresh
        en_raw = parse_kos(row.get("eggnog_candidate_ko", "-"))

        num_calling_tools = len(tool_calls)

        # Case 0: No tool made an independently confident threshold-passing call
        if num_calling_tools == 0:
            agreeing_cands = []
            cand_agree = set()
            if "deepkoala" in active_tools and "eggnog" in active_tools:
                cand_agree = dk_cand & en_cand
                if cand_agree:
                    agreeing_cands = ["deepkoala(candidate)", "eggnog(candidate)"]
            if not cand_agree and "kofam" in active_tools and "deepkoala" in active_tools:
                cand_agree = kf_rescued & dk_cand
                if cand_agree:
                    agreeing_cands = ["deepkoala(candidate)", "kofam(rescued)"]
            if not cand_agree and "kofam" in active_tools and "eggnog" in active_tools:
                cand_agree = kf_rescued & en_cand
                if cand_agree:
                    agreeing_cands = ["eggnog(candidate)", "kofam(rescued)"]

            if cand_agree:
                dropped_en = (en_raw - cand_agree) if (len(en_raw) > 1 and bool(cand_agree & en_raw)) else set()
                acc, alt, defn, alt_defn = _assign_accepted_or_alt(cand_agree, row, ko_definitions, additional_alternatives=dropped_en)
                accepted_kos.append(acc)
                alternative_kos_list.append(alt)
                consensus_levels.append("orthogonal_dual_candidate")
                evidence_list.append(",".join(sorted(agreeing_cands)))
                definitions.append(defn)
                alt_definitions.append(alt_defn)
            elif kf_rescued:
                dropped_en = (en_raw - kf_rescued) if (len(en_raw) > 1 and bool(kf_rescued & en_raw)) else set()
                acc, alt, defn, alt_defn = _assign_accepted_or_alt(kf_rescued, row, ko_definitions, additional_alternatives=dropped_en)
                accepted_kos.append(acc)
                alternative_kos_list.append(alt)
                consensus_levels.append("single_tool")
                evidence_list.append("kofam(rescued)")
                definitions.append(defn)
                alt_definitions.append(alt_defn)
            else:
                accepted_kos.append("-")
                alternative_kos_list.append("-")
                consensus_levels.append("unannotated")
                evidence_list.append("-")
                definitions.append("-")
                alt_definitions.append("-")
            continue

        # Case 1: Exactly 1 confident tool call
        if num_calling_tools == 1:
            tool_name = next(iter(tool_calls))
            kos = tool_calls[tool_name]

            rescued_by = []
            if tool_name != "deepkoala" and "deepkoala" in active_tools and (kos & dk_cand):
                rescued_by.append("deepkoala(candidate)")
            if tool_name != "eggnog" and "eggnog" in active_tools and (kos & en_cand):
                rescued_by.append("eggnog(candidate)")
            if tool_name != "kofam" and "kofam" in active_tools and (kos & kf_rescued):
                rescued_by.append("kofam(rescued)")

            dropped_en = (en_raw - kos) if (len(en_raw) > 1 and bool(kos & en_raw)) else set()
            acc, alt, defn, alt_defn = _assign_accepted_or_alt(kos, row, ko_definitions, additional_alternatives=dropped_en)
            accepted_kos.append(acc)
            alternative_kos_list.append(alt)
            definitions.append(defn)
            alt_definitions.append(alt_defn)

            if rescued_by:
                consensus_level = "single_tool_with_candidate"
                evidence_list.append(",".join([tool_name] + rescued_by))
            else:
                consensus_level = "single_tool"
                evidence_list.append(tool_name)
            consensus_levels.append(consensus_level)
            continue

        # Case 2: 2 or more confident tools made calls
        all_sets = list(tool_calls.values())
        union_all_calls = set.union(*all_sets)
        common_all = set.intersection(*all_sets)

        if common_all:
            is_unanimous = (num_calling_tools == len(active_tools))
            dropped_en = (en_raw - common_all) if (len(en_raw) > 1 and bool(common_all & en_raw)) else set()
            unselected_minority = union_all_calls - common_all
            all_extra_alts = unselected_minority | dropped_en

            acc, alt, defn, alt_defn = _assign_accepted_or_alt(common_all, row, ko_definitions, additional_alternatives=all_extra_alts)
            accepted_kos.append(acc)
            alternative_kos_list.append(alt)
            consensus_levels.append("unanimous" if is_unanimous else "majority")
            evidence_list.append(",".join(sorted(tool_calls.keys())))
            definitions.append(defn)
            alt_definitions.append(alt_defn)
            continue

        # Pairwise check for majority (e.g. 2 out of 3 confident tools agree)
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
            dropped_en = (en_raw - pairwise_agreed) if (len(en_raw) > 1 and bool(pairwise_agreed & en_raw)) else set()
            unselected_minority = union_all_calls - pairwise_agreed
            all_extra_alts = unselected_minority | dropped_en

            acc, alt, defn, alt_defn = _assign_accepted_or_alt(pairwise_agreed, row, ko_definitions, additional_alternatives=all_extra_alts)
            accepted_kos.append(acc)
            alternative_kos_list.append(alt)
            consensus_levels.append("majority")
            evidence_list.append(",".join(sorted(agreeing_tools)))
            definitions.append(defn)
            alt_definitions.append(alt_defn)
            continue

        # Disjoint conflict: zero overlap between confident tools
        union_kos = union_all_calls
        dropped_en = (en_raw - union_kos) if (len(en_raw) > 1 and bool(union_kos & en_raw)) else set()
        all_alts = union_kos | dropped_en

        if conflict_strategy == "priority":
            top_tool = next((t for t in priority_order if t in tool_calls), tools_list[0])
            top_kos = tool_calls[top_tool]
            unselected = union_kos - top_kos
            dropped_en_top = (en_raw - top_kos) if (len(en_raw) > 1 and bool(top_kos & en_raw)) else set()
            extra_alts = unselected | dropped_en_top
            if len(top_kos) == 1:
                accepted_kos.append(format_kos(top_kos))
                alternative_kos_list.append(format_kos(extra_alts) if extra_alts else "-")
                definitions.append(_resolve_definition(top_kos, row, ko_definitions))
                alt_definitions.append(_resolve_definition(extra_alts, row, ko_definitions) if extra_alts else "-")
            else:
                accepted_kos.append("-")
                alternative_kos_list.append(format_kos(all_alts))
                definitions.append("-")
                alt_definitions.append(_resolve_definition(all_alts, row, ko_definitions))
            consensus_levels.append("conflict_priority")
            evidence_list.append(top_tool)
        else:  # "multiple" / "union" (default)
            accepted_kos.append("-")
            alternative_kos_list.append(format_kos(all_alts))
            consensus_levels.append("conflict")
            evidence_list.append(",".join(sorted(tool_calls.keys())))
            definitions.append("-")
            alt_definitions.append(_resolve_definition(all_alts, row, ko_definitions))

    merged_df["accepted_ko"] = accepted_kos
    merged_df["ko"] = accepted_kos  # Compatibility alias
    merged_df["alternative_kos"] = alternative_kos_list
    merged_df["definition"] = definitions
    merged_df["alternative_definition"] = alt_definitions
    merged_df["consensus_level"] = consensus_levels
    merged_df["evidence"] = evidence_list

    return merged_df


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
) -> pd.DataFrame:
    """Main evidence-first integration pipeline for kolach annotation outputs.

    Parses each method once into normalized per-gene, per-KO records.
    Adjudicates filtering, disambiguation, rescue, and consensus directly from
    those records, guaranteeing consistency between the gene summary table
    and the long-form evidence table.
    """
    active_tools = []
    all_evidence_records: list[dict] = []
    evidence_defs: dict[str, str] = {}  # KO-specific definitions from KOfam

    # Validate explicit input paths
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
                    f"Explicitly supplied {arg_name} file does not exist or is not a regular file: '{p}' (resolved: '{norm_p}')"
                )

    # 1. Ingest KOfam
    if kofam_tsv is not None:
        kf_records = extract_kofam_records(kofam_tsv)
        all_evidence_records.extend(kf_records)
        for r in kf_records:
            if r.get("definition") and r["definition"] != "-":
                evidence_defs[r["ko"]] = r["definition"]
        active_tools.append("kofam")

    # 2. Ingest DeepKOALA
    if deepkoala_tsv is not None:
        dk_records = extract_deepkoala_records(deepkoala_tsv)
        all_evidence_records.extend(dk_records)
        active_tools.append("deepkoala")

    # 3. Ingest eggNOG
    if eggnog_tsv is not None:
        en_records = extract_eggnog_records(
            eggnog_tsv,
            min_bitscore=eggnog_min_bitscore,
            max_evalue=eggnog_max_evalue,
        )
        all_evidence_records.extend(en_records)
        active_tools.append("eggnog")

    # Determine universe and order of gene IDs
    if protein_fasta is not None:
        all_gene_ids = read_fasta_ids(protein_fasta)
    else:
        seen = set()
        ordered_ids = []
        for r in all_evidence_records:
            gid = r["gene_id"]
            if gid not in seen:
                seen.add(gid)
                ordered_ids.append(gid)
        all_gene_ids = ordered_ids

    # Load master KO definitions from ko_list if available
    ko_list_file = None
    if database_dir:
        db_path = Path(database_dir).expanduser().resolve()
        if (db_path / "kofam" / "ko_list").is_file():
            ko_list_file = db_path / "kofam" / "ko_list"
        elif (db_path / "ko_list").is_file():
            ko_list_file = db_path / "ko_list"
    master_ko_defs = load_ko_definitions(ko_list_file)

    # Group evidence records by gene_id
    records_by_gene: dict[str, list[dict]] = {gid: [] for gid in all_gene_ids}
    for r in all_evidence_records:
        gid = r["gene_id"]
        if gid not in records_by_gene:
            records_by_gene[gid] = []
            all_gene_ids.append(gid)
        records_by_gene[gid].append(r)

    # Adjudicate consensus and cross-method evidence per gene directly from records
    gene_summary_rows = []
    adjudicated_records = []

    priority_order = [t for t in ["kofam", "eggnog", "deepkoala"] if t in active_tools]

    for gid in all_gene_ids:
        recs = records_by_gene.get(gid, [])

        # Categorize records by tool and status
        tool_recs: dict[str, list[dict]] = {t: [] for t in active_tools}
        for r in recs:
            t = r["method"]
            if t in tool_recs:
                tool_recs[t].append(r)

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

        # eggNOG multi-KO disambiguation
        en_agreed = set(confident_kos.get("eggnog", set()))
        en_dropped = set()
        en_raw_cands = {r["ko"] for r in tool_recs.get("eggnog", []) if r["original_status"] in ("threshold_passing", "below_threshold")}

        if "eggnog" in active_tools and en_agreed:
            # Check if any eggNOG record is multi-KO
            has_multi = any(r.get("shared_seed_hit", False) for r in tool_recs.get("eggnog", []))
            if has_multi:
                trusted_other = set().union(*[confident_kos.get(t, set()) for t in active_tools if t != "eggnog"])
                cand_other = set().union(*[rescued_kos.get(t, set()) | candidate_kos.get(t, set()) for t in active_tools if t != "eggnog"])

                if eggnog_filter_multi == "disambiguate":
                    if en_agreed & trusted_other:
                        overlap = en_agreed & trusted_other
                        en_dropped = en_agreed - overlap
                        en_agreed = overlap
                    elif en_agreed & cand_other:
                        overlap = en_agreed & cand_other
                        en_dropped = en_agreed - overlap
                        en_agreed = overlap
                elif eggnog_filter_multi == "strict":
                    if en_agreed & trusted_other:
                        overlap = en_agreed & trusted_other
                        en_dropped = en_agreed - overlap
                        en_agreed = overlap
                    else:
                        en_dropped = en_agreed
                        en_agreed = set()

                if en_agreed:
                    confident_kos["eggnog"] = en_agreed
                else:
                    confident_kos.pop("eggnog", None)

        # Adjudicate consensus directly
        tool_calls = {t: confident_kos[t] for t in active_tools if t in confident_kos}
        num_calling_tools = len(tool_calls)

        accepted_ko = "-"
        alt_set = set()
        consensus_level = "unannotated"
        evidence_str = "-"

        if num_calling_tools == 0:
            dk_cand = candidate_kos.get("deepkoala", set())
            en_cand = candidate_kos.get("eggnog", set())
            kf_resc = rescued_kos.get("kofam", set())

            cand_agree = set()
            agreeing_cands = []
            if "deepkoala" in active_tools and "eggnog" in active_tools and (dk_cand & en_cand):
                cand_agree = dk_cand & en_cand
                agreeing_cands = ["deepkoala(candidate)", "eggnog(candidate)"]
            if not cand_agree and "kofam" in active_tools and "deepkoala" in active_tools and (kf_resc & dk_cand):
                cand_agree = kf_resc & dk_cand
                agreeing_cands = ["deepkoala(candidate)", "kofam(rescued)"]
            if not cand_agree and "kofam" in active_tools and "eggnog" in active_tools and (kf_resc & en_cand):
                cand_agree = kf_resc & en_cand
                agreeing_cands = ["eggnog(candidate)", "kofam(rescued)"]

            if cand_agree:
                consensus_level = "orthogonal_dual_candidate"
                evidence_str = ",".join(sorted(agreeing_cands))
                dropped_en = (en_raw_cands - cand_agree) if (len(en_raw_cands) > 1 and bool(cand_agree & en_raw_cands)) else en_dropped
                if len(cand_agree) == 1:
                    accepted_ko = next(iter(cand_agree))
                    alt_set = dropped_en - {accepted_ko}
                else:
                    accepted_ko = "-"
                    alt_set = cand_agree | dropped_en
            elif kf_resc:
                consensus_level = "single_tool"
                evidence_str = "kofam(rescued)"
                dropped_en = (en_raw_cands - kf_resc) if (len(en_raw_cands) > 1 and bool(kf_resc & en_raw_cands)) else en_dropped
                if len(kf_resc) == 1:
                    accepted_ko = next(iter(kf_resc))
                    alt_set = dropped_en - {accepted_ko}
                else:
                    accepted_ko = "-"
                    alt_set = kf_resc | dropped_en
            else:
                accepted_ko = "-"
                alt_set = en_dropped
                consensus_level = "unannotated"
                evidence_str = "-"

        elif num_calling_tools == 1:
            tool_name = next(iter(tool_calls))
            called_kos = tool_calls[tool_name]

            dk_cand = candidate_kos.get("deepkoala", set())
            en_cand = candidate_kos.get("eggnog", set())
            kf_resc = rescued_kos.get("kofam", set())

            rescued_by = []
            if tool_name != "deepkoala" and "deepkoala" in active_tools and (called_kos & dk_cand):
                rescued_by.append("deepkoala(candidate)")
            if tool_name != "eggnog" and "eggnog" in active_tools and (called_kos & en_cand):
                rescued_by.append("eggnog(candidate)")
            if tool_name != "kofam" and "kofam" in active_tools and (called_kos & kf_resc):
                rescued_by.append("kofam(rescued)")

            dropped_en = (en_raw_cands - called_kos) if (len(en_raw_cands) > 1 and bool(called_kos & en_raw_cands)) else en_dropped

            if len(called_kos) == 1:
                accepted_ko = next(iter(called_kos))
                alt_set = dropped_en - {accepted_ko}
            else:
                accepted_ko = "-"
                alt_set = called_kos | dropped_en

            if rescued_by:
                consensus_level = "single_tool_with_candidate"
                evidence_str = ",".join([tool_name] + rescued_by)
            else:
                consensus_level = "single_tool"
                evidence_str = tool_name

        else:  # num_calling_tools >= 2
            all_sets = list(tool_calls.values())
            union_all_calls = set.union(*all_sets)
            common_all = set.intersection(*all_sets)

            if common_all:
                is_unanimous = (num_calling_tools == len(active_tools))
                consensus_level = "unanimous" if is_unanimous else "majority"
                evidence_str = ",".join(sorted(tool_calls.keys()))
                dropped_en = (en_raw_cands - common_all) if (len(en_raw_cands) > 1 and bool(common_all & en_raw_cands)) else en_dropped
                unselected_minority = union_all_calls - common_all
                extra_alts = unselected_minority | dropped_en

                if len(common_all) == 1:
                    accepted_ko = next(iter(common_all))
                    alt_set = extra_alts - {accepted_ko}
                else:
                    accepted_ko = "-"
                    alt_set = common_all | extra_alts
            else:
                # Pairwise majority check
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
                    consensus_level = "majority"
                    evidence_str = ",".join(sorted(agreeing_tools))
                    dropped_en = (en_raw_cands - pairwise_agreed) if (len(en_raw_cands) > 1 and bool(pairwise_agreed & en_raw_cands)) else en_dropped
                    unselected_minority = union_all_calls - pairwise_agreed
                    extra_alts = unselected_minority | dropped_en

                    if len(pairwise_agreed) == 1:
                        accepted_ko = next(iter(pairwise_agreed))
                        alt_set = extra_alts - {accepted_ko}
                    else:
                        accepted_ko = "-"
                        alt_set = pairwise_agreed | extra_alts
                else:
                    # Disjoint conflict
                    all_alts = union_all_calls | en_dropped
                    if conflict_strategy == "priority":
                        top_tool = next((t for t in priority_order if t in tool_calls), tools_list[0])
                        top_kos = tool_calls[top_tool]
                        dropped_en_top = (en_raw_cands - top_kos) if (len(en_raw_cands) > 1 and bool(top_kos & en_raw_cands)) else en_dropped
                        if len(top_kos) == 1:
                            accepted_ko = next(iter(top_kos))
                            alt_set = ((union_all_calls - top_kos) | dropped_en_top) - {accepted_ko}
                        else:
                            accepted_ko = "-"
                            alt_set = all_alts
                        consensus_level = "conflict_priority"
                        evidence_str = top_tool
                    else:  # "multiple" / "union"
                        accepted_ko = "-"
                        alt_set = all_alts
                        consensus_level = "conflict"
                        evidence_str = ",".join(sorted(tool_calls.keys()))

        alternative_kos = format_kos(alt_set)

        # Definitions resolved strictly by KO
        accepted_set = {accepted_ko} if accepted_ko != "-" else set()
        definition = resolve_definition_for_kos(accepted_set, master_ko_defs, evidence_defs)
        alternative_definition = resolve_definition_for_kos(alt_set, master_ko_defs, evidence_defs)

        # Populate method-specific summary fields
        kf_recs = tool_recs.get("kofam", [])
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
                kf_assign_val = ",".join(f"{r['ko']}:{r['assignment']}" for r in kf_recs)
                kf_bs_val = ",".join(
                    f"{r['ko']}:{select_kofam_scores(r.get('score_type'), r.get('bit_score'), r.get('e_value'), r.get('domain_bit_score'), r.get('domain_e_value'))[0]}"
                    for r in kf_recs
                )
                kf_ev_val = ",".join(
                    f"{r['ko']}:{select_kofam_scores(r.get('score_type'), r.get('bit_score'), r.get('e_value'), r.get('domain_bit_score'), r.get('domain_e_value'))[1]}"
                    for r in kf_recs
                )
                kf_th_val = ",".join(f"{r['ko']}:{r['threshold']}" for r in kf_recs)
        else:
            kofam_ko_val = "-"
            kf_st_val = "-"
            kf_assign_val = "-"
            kf_bs_val = np.nan
            kf_ev_val = np.nan
            kf_th_val = np.nan

        dk_recs = tool_recs.get("deepkoala", [])
        if dk_recs:
            dk_thresh_kos = {r["ko"] for r in dk_recs if r["original_status"] == "threshold_passing"}
            dk_cands_kos = {r["ko"] for r in dk_recs if r["original_status"] in ("threshold_passing", "below_threshold")}
            dk_ko_val = format_kos(dk_thresh_kos)
            dk_cand_val = format_kos(dk_cands_kos)
            best_dk = max(dk_recs, key=lambda r: (r["deepkoala_probability"] if pd.notna(r["deepkoala_probability"]) else -1.0))
            dk_score_val = best_dk["deepkoala_probability"]
            dk_th_val = best_dk["deepkoala_threshold"]
        else:
            dk_ko_val = "-"
            dk_cand_val = "-"
            dk_score_val = np.nan
            dk_th_val = np.nan

        en_recs = tool_recs.get("eggnog", [])
        if en_recs:
            en_cands_kos = {r["ko"] for r in en_recs if r["original_status"] in ("threshold_passing", "below_threshold")}
            en_ko_val = format_kos(en_agreed)
            en_cand_val = format_kos(en_cands_kos)
            valid_en = [r for r in en_recs if r["original_status"] in ("threshold_passing", "below_threshold")]
            best_en = max(valid_en, key=lambda r: (r["bit_score"], -r["e_value"])) if valid_en else en_recs[0]
            en_bs_val = best_en["bit_score"]
            en_ev_val = best_en["e_value"]
        else:
            en_ko_val = "-"
            en_cand_val = "-"
            en_bs_val = np.nan
            en_ev_val = np.nan

        gene_summary_rows.append({
            "gene_id": gid,
            "accepted_ko": accepted_ko,
            "ko": accepted_ko,
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
        })

        # Determine cross_method_status and supporting_methods on each per-KO record
        for r in recs:
            ko = r["ko"]
            method = r["method"]
            orig_status = r["original_status"]

            # Supporting methods for this KO
            supporting = []
            for other_t in active_tools:
                if other_t != method:
                    other_t_recs = tool_recs.get(other_t, [])
                    other_matching = [
                        o for o in other_t_recs
                        if o["ko"] == ko and o["original_status"] in ("threshold_passing", "heuristic_rescued", "below_threshold")
                    ]
                    if other_matching:
                        # Only include if other method actively supported this KO in consensus
                        if ko in accepted_set:
                            supporting.append(other_t)
                        elif any(o["original_status"] == "threshold_passing" for o in other_matching):
                            supporting.append(other_t)

            supp_str = ",".join(sorted(supporting)) if supporting else "-"

            # Cross method status
            if orig_status == "unknown":
                cross_status = "unknown"
            elif method == "eggnog" and ko in en_dropped:
                cross_status = "disambiguated_dropped"
            elif ko in accepted_set:
                if orig_status == "threshold_passing":
                    if method == "eggnog" and r.get("shared_seed_hit", False):
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
                if consensus_level.startswith("conflict"):
                    cross_status = "conflict"
                elif r.get("shared_seed_hit", False) and len(alt_set) > 1 and accepted_ko == "-":
                    cross_status = "unresolved_multi_ko"
                elif orig_status == "threshold_passing":
                    cross_status = "alternative"
                elif orig_status == "heuristic_rescued":
                    cross_status = "heuristic_rescued"
                elif orig_status == "below_threshold":
                    if consensus_level == "orthogonal_dual_candidate":
                        cross_status = "rescued"
                    else:
                        cross_status = "below_threshold_unrescued"
                else:
                    cross_status = "unknown"
            else:
                if orig_status == "below_threshold":
                    cross_status = "below_threshold_unrescued"
                elif orig_status == "threshold_passing":
                    cross_status = "unselected"
                else:
                    cross_status = "unselected"

            adjudicated_records.append({
                "gene_id": gid,
                "ko": ko,
                "method": method,
                "original_status": orig_status,
                "bit_score": r.get("bit_score", np.nan),
                "e_value": r.get("e_value", np.nan),
                "domain_bit_score": r.get("domain_bit_score", np.nan),
                "domain_e_value": r.get("domain_e_value", np.nan),
                "score_type": r.get("score_type", "-"),
                "threshold": r.get("threshold", np.nan),
                "deepkoala_probability": r.get("deepkoala_probability", np.nan),
                "deepkoala_threshold": r.get("deepkoala_threshold", np.nan),
                "shared_seed_hit": bool(r.get("shared_seed_hit", False)),
                "cross_method_status": cross_status,
                "supporting_methods": supp_str,
            })

    # Build long-form evidence table
    if adjudicated_records:
        evidence_df = pd.DataFrame(adjudicated_records)[EVIDENCE_COLUMNS]
        evidence_df = evidence_df.sort_values(by=["gene_id", "ko", "method"]).reset_index(drop=True)
    else:
        evidence_df = pd.DataFrame(columns=EVIDENCE_COLUMNS)

    # Define final ordered gene-level schema
    final_cols = [
        "gene_id",
        "accepted_ko",
        "ko",
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

    if gene_summary_rows:
        merged_summary = pd.DataFrame(gene_summary_rows)
        output_cols = [c for c in final_cols if c in merged_summary.columns]
        result_df = merged_summary[output_cols].copy()
    else:
        result_df = pd.DataFrame(columns=final_cols)

    result_df.attrs["evidence"] = evidence_df

    # Export gene-level TSV
    if output_tsv:
        out_path = Path(output_tsv).expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_export = result_df.copy()
        for col in tsv_export.columns:
            if pd.api.types.is_float_dtype(tsv_export[col]):
                tsv_export[col] = tsv_export[col].apply(lambda x: "-" if pd.isna(x) else str(x))
            else:
                tsv_export[col] = tsv_export[col].fillna("-")
        tsv_export.to_csv(out_path, sep="\t", index=False)

        # Automatically export evidence TSV alongside annotations if not separately specified
        if evidence_tsv is None:
            evidence_tsv = out_path.parent / "kolach_evidence.tsv"

    # Export evidence TSV
    if evidence_tsv:
        ev_path = Path(evidence_tsv).expanduser().resolve()
        ev_path.parent.mkdir(parents=True, exist_ok=True)
        ev_export = evidence_df.copy()
        for col in ev_export.columns:
            if pd.api.types.is_float_dtype(ev_export[col]):
                ev_export[col] = ev_export[col].apply(lambda x: "-" if pd.isna(x) else str(x))
            elif pd.api.types.is_bool_dtype(ev_export[col]):
                ev_export[col] = ev_export[col].apply(lambda x: "True" if x else "False")
            else:
                ev_export[col] = ev_export[col].fillna("-")
        ev_export.to_csv(ev_path, sep="\t", index=False)

    return result_df


def main():
    """CLI entrypoint for standalone integration."""
    parser = argparse.ArgumentParser(
        prog="kolach-integrate",
        description="Integrate KOfam, DeepKOALA, and eggNOG annotation tables into a unified consensus table.",
    )
    parser.add_argument("--protein-fasta", type=str, default=None, help="Input protein FASTA file.")
    parser.add_argument("--kofam-table", type=str, default=None, help="KOfam annotations TSV.")
    parser.add_argument("--deepkoala-table", type=str, default=None, help="DeepKOALA annotations TSV.")
    parser.add_argument("--eggnog-table", type=str, default=None, help="eggNOG annotations TSV.")
    parser.add_argument("--output-file", type=str, required=True, help="Path for integrated output TSV.")
    parser.add_argument(
        "--evidence-file",
        type=str,
        default=None,
        help="Path for long-form evidence table TSV (default: kolach_evidence.tsv alongside output-file).",
    )
    parser.add_argument("--database-dir", type=str, default=None, help="Database directory (for ko_list definitions).")
    parser.add_argument("--eggnog-min-bitscore", type=float, default=60.0, help="Minimum eggNOG bitscore (default: 60.0).")
    parser.add_argument("--eggnog-max-evalue", type=float, default=1e-5, help="Maximum eggNOG e-value (default: 1e-5).")
    parser.add_argument(
        "--eggnog-filter-multi",
        choices=["disambiguate", "strict", "none"],
        default="disambiguate",
        help="eggNOG multi-KO filtering strategy (default: disambiguate).",
    )
    parser.add_argument(
        "--conflict-strategy",
        choices=["multiple", "priority", "drop", "union"],
        default="multiple",
        help=(
            "Consensus conflict strategy for disjoint calls: multiple (default: sets accepted_ko and ko to '-', "
            "records conflicting alternatives in alternative_kos; recommended for downstream pathway tools to avoid "
            "false multifunctional enzyme inference), priority (selects top method in hierarchy), or drop (discards "
            "conflicting calls, setting accepted_ko and ko to '-', retaining alternatives in alternative_kos)."
        ),
    )

    args = parser.parse_args()
    integrate_annotations(
        protein_fasta=args.protein_fasta,
        kofam_tsv=args.kofam_table,
        deepkoala_tsv=args.deepkoala_table,
        eggnog_tsv=args.eggnog_table,
        output_tsv=args.output_file,
        evidence_tsv=args.evidence_file,
        database_dir=args.database_dir,
        eggnog_min_bitscore=args.eggnog_min_bitscore,
        eggnog_max_evalue=args.eggnog_max_evalue,
        eggnog_filter_multi=args.eggnog_filter_multi,
        conflict_strategy=args.conflict_strategy,
    )


if __name__ == "__main__":
    main()
