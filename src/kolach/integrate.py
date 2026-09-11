"""Consensus and integration of annotation tables for kolach.

This module integrates results from KOfam, DeepKOALA, and eggNOG into a
unified, high-confidence annotation table using pandas, preserving all
method-specific metrics (bit scores, E-values, probabilities, and thresholds)
and evidence provenance.
"""
import argparse
from pathlib import Path
import re
from typing import Optional, Union

import numpy as np
import pandas as pd

KO_REGEX = re.compile(r"K\d{5}")

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
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                gene_ids.append(line[1:].split()[0])
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
                if len(parts) >= 2 and parts[0].startswith("K"):
                    ko = parts[0].strip()
                    definition = parts[-1].strip() if len(parts) > 1 else ""
                    definitions[ko] = definition
    except Exception:
        pass
    return definitions


def extract_kofam_records(tsv_path: Union[str, Path]) -> list[dict]:
    """Extract per-hit KO records from KOfam annotations TSV.

    Duplicate hits for a (gene_id, ko) pair are deterministically resolved by
    selecting the hit with the highest bit_score (and lowest e_value).
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
    df = df[df["ko"].str.match(r"^K\d{5}$", na=False)].copy()
    if df.empty:
        return []

    # Deterministic deduplication: highest bit_score, then lowest e_value
    df = df.sort_values(by=["gene_id", "ko", "bit_score", "e_value"], ascending=[True, True, False, True])
    df = df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in df.iterrows():
        assign = str(row["assignment"]).strip().lower()
        if assign == "threshold":
            orig_status = "threshold_passing"
        elif assign == "rescued":
            orig_status = "heuristic_rescued"
        else:
            orig_status = "below_threshold"

        records.append({
            "gene_id": row["gene_id"],
            "ko": row["ko"],
            "method": "kofam",
            "original_status": orig_status,
            "bit_score": row["bit_score"],
            "e_value": row["e_value"],
            "domain_bit_score": row["domain_bit_score"],
            "domain_e_value": row["domain_e_value"],
            "score_type": row["score_type"] if pd.notna(row["score_type"]) else "-",
            "threshold": row["threshold"],
            "deepkoala_probability": np.nan,
            "deepkoala_threshold": np.nan,
            "shared_seed_hit": False,
            "definition": row["definition"] if pd.notna(row["definition"]) else "-",
            "assignment": row["assignment"],
        })
    return records


def load_kofam(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load and aggregate KOfam annotations table into a gene-level summary.

    Expected columns in input: gene_id, ko, assignment, score_type, threshold,
                              bit_score, e_value, domain_bit_score, domain_e_value, definition.
    Per-KO statuses and scores are preserved without attaching one KO's status
    or score to other KOs.
    """
    records = extract_kofam_records(tsv_path)
    keep_cols = [
        "gene_id", "kofam_ko", "kofam_assignment", "kofam_bit_score",
        "kofam_evalue", "kofam_threshold", "kofam_definition"
    ]
    if not records:
        return pd.DataFrame(columns=keep_cols)

    rec_df = pd.DataFrame(records)
    gene_rows = []
    for gene_id, group in rec_df.groupby("gene_id", sort=False):
        group = group.sort_values(by=["bit_score"], ascending=False)
        kos = sorted(set(group["ko"]))
        ko_str = format_kos(set(kos))

        if len(group) == 1:
            row = group.iloc[0]
            gene_rows.append({
                "gene_id": gene_id,
                "kofam_ko": ko_str,
                "kofam_assignment": str(row["assignment"]),
                "kofam_bit_score": row["bit_score"],
                "kofam_evalue": row["e_value"],
                "kofam_threshold": row["threshold"],
                "kofam_definition": str(row["definition"]),
            })
        else:
            # Multi-KO gene: map each KO unambiguously to its own assignment and scores
            assign_str = ",".join(f"{r['ko']}:{r['assignment']}" for _, r in group.iterrows())
            score_str = ",".join(f"{r['ko']}:{r['bit_score']}" for _, r in group.iterrows())
            evalue_str = ",".join(f"{r['ko']}:{r['e_value']}" for _, r in group.iterrows())
            thresh_str = ",".join(f"{r['ko']}:{r['threshold']}" for _, r in group.iterrows())
            defs = [str(r["definition"]) for _, r in group.iterrows() if str(r["definition"]) not in ("", "-")]
            def_str = "; ".join(defs) if defs else "-"

            gene_rows.append({
                "gene_id": gene_id,
                "kofam_ko": ko_str,
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
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return []

    df.columns = [c.lstrip("#").strip() for c in df.columns]
    cols_lower = [c.lower() for c in df.columns]
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

    # Determine validity mask
    valid_mask = pd.Series(True, index=df.index)
    if "deepkoala_annotate" in df.columns:
        valid_mask = valid_mask & df["deepkoala_annotate"].astype(str).str.contains(r"\*")
    elif df["deepkoala_score"].notna().any() and df["deepkoala_threshold"].notna().any():
        has_scores = df["deepkoala_score"].notna() & df["deepkoala_threshold"].notna()
        valid_mask = valid_mask & (~has_scores | (df["deepkoala_score"] >= df["deepkoala_threshold"]))

    df["is_threshold_passing"] = valid_mask

    exploded_rows = []
    for _, row in df.iterrows():
        raw_val = row.get("predict_label", "-")
        kos = parse_kos(raw_val)
        for ko in kos:
            exploded_rows.append({
                "gene_id": row["gene_id"],
                "ko": ko,
                "deepkoala_score": row["deepkoala_score"],
                "deepkoala_threshold": row["deepkoala_threshold"],
                "is_threshold_passing": row["is_threshold_passing"],
            })

    if not exploded_rows:
        return []

    exp_df = pd.DataFrame(exploded_rows)
    # Deduplicate by highest score
    exp_df = exp_df.sort_values(by=["gene_id", "ko", "deepkoala_score"], ascending=[True, True, False])
    exp_df = exp_df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in exp_df.iterrows():
        orig_status = "threshold_passing" if row["is_threshold_passing"] else "below_threshold"
        records.append({
            "gene_id": row["gene_id"],
            "ko": row["ko"],
            "method": "deepkoala",
            "original_status": orig_status,
            "bit_score": np.nan,
            "e_value": np.nan,
            "domain_bit_score": np.nan,
            "domain_e_value": np.nan,
            "score_type": "-",
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
        cands = sorted(set(group["ko"]))
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


def extract_eggnog_records(
    tsv_path: Union[str, Path],
    min_bitscore: float = 60.0,
    max_evalue: float = 1e-5,
) -> list[dict]:
    """Extract per-hit KO records from eggNOG annotations TSV.

    Shared seed-hit metrics for multi-KO assignments are explicitly tagged.
    Duplicate hits for a (gene_id, ko) pair are deterministically resolved by
    selecting the hit with the highest bit_score (and lowest e_value).
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []

    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if df.empty:
        return []

    df.columns = [c.lstrip("#").strip() for c in df.columns]
    cols_lower = [c.lower() for c in df.columns]
    col_map = {}
    if "query" in cols_lower:
        col_map[df.columns[cols_lower.index("query")]] = "gene_id"
    elif "gene_id" in cols_lower:
        col_map[df.columns[cols_lower.index("gene_id")]] = "gene_id"
    else:
        col_map[df.columns[0]] = "gene_id"

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
        for ko in kos:
            exploded_rows.append({
                "gene_id": row["gene_id"],
                "ko": ko,
                "bit_score": row["eggnog_bit_score"],
                "e_value": row["eggnog_evalue"],
                "shared_seed_hit": is_multi,
                "definition": row["eggnog_description"] if pd.notna(row["eggnog_description"]) else "-",
            })

    if not exploded_rows:
        return []

    exp_df = pd.DataFrame(exploded_rows)
    # Deduplicate: highest bit_score, then lowest e_value
    exp_df = exp_df.sort_values(by=["gene_id", "ko", "bit_score", "e_value"], ascending=[True, True, False, True])
    exp_df = exp_df.drop_duplicates(subset=["gene_id", "ko"], keep="first")

    records = []
    for _, row in exp_df.iterrows():
        bs = row["bit_score"]
        ev = row["e_value"]
        is_thresh = (pd.notna(bs) and bs >= min_bitscore) and (pd.isna(ev) or ev <= max_evalue)
        orig_status = "threshold_passing" if is_thresh else "below_threshold"

        records.append({
            "gene_id": row["gene_id"],
            "ko": row["ko"],
            "method": "eggnog",
            "original_status": orig_status,
            "bit_score": bs,
            "e_value": ev,
            "domain_bit_score": np.nan,
            "domain_e_value": np.nan,
            "score_type": "seed_hit",
            "threshold": min_bitscore,
            "deepkoala_probability": np.nan,
            "deepkoala_threshold": np.nan,
            "shared_seed_hit": bool(row["shared_seed_hit"]),
            "definition": row["definition"],
        })
    return records


def load_eggnog(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load raw eggNOG-mapper annotations table."""
    path = Path(tsv_path).expanduser().resolve()
    keep_cols = ["gene_id", "eggnog_raw_ko", "eggnog_candidate_ko", "eggnog_bit_score", "eggnog_evalue", "eggnog_description"]
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=keep_cols)

    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if df.empty:
        return pd.DataFrame(columns=keep_cols)

    df.columns = [c.lstrip("#").strip() for c in df.columns]
    cols_lower = [c.lower() for c in df.columns]
    col_map = {}
    if "query" in cols_lower:
        col_map[df.columns[cols_lower.index("query")]] = "gene_id"
    elif "gene_id" in cols_lower:
        col_map[df.columns[cols_lower.index("gene_id")]] = "gene_id"
    else:
        col_map[df.columns[0]] = "gene_id"

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

    df["eggnog_candidate_ko"] = df["eggnog_raw_ko"].apply(lambda v: format_kos(parse_kos(v)))
    df = df.sort_values(by=["gene_id", "eggnog_bit_score"], ascending=[True, False])
    return df[keep_cols].drop_duplicates(subset=["gene_id"])


def extract_kofam_status_sets(kofam_ko: str, kofam_assignment: str) -> tuple[set[str], set[str]]:
    """Separate confident threshold-passing KOs from heuristic-rescued KOs."""
    kos = parse_kos(kofam_ko)
    if not kos:
        return set(), set()
    assign_str = str(kofam_assignment).strip()
    if not assign_str or assign_str == "-":
        return kos, set()

    if ":" in assign_str:
        thresh_kos = set()
        rescued_kos = set()
        for part in assign_str.split(","):
            if ":" in part:
                k, a = part.split(":", 1)
                k = k.strip()
                a = a.strip().lower()
                if a == "threshold":
                    thresh_kos.add(k)
                elif a == "rescued":
                    rescued_kos.add(k)
                else:
                    thresh_kos.add(k)
        return thresh_kos, rescued_kos

    if assign_str.lower() == "threshold":
        return kos, set()
    if assign_str.lower() == "rescued":
        return set(), kos
    return kos, set()


def filter_and_disambiguate_eggnog(
    eggnog_df: pd.DataFrame,
    merged_df: pd.DataFrame,
    min_bitscore: float = 60.0,
    max_evalue: float = 1e-5,
    filter_multi_mode: str = "disambiguate",
) -> pd.DataFrame:
    """Filter eggNOG hits by bitscore and disambiguate multi-KO assignments.

    Crucially, below-threshold eggNOG hits are NEVER promoted to confident calls
    in eggnog_ko. Instead, they remain in eggnog_candidate_ko and are evaluated
    strictly as candidate evidence during consensus adjudication.
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

        # 1. Below threshold: never promote to confident eggnog_ko!
        if not en_kos or pd.isna(bitscore) or bitscore < min_bitscore or (pd.notna(evalue) and evalue > max_evalue):
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
    Any additional alternatives (e.g. sibling candidate hits dropped when other
    tools disambiguate eggNOG multi-KO hits) are preserved in alternative_kos.
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


def adjudicate_consensus(
    merged_df: pd.DataFrame,
    active_tools: list[str],
    conflict_strategy: str = "multiple",
    ko_definitions: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    """Evaluate consensus KO, consensus level, evidence, and functional definition.

    - Computes confident agreement using ONLY calls that independently passed
      thresholds (KOfam heuristic rescues and sub-threshold eggNOG/DeepKOALA
      candidates do NOT count as confident votes).
    - Labeled 'single_tool_with_candidate' when 1 confident tool is corroborated
      by sub-threshold candidate evidence (never 'majority' or 'unanimous').
    - Labeled 'orthogonal_dual_candidate' when sub-threshold candidates from two
      methods agree.
    - Separates accepted single KOs (accepted_ko, ko) from conflicting or multi-KO
      alternatives (alternative_kos) and functional definitions (definition, alternative_definition).
      Strict single-hit policy: accepted_ko contains single hits only.
    """
    if ko_definitions is None:
        ko_definitions = {}

    accepted_kos = []
    alternative_kos_list = []
    consensus_levels = []
    evidence_list = []
    definitions = []
    alt_definitions = []

    priority_order = [t for t in ["kofam", "deepkoala", "eggnog"] if t in active_tools]

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

        # Sub-threshold candidates
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
                consensus_levels.append("single_tool_with_candidate")
                evidence_list.append(",".join([tool_name] + rescued_by))
            else:
                consensus_levels.append("single_tool")
                evidence_list.append(tool_name)
            continue

        # Case 2: 2 or more confident tools made calls
        all_sets = list(tool_calls.values())
        common_all = set.intersection(*all_sets)

        if common_all:
            is_unanimous = (num_calling_tools == len(active_tools))
            dropped_en = (en_raw - common_all) if (len(en_raw) > 1 and bool(common_all & en_raw)) else set()
            acc, alt, defn, alt_defn = _assign_accepted_or_alt(common_all, row, ko_definitions, additional_alternatives=dropped_en)
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
            acc, alt, defn, alt_defn = _assign_accepted_or_alt(pairwise_agreed, row, ko_definitions, additional_alternatives=dropped_en)
            accepted_kos.append(acc)
            alternative_kos_list.append(alt)
            consensus_levels.append("majority")
            evidence_list.append(",".join(sorted(agreeing_tools)))
            definitions.append(defn)
            alt_definitions.append(alt_defn)
            continue

        # Disjoint conflict: zero overlap between confident tools
        union_kos = set.union(*all_sets)
        dropped_en = (en_raw - union_kos) if (len(en_raw) > 1 and bool(union_kos & en_raw)) else set()
        all_alts = union_kos | dropped_en

        if conflict_strategy == "drop":
            accepted_kos.append("-")
            alternative_kos_list.append(format_kos(all_alts))
            consensus_levels.append("conflict_dropped")
            evidence_list.append("-")
            definitions.append("-")
            alt_definitions.append(_resolve_definition(all_alts, row, ko_definitions))
        elif conflict_strategy == "priority":
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


def _resolve_definition(
    kos: set[str],
    row: pd.Series,
    ko_definitions: dict[str, str],
) -> str:
    """Find the best definition text for the specified KO(s)."""
    if not kos:
        return "-"
    found_defs = []
    for ko in sorted(kos):
        if ko in ko_definitions and ko_definitions[ko]:
            found_defs.append(ko_definitions[ko])
    if found_defs:
        return "; ".join(found_defs)

    kf_def = row.get("kofam_definition")
    if pd.notna(kf_def) and kf_def and kf_def != "-":
        return str(kf_def)
    en_def = row.get("eggnog_description")
    if pd.notna(en_def) and en_def and en_def != "-":
        return str(en_def)
    return "-"


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
    """Main pandas integration pipeline for kolach annotation outputs.

    Produces:
    1. Gene-level annotations table (output_tsv, default kolach_annotations.tsv)
    2. Long-form per-KO evidence table (evidence_tsv, default kolach_evidence.tsv)
    """
    active_tools = []
    dfs_to_merge = {}
    all_evidence_records = []

    # 1. Ingest KOfam
    if kofam_tsv and Path(kofam_tsv).exists():
        kf_records = extract_kofam_records(kofam_tsv)
        all_evidence_records.extend(kf_records)
        kf_df = load_kofam(kofam_tsv)
        dfs_to_merge["kofam"] = kf_df
        active_tools.append("kofam")

    # 2. Ingest DeepKOALA
    if deepkoala_tsv and Path(deepkoala_tsv).exists():
        dk_records = extract_deepkoala_records(deepkoala_tsv)
        all_evidence_records.extend(dk_records)
        dk_df = load_deepkoala(deepkoala_tsv)
        dfs_to_merge["deepkoala"] = dk_df
        active_tools.append("deepkoala")

    # 3. Ingest eggNOG
    en_df = None
    if eggnog_tsv and Path(eggnog_tsv).exists():
        en_records = extract_eggnog_records(
            eggnog_tsv,
            min_bitscore=eggnog_min_bitscore,
            max_evalue=eggnog_max_evalue,
        )
        all_evidence_records.extend(en_records)
        en_df = load_eggnog(eggnog_tsv)
        active_tools.append("eggnog")

    # Determine universe and order of gene IDs
    if protein_fasta and Path(protein_fasta).exists():
        all_gene_ids = read_fasta_ids(protein_fasta)
        base_df = pd.DataFrame({"gene_id": all_gene_ids})
    else:
        all_ids = set()
        for df_item in dfs_to_merge.values():
            all_ids.update(df_item["gene_id"])
        if en_df is not None:
            all_ids.update(en_df["gene_id"])
        base_df = pd.DataFrame({"gene_id": sorted(all_ids)})

    # Merge KOfam & DeepKOALA
    merged = base_df.copy()
    if "kofam" in dfs_to_merge:
        merged = merged.merge(dfs_to_merge["kofam"], on="gene_id", how="left")
    else:
        merged["kofam_ko"] = "-"
        merged["kofam_assignment"] = "-"
        merged["kofam_bit_score"] = np.nan
        merged["kofam_evalue"] = np.nan
        merged["kofam_threshold"] = np.nan
        merged["kofam_definition"] = "-"

    if "deepkoala" in dfs_to_merge:
        merged = merged.merge(dfs_to_merge["deepkoala"], on="gene_id", how="left")
    else:
        merged["deepkoala_ko"] = "-"
        merged["deepkoala_candidate_ko"] = "-"
        merged["deepkoala_score"] = np.nan
        merged["deepkoala_threshold"] = np.nan

    merged["kofam_ko"] = merged["kofam_ko"].fillna("-")
    merged["kofam_assignment"] = merged["kofam_assignment"].fillna("-")
    merged["kofam_definition"] = merged["kofam_definition"].fillna("-")
    merged["deepkoala_ko"] = merged["deepkoala_ko"].fillna("-")
    merged["deepkoala_candidate_ko"] = merged["deepkoala_candidate_ko"].fillna("-")

    # Merge & filter eggNOG
    if en_df is not None:
        merged = filter_and_disambiguate_eggnog(
            en_df,
            merged,
            min_bitscore=eggnog_min_bitscore,
            max_evalue=eggnog_max_evalue,
            filter_multi_mode=eggnog_filter_multi,
        )
    else:
        merged["eggnog_ko"] = "-"
        merged["eggnog_candidate_ko"] = "-"
        merged["eggnog_bit_score"] = np.nan
        merged["eggnog_evalue"] = np.nan
        merged["eggnog_description"] = "-"

    merged["eggnog_ko"] = merged["eggnog_ko"].fillna("-")
    merged["eggnog_candidate_ko"] = merged["eggnog_candidate_ko"].fillna("-")
    merged["eggnog_description"] = merged["eggnog_description"].fillna("-")

    # Load external KO definitions if available
    ko_list_file = None
    if database_dir:
        candidate = Path(database_dir).expanduser().resolve() / "kofam" / "ko_list"
        if candidate.is_file():
            ko_list_file = candidate
    ko_defs = load_ko_definitions(ko_list_file)

    # Adjudicate consensus
    merged = adjudicate_consensus(
        merged,
        active_tools=active_tools,
        conflict_strategy=conflict_strategy,
        ko_definitions=ko_defs,
    )

    # Build long-form evidence table
    evidence_rows = []
    # Map consensus results by gene_id for cross-method status evaluation
    consensus_map = merged.set_index("gene_id").to_dict(orient="index")

    for rec in all_evidence_records:
        gid = rec["gene_id"]
        ko = rec["ko"]
        method = rec["method"]
        orig_status = rec["original_status"]
        gene_info = consensus_map.get(gid, {})

        accepted_set = parse_kos(gene_info.get("accepted_ko", "-"))
        alt_set = parse_kos(gene_info.get("alternative_kos", "-"))
        consensus_level = gene_info.get("consensus_level", "unannotated")

        # Determine supporting methods for this KO
        supporting = []
        if orig_status == "threshold_passing":
            if ko in accepted_set:
                cross_status = "accepted"
                if method == "eggnog" and rec.get("shared_seed_hit", False):
                    cross_status = "disambiguated_retained"
                for other_tool in active_tools:
                    if other_tool != method:
                        other_kos = parse_kos(gene_info.get(f"{other_tool}_ko", "-"))
                        if ko in other_kos:
                            supporting.append(other_tool)
            elif method == "eggnog" and ko not in parse_kos(gene_info.get("eggnog_ko", "-")):
                cross_status = "disambiguated_dropped"
                for other_tool in active_tools:
                    if other_tool != "eggnog":
                        other_kos = parse_kos(gene_info.get(f"{other_tool}_ko", "-"))
                        if other_kos:
                            supporting.append(other_tool)
            elif ko in alt_set or consensus_level.startswith("conflict"):
                if consensus_level.startswith("conflict"):
                    cross_status = "conflict"
                elif rec.get("shared_seed_hit", False):
                    cross_status = "unresolved_multi_ko"
                else:
                    cross_status = "alternative"
                for other_tool in active_tools:
                    if other_tool != method:
                        other_kos = parse_kos(gene_info.get(f"{other_tool}_ko", "-"))
                        if ko in other_kos:
                            supporting.append(other_tool)
            else:
                cross_status = "unselected"
        elif orig_status == "heuristic_rescued":
            if ko in accepted_set:
                cross_status = "heuristic_rescued"
            elif ko in alt_set:
                cross_status = "unresolved_multi_ko" if len(alt_set) > 1 else "heuristic_rescued"
            else:
                cross_status = "heuristic_rescued"
            for other_tool in active_tools:
                if other_tool != method:
                    other_kos = parse_kos(gene_info.get(f"{other_tool}_ko", "-"))
                    if ko in other_kos:
                        supporting.append(other_tool)
        else:  # below_threshold
            if ko in accepted_set:
                cross_status = "rescued"
            elif ko in alt_set and consensus_level == "orthogonal_dual_candidate":
                cross_status = "rescued"
            else:
                cross_status = "below_threshold_unrescued"
            for other_tool in active_tools:
                if other_tool != method:
                    other_kos = parse_kos(gene_info.get(f"{other_tool}_ko", "-"))
                    if ko in other_kos:
                        supporting.append(other_tool)

        supp_str = ",".join(sorted(supporting)) if supporting else "-"

        evidence_rows.append({
            "gene_id": gid,
            "ko": ko,
            "method": method,
            "original_status": orig_status,
            "bit_score": rec.get("bit_score", np.nan),
            "e_value": rec.get("e_value", np.nan),
            "domain_bit_score": rec.get("domain_bit_score", np.nan),
            "domain_e_value": rec.get("domain_e_value", np.nan),
            "score_type": rec.get("score_type", "-"),
            "threshold": rec.get("threshold", np.nan),
            "deepkoala_probability": rec.get("deepkoala_probability", np.nan),
            "deepkoala_threshold": rec.get("deepkoala_threshold", np.nan),
            "shared_seed_hit": bool(rec.get("shared_seed_hit", False)),
            "cross_method_status": cross_status,
            "supporting_methods": supp_str,
        })

    if evidence_rows:
        evidence_df = pd.DataFrame(evidence_rows)[EVIDENCE_COLUMNS]
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

    output_cols = [c for c in final_cols if c in merged.columns]
    result_df = merged[output_cols].copy()
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
