"""Consensus and integration of annotation tables for kolach.

This module integrates results from KOfam, DeepKOALA, and eggNOG into a
unified, high-confidence annotation table using pandas, preserving all
method-specific metrics (bit scores, E-values, probabilities, and thresholds).
"""
import argparse
from pathlib import Path
import re
from typing import Optional, Union

import numpy as np
import pandas as pd

KO_REGEX = re.compile(r"K\d{5}")


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
                    # Definition is typically the last non-threshold field
                    definition = parts[-1].strip() if len(parts) > 1 else ""
                    definitions[ko] = definition
    except Exception:
        pass
    return definitions


def load_kofam(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load and aggregate KOfam annotations table.

    Expected columns: gene_id, ko, assignment, score_type, threshold,
                      bit_score, e_value, domain_bit_score, domain_e_value, definition
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=[
            "gene_id", "kofam_ko", "kofam_assignment", "kofam_bit_score",
            "kofam_evalue", "kofam_threshold", "kofam_definition"
        ])

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return pd.DataFrame(columns=[
            "gene_id", "kofam_ko", "kofam_assignment", "kofam_bit_score",
            "kofam_evalue", "kofam_threshold", "kofam_definition"
        ])

    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df["ko"] = df["ko"].astype(str).str.strip()
    df["bit_score"] = pd.to_numeric(df["bit_score"], errors="coerce")
    df["e_value"] = pd.to_numeric(df["e_value"], errors="coerce")
    df["threshold"] = pd.to_numeric(df.get("threshold", np.nan), errors="coerce")
    if "definition" not in df.columns:
        df["definition"] = "-"
    if "assignment" not in df.columns:
        df["assignment"] = "-"

    # Sort so highest bit_score comes first
    df = df.sort_values(by=["gene_id", "bit_score"], ascending=[True, False])

    # Top hit per gene provides scalar scores and definitions
    best = df.drop_duplicates(subset=["gene_id"], keep="first").copy()
    best = best.rename(columns={
        "bit_score": "kofam_bit_score",
        "e_value": "kofam_evalue",
        "threshold": "kofam_threshold",
        "assignment": "kofam_assignment",
        "definition": "kofam_definition",
    })

    # Aggregate all unique KOs for each gene
    ko_map = df.groupby("gene_id")["ko"].agg(lambda s: format_kos({k for k in s if k and k != "-"}))
    best["kofam_ko"] = best["gene_id"].map(ko_map).fillna("-")

    # If any row achieved 'threshold', prioritize 'threshold' over 'rescued'
    thresh_genes = set(df.loc[df["assignment"] == "threshold", "gene_id"])
    best.loc[best["gene_id"].isin(thresh_genes), "kofam_assignment"] = "threshold"

    keep_cols = [
        "gene_id", "kofam_ko", "kofam_assignment", "kofam_bit_score",
        "kofam_evalue", "kofam_threshold", "kofam_definition"
    ]
    return best[keep_cols]


def load_deepkoala(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load and format DeepKOALA annotations table.

    Expected columns:
      - Detail mode: name, predict_label, probability, threshold, annotate
      - Simple mode: name, predict_label
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=[
            "gene_id", "deepkoala_ko", "deepkoala_score", "deepkoala_threshold"
        ])

    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return pd.DataFrame(columns=[
            "gene_id", "deepkoala_ko", "deepkoala_score", "deepkoala_threshold"
        ])

    df.columns = [c.lstrip("#").strip() for c in df.columns]

    # Standardize column names
    col_map = {}
    cols_lower = [c.lower() for c in df.columns]
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

    if "deepkoala_score" in df.columns:
        df["deepkoala_score"] = pd.to_numeric(df["deepkoala_score"], errors="coerce")
    else:
        df["deepkoala_score"] = np.nan

    if "deepkoala_threshold" in df.columns:
        df["deepkoala_threshold"] = pd.to_numeric(df["deepkoala_threshold"], errors="coerce")
    else:
        df["deepkoala_threshold"] = np.nan

    # A DeepKOALA prediction is only trusted if probability >= threshold (or annotate contains '*')
    valid_mask = pd.Series(True, index=df.index)
    if "deepkoala_annotate" in df.columns:
        valid_mask = valid_mask & df["deepkoala_annotate"].astype(str).str.contains(r"\*")
    elif df["deepkoala_score"].notna().any() and df["deepkoala_threshold"].notna().any():
        has_scores = df["deepkoala_score"].notna() & df["deepkoala_threshold"].notna()
        valid_mask = valid_mask & (~has_scores | (df["deepkoala_score"] >= df["deepkoala_threshold"]))

    if "predict_label" in df.columns:
        raw_kos = df["predict_label"].apply(lambda v: format_kos(parse_kos(v)))
        df["deepkoala_candidate_ko"] = raw_kos
        df["deepkoala_ko"] = np.where(valid_mask, raw_kos, "-")
    else:
        df["deepkoala_candidate_ko"] = "-"
        df["deepkoala_ko"] = "-"

    keep_cols = ["gene_id", "deepkoala_ko", "deepkoala_candidate_ko", "deepkoala_score", "deepkoala_threshold"]
    return df[keep_cols].drop_duplicates(subset=["gene_id"])


def load_eggnog(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Load raw eggNOG-mapper annotations table.

    Expected columns: query, score (bitscore), evalue, KEGG_ko, Description...
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=[
            "gene_id", "eggnog_raw_ko", "eggnog_candidate_ko", "eggnog_bit_score", "eggnog_evalue", "eggnog_description"
        ])

    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if df.empty:
        return pd.DataFrame(columns=[
            "gene_id", "eggnog_raw_ko", "eggnog_candidate_ko", "eggnog_bit_score", "eggnog_evalue", "eggnog_description"
        ])

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
    if "eggnog_bit_score" in df.columns:
        df["eggnog_bit_score"] = pd.to_numeric(df["eggnog_bit_score"], errors="coerce")
    else:
        df["eggnog_bit_score"] = np.nan

    if "eggnog_evalue" in df.columns:
        df["eggnog_evalue"] = pd.to_numeric(df["eggnog_evalue"], errors="coerce")
    else:
        df["eggnog_evalue"] = np.nan

    if "eggnog_description" not in df.columns:
        df["eggnog_description"] = "-"

    df["eggnog_candidate_ko"] = df["eggnog_raw_ko"].apply(lambda v: format_kos(parse_kos(v)))

    # In case of duplicate query hits, sort by bit_score descending
    df = df.sort_values(by=["gene_id", "eggnog_bit_score"], ascending=[True, False])
    keep_cols = ["gene_id", "eggnog_raw_ko", "eggnog_candidate_ko", "eggnog_bit_score", "eggnog_evalue", "eggnog_description"]
    return df[keep_cols].drop_duplicates(subset=["gene_id"])


def filter_and_disambiguate_eggnog(
    eggnog_df: pd.DataFrame,
    merged_df: pd.DataFrame,
    min_bitscore: float = 60.0,
    filter_multi_mode: str = "disambiguate",
) -> pd.DataFrame:
    """Filter eggNOG hits by bitscore and disambiguate multi-KO assignments."""
    if eggnog_df.empty:
        merged_df["eggnog_ko"] = "-"
        merged_df["eggnog_candidate_ko"] = "-"
        merged_df["eggnog_bit_score"] = np.nan
        merged_df["eggnog_evalue"] = np.nan
        merged_df["eggnog_description"] = "-"
        return merged_df

    # Join raw eggnog cols
    merged_df = merged_df.merge(eggnog_df, on="gene_id", how="left")

    if "eggnog_raw_ko" not in merged_df.columns:
        merged_df["eggnog_raw_ko"] = "-"
    merged_df["eggnog_raw_ko"] = merged_df["eggnog_raw_ko"].fillna("-")
    if "eggnog_candidate_ko" not in merged_df.columns:
        merged_df["eggnog_candidate_ko"] = "-"
    merged_df["eggnog_candidate_ko"] = merged_df["eggnog_candidate_ko"].fillna("-")

    filtered_kos = []
    for _, row in merged_df.iterrows():
        raw_val = row.get("eggnog_raw_ko", "-")
        bitscore = row.get("eggnog_bit_score", np.nan)
        en_kos = parse_kos(raw_val)

        kf_kos = parse_kos(row.get("kofam_ko", "-"))
        dk_kos = parse_kos(row.get("deepkoala_ko", "-"))
        dk_cand = parse_kos(row.get("deepkoala_candidate_ko", "-"))
        trusted_other = kf_kos | dk_kos

        # 1. Below bitscore threshold: check candidate rescue
        if not en_kos or pd.isna(bitscore) or bitscore < min_bitscore:
            overlap = en_kos & trusted_other
            if overlap:
                # Sub-threshold eggNOG hit rescued by confident KOfam or DeepKOALA
                filtered_kos.append(format_kos(overlap))
            else:
                filtered_kos.append("-")
            continue

        # 2. Confident single-KO
        if len(en_kos) == 1 or filter_multi_mode == "none":
            filtered_kos.append(format_kos(en_kos))
            continue

        # 3. Multi-KO handling (bitscore >= min_bitscore)
        overlap = trusted_other & en_kos
        if overlap:
            # Disambiguated by confident KOfam or DeepKOALA
            filtered_kos.append(format_kos(overlap))
        elif dk_cand & en_kos:
            # Disambiguated by DeepKOALA candidate
            filtered_kos.append(format_kos(dk_cand & en_kos))
        elif trusted_other:
            # Disagreeing call against other trusted tools
            if filter_multi_mode == "strict":
                filtered_kos.append("-")
            else:
                filtered_kos.append(format_kos(en_kos))
        else:
            # No other tool called anything
            if filter_multi_mode == "strict":
                filtered_kos.append("-")
            else:
                filtered_kos.append(format_kos(en_kos))

    merged_df["eggnog_ko"] = filtered_kos
    merged_df.drop(columns=["eggnog_raw_ko"], inplace=True)
    return merged_df


def adjudicate_consensus(
    merged_df: pd.DataFrame,
    active_tools: list[str],
    conflict_strategy: str = "priority",
    ko_definitions: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    """Evaluate consensus KO, consensus level, evidence, and functional definition."""
    if ko_definitions is None:
        ko_definitions = {}

    consensus_kos = []
    consensus_levels = []
    evidence_list = []
    definitions = []

    priority_order = [t for t in ["kofam", "deepkoala", "eggnog"] if t in active_tools]

    for _, row in merged_df.iterrows():
        tool_calls: dict[str, set[str]] = {}
        for tool in active_tools:
            col = f"{tool}_ko"
            if col in row:
                kos = parse_kos(row[col])
                if kos:
                    tool_calls[tool] = kos

        num_calling_tools = len(tool_calls)

        dk_cand = parse_kos(row.get("deepkoala_candidate_ko", "-"))
        en_cand = parse_kos(row.get("eggnog_candidate_ko", "-"))

        # Case 0: No confident tool calls
        if num_calling_tools == 0:
            # Candidate rescue: Check if DeepKOALA candidate and eggNOG candidate agree
            if "deepkoala" in active_tools and "eggnog" in active_tools:
                cand_agree = dk_cand & en_cand
                if cand_agree:
                    consensus_kos.append(format_kos(cand_agree))
                    consensus_levels.append("rescued_dual_candidate")
                    evidence_list.append("deepkoala(candidate),eggnog(candidate)")
                    definitions.append(_resolve_definition(cand_agree, row, ko_definitions))
                    continue

            consensus_kos.append("-")
            consensus_levels.append("unannotated")
            evidence_list.append("-")
            definitions.append("-")
            continue

        # Case 1: Exactly 1 confident tool call
        if num_calling_tools == 1:
            tool_name = next(iter(tool_calls))
            kos = tool_calls[tool_name]

            # Check if an unconfident candidate from another active tool corroborates this KO
            rescued_by = []
            if tool_name != "deepkoala" and "deepkoala" in active_tools and (kos & dk_cand):
                rescued_by.append("deepkoala(candidate)")
            if tool_name != "eggnog" and "eggnog" in active_tools and (kos & en_cand):
                rescued_by.append("eggnog(candidate)")

            if rescued_by:
                final_ko_str = format_kos(kos)
                consensus_kos.append(final_ko_str)
                consensus_levels.append("majority_rescued")
                evidence_list.append(",".join([tool_name] + rescued_by))
                definitions.append(_resolve_definition(kos, row, ko_definitions))
            else:
                final_ko_str = format_kos(kos)
                consensus_kos.append(final_ko_str)
                consensus_levels.append("single_tool")
                evidence_list.append(tool_name)
                definitions.append(_resolve_definition(kos, row, ko_definitions))
            continue

        # Case 2: 2 or more confident tools made calls: evaluate agreement
        all_sets = list(tool_calls.values())
        common_all = set.intersection(*all_sets)

        if common_all:
            # All calling tools agree on common_all
            is_unanimous = (num_calling_tools == len(active_tools))
            consensus_kos.append(format_kos(common_all))
            consensus_levels.append("unanimous" if is_unanimous else "majority")
            evidence_list.append(",".join(sorted(tool_calls.keys())))
            definitions.append(_resolve_definition(common_all, row, ko_definitions))
            continue

        # Pairwise check for majority (e.g. 2 out of 3 agree)
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
            consensus_kos.append(format_kos(pairwise_agreed))
            consensus_levels.append("majority")
            evidence_list.append(",".join(sorted(agreeing_tools)))
            definitions.append(_resolve_definition(pairwise_agreed, row, ko_definitions))
            continue

        # Conflict: Zero overlap between calling tools
        if conflict_strategy == "drop":
            consensus_kos.append("-")
            consensus_levels.append("conflict_dropped")
            evidence_list.append("-")
            definitions.append("-")
        elif conflict_strategy == "priority":
            top_tool = next((t for t in priority_order if t in tool_calls), tools_list[0])
            top_kos = tool_calls[top_tool]
            consensus_kos.append(format_kos(top_kos))
            consensus_levels.append("conflict_priority")
            evidence_list.append(top_tool)
            definitions.append(_resolve_definition(top_kos, row, ko_definitions))
        else:  # "multiple" / "union" (default)
            union_kos = set.union(*all_sets)
            consensus_kos.append(format_kos(union_kos))
            consensus_levels.append("conflict")
            evidence_list.append(",".join(sorted(tool_calls.keys())))
            definitions.append(_resolve_definition(union_kos, row, ko_definitions))

    merged_df["ko"] = consensus_kos
    merged_df["definition"] = definitions
    merged_df["consensus_level"] = consensus_levels
    merged_df["evidence"] = evidence_list

    return merged_df


def _resolve_definition(
    kos: set[str],
    row: pd.Series,
    ko_definitions: dict[str, str]
) -> str:
    """Find the best definition text for the assigned KO(s)."""
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
    database_dir: Optional[Union[str, Path]] = None,
    eggnog_min_bitscore: float = 60.0,
    eggnog_filter_multi: str = "disambiguate",
    conflict_strategy: str = "multiple",
) -> pd.DataFrame:
    """Main pandas integration pipeline for kolach annotation outputs."""
    active_tools = []
    dfs_to_merge = {}

    # Ingest KOfam
    if kofam_tsv and Path(kofam_tsv).exists():
        kf_df = load_kofam(kofam_tsv)
        dfs_to_merge["kofam"] = kf_df
        active_tools.append("kofam")

    # Ingest DeepKOALA
    if deepkoala_tsv and Path(deepkoala_tsv).exists():
        dk_df = load_deepkoala(deepkoala_tsv)
        dfs_to_merge["deepkoala"] = dk_df
        active_tools.append("deepkoala")

    # Ingest eggNOG
    en_df = None
    if eggnog_tsv and Path(eggnog_tsv).exists():
        en_df = load_eggnog(eggnog_tsv)
        active_tools.append("eggnog")

    # Determine universe and order of gene IDs
    if protein_fasta and Path(protein_fasta).exists():
        all_gene_ids = read_fasta_ids(protein_fasta)
        base_df = pd.DataFrame({"gene_id": all_gene_ids})
    else:
        # Collect union of all gene IDs from input tables
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

    # Fill NaNs for tool KO columns before eggNOG disambiguation
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

    # Define final ordered column schema
    final_cols = [
        "gene_id",
        "ko",
        "definition",
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

    # Reorder columns that exist
    output_cols = [c for c in final_cols if c in merged.columns]
    result_df = merged[output_cols].copy()

    # Format output for TSV: replace float NaN with '-' or clean string
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
    parser.add_argument("--database-dir", type=str, default=None, help="Database directory (for ko_list definitions).")
    parser.add_argument("--eggnog-min-bitscore", type=float, default=60.0, help="Minimum eggNOG bitscore (default: 60.0).")
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
        help="Conflict resolution strategy for disjoint calls: multiple (report as conflict and list all KOs, default), priority, or drop.",
    )

    args = parser.parse_args()
    integrate_annotations(
        protein_fasta=args.protein_fasta,
        kofam_tsv=args.kofam_table,
        deepkoala_tsv=args.deepkoala_table,
        eggnog_tsv=args.eggnog_table,
        output_tsv=args.output_file,
        database_dir=args.database_dir,
        eggnog_min_bitscore=args.eggnog_min_bitscore,
        eggnog_filter_multi=args.eggnog_filter_multi,
        conflict_strategy=args.conflict_strategy,
    )


if __name__ == "__main__":
    main()
