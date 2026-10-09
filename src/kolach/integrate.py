"""Integrate KOfam, DeepKOALA and eggNOG evidence into gene annotations."""

from collections import Counter
from dataclasses import dataclass, field
import gzip
from pathlib import Path
import re
from typing import Any, Optional, Union

import numpy as np
import pandas as pd

KO_REGEX = re.compile(r"\bK\d{5}\b")

EVIDENCE_COLUMNS = [
    "gene_id", "ko", "method", "bit_score", "e_value",
    "domain_bit_score", "domain_e_value", "score_type",
    "bit_score_threshold", "e_value_threshold",
    "deepkoala_probability", "deepkoala_threshold", "eggnog_shared_seed_hit",
    "original_status", "cross_method_status", "consensus_level", "supporting_methods",
]

def _validate_gene_ids(values: pd.Series, path: Path, method: str) -> None:
    invalid = values.isna() | values.astype("string").str.strip().eq("")
    if invalid.any():
        rows = ", ".join(str(i + 2) for i in values.index[invalid][:5])
        raise ValueError(f"{method} file '{path}' contains missing or blank gene IDs at row(s): {rows}")

def parse_kos(val: Any) -> set[str]:
    """Extract KO tokens from common tool delimiters; missing values give no KOs."""
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
    """Sort KO identifiers for deterministic output."""
    return ",".join(sorted(kos)) if kos else "-"

def format_metric_value(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "-"
    return f"{numeric:.6g}" if np.isfinite(numeric) else "-"

def read_fasta_ids(fasta_path: Union[str, Path]) -> list[str]:
    """Read the first token of each header, preserving order and duplicates."""
    fasta_path = Path(fasta_path).expanduser().resolve()
    gene_ids = []

    with _open_text(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split()
                if parts:
                    gene_ids.append(parts[0])
    return gene_ids

def _open_text(path):
    if path.suffix.lower() in (".gz", ".gzip"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def load_ko_definitions(ko_list_path: Optional[Union[str, Path]]) -> dict[str, str]:
    """Optional master definitions; unavailable files preserve the historical fallback."""
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
    score_type: Optional[str], bit_score: Any, e_value: Any,
    domain_bit_score: Any, domain_e_value: Any,
) -> tuple[float, float, str]:
    """Domain profiles never fall back to full-sequence metrics."""
    kind = str(score_type).strip().lower() if pd.notna(score_type) else ""
    if kind == "full":
        score, evalue = bit_score, e_value
    elif kind == "domain":
        score, evalue = domain_bit_score, domain_e_value
    else:
        return np.nan, np.nan, "-"
    selected_score = float(score) if _finite(score) else np.nan
    selected_evalue = float(evalue) if _finite(evalue) else np.nan
    return selected_score, selected_evalue, kind

def make_evidence_record(
    gene_id: str, ko: str, method: str, original_status: str,
    bit_score: float = np.nan, e_value: float = np.nan,
    domain_bit_score: float = np.nan, domain_e_value: float = np.nan,
    score_type: str = "-", threshold: float = np.nan,
    bit_score_threshold: float = np.nan, e_value_threshold: float = np.nan,
    deepkoala_probability: float = np.nan, deepkoala_threshold: float = np.nan,
    eggnog_shared_seed_hit: bool = False, shared_seed_hit: bool = False,
    definition: str = "-", **tool_kwargs: Any,
) -> dict[str, Any]:
    """Normalize evidence defaults while preserving tool-specific metadata."""
    record = {
        "gene_id": gene_id, "ko": ko, "method": method, "original_status": original_status,
        "bit_score": bit_score, "e_value": e_value,
        "domain_bit_score": domain_bit_score, "domain_e_value": domain_e_value,
        "score_type": score_type or "-", "threshold": threshold,
        "bit_score_threshold": bit_score_threshold, "e_value_threshold": e_value_threshold,
        "deepkoala_probability": deepkoala_probability, "deepkoala_threshold": deepkoala_threshold,
        "eggnog_shared_seed_hit": eggnog_shared_seed_hit, "shared_seed_hit": shared_seed_hit,
        "definition": definition,
    }
    record.update(tool_kwargs)
    return record

def _numeric_columns(frame, columns):
    for column in columns:
        frame[column] = pd.to_numeric(frame.get(column), errors="coerce")


def _deduplicate(frame, status_column, score_columns, ascending):
    """Prefer status before metrics, breaking ties by original input order."""
    frame["_status_rank"] = frame[status_column].map({
        "threshold_passing": 3, "heuristic_rescued": 2,
        "below_threshold": 1, "unknown": 0,
    })
    frame["_orig_idx"] = np.arange(len(frame))
    frame = frame.sort_values(
        ["gene_id", "ko", "_status_rank"] + score_columns + ["_orig_idx"],
        ascending=[True, True, False] + ascending + [True], na_position="last",
    )
    return frame.drop_duplicates(["gene_id", "ko"], keep="first")


def _rename_columns(frame, aliases, strip_hash=False):
    """Use the first recognized alias, preserving each adapter's precedence."""
    names = []
    for column in frame.columns:
        name = column.lstrip("#").strip() if strip_hash else column
        names.append(name.lower())
    renamed = {}
    for target, candidates in aliases.items():
        for candidate in candidates:
            if candidate in names:
                renamed[frame.columns[names.index(candidate)]] = target
                break
    return frame.rename(columns=renamed)


def _read_tool_table(path):
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str)


def _finite(value):
    return pd.notna(value) and np.isfinite(value)


def extract_kofam_records(
    tsv_path: Union[str, Path], heuristic_bitscore_fraction: float = 0.75,
    heuristic_e_value: float = 1e-5,
) -> list[dict]:
    """Explicit assignment wins over inferred full/domain profile status."""
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []
    frame = _read_tool_table(path)
    missing = {"gene_id", "ko"} - set(frame.columns)
    if missing:
        raise ValueError(f"KOfam file '{path}' missing required column(s): {', '.join(sorted(missing))}")
    if frame.empty:
        return []
    _validate_gene_ids(frame["gene_id"], path, "KOfam")
    for column in ("gene_id", "ko"):
        frame[column] = frame[column].astype(str).str.strip()
    _numeric_columns(frame, ["bit_score", "e_value", "domain_bit_score", "domain_e_value", "threshold"])
    for column in ("score_type", "definition", "assignment"):
        if column not in frame.columns:
            frame[column] = "-"
    frame = frame[frame["ko"].str.fullmatch(KO_REGEX, na=False)].copy()
    if frame.empty:
        return []
    selected = []
    statuses = []
    assignments = {"threshold": "threshold_passing", "*": "threshold_passing",
                   "rescued": "heuristic_rescued", "below_threshold": "below_threshold",
                   "below": "below_threshold", "none": "below_threshold"}
    # Iterating columns avoids constructing a pandas Series for each hit.
    columns = zip(
        frame["score_type"], frame["bit_score"], frame["e_value"],
        frame["domain_bit_score"], frame["domain_e_value"], frame["threshold"],
        frame["assignment"],
    )
    for score_type, bit_score, e_value, domain_score, domain_evalue, threshold, assignment in columns:
        score, evalue, kind = select_kofam_scores(
            score_type, bit_score, e_value, domain_score, domain_evalue,
        )
        selected.append((score, evalue))
        assignment = str(assignment).strip().lower() if pd.notna(assignment) else ""
        status = assignments.get(assignment)
        if status is None:
            if kind != "-" and _finite(score) and _finite(threshold):
                status = "threshold_passing" if score >= threshold else "below_threshold"
            else:
                status = "unknown"
        statuses.append(status)
    frame["_sel_bit_score"] = [score for score, _ in selected]
    frame["_sel_e_value"] = [evalue for _, evalue in selected]
    frame["_original_status"] = statuses
    frame = _deduplicate(frame, "_original_status", ["_sel_bit_score", "_sel_e_value"], [False, True])
    records = []
    for _, row in frame.iterrows():
        threshold = row["threshold"]
        rescued = row["_original_status"] == "heuristic_rescued"
        cutoff = heuristic_bitscore_fraction * threshold if rescued and _finite(threshold) else threshold
        # A nonfinite rescued cutoff remains missing, as in the upstream adapter.
        if rescued and not _finite(threshold):
            cutoff = np.nan
        records.append(make_evidence_record(
            row["gene_id"], row["ko"], "kofam", row["_original_status"],
            bit_score=row["bit_score"], e_value=row["e_value"],
            domain_bit_score=row["domain_bit_score"], domain_e_value=row["domain_e_value"],
            score_type=str(row["score_type"]).strip() if pd.notna(row["score_type"]) else "-",
            threshold=threshold, bit_score_threshold=cutoff,
            e_value_threshold=float(heuristic_e_value) if rescued else np.nan,
            definition=str(row["definition"]).strip() if pd.notna(row["definition"]) else "-",
            assignment=(str(row["assignment"]).strip() if pd.notna(row["assignment"]) else "") or "-",
        ))
        if not records[-1]["definition"]:
            records[-1]["definition"] = "-"
    return records

def extract_deepkoala_records(tsv_path: Union[str, Path]) -> list[dict]:
    """Read upstream-filtered labels or predictions with threshold evidence."""
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []
    frame = _read_tool_table(path)
    columns = [column.lstrip("#").strip().lower() for column in frame.columns]
    simple = len(columns) == 2 and "name" in columns and "predict_label" in columns
    if not any(column in columns for column in ("name", "gene_id")):
        raise ValueError(f"DeepKOALA file '{path}' missing required 'name' or 'gene_id' column")
    if not any(column in columns for column in ("predict_label", "ko")):
        raise ValueError(f"DeepKOALA file '{path}' missing required 'predict_label' or 'ko' column")
    frame = _rename_columns(frame, {
        "gene_id": ("name", "gene_id"), "predict_label": ("predict_label", "ko"),
        "deepkoala_score": ("probability", "score"),
        "deepkoala_threshold": ("threshold",), "deepkoala_annotate": ("annotate",),
    }, strip_hash=True)
    if frame.empty:
        return []
    _validate_gene_ids(frame["gene_id"], path, "DeepKOALA")
    frame["gene_id"] = frame["gene_id"].astype(str).str.strip()
    _numeric_columns(frame, ["deepkoala_score", "deepkoala_threshold"])
    rows = []
    for _, row in frame.iterrows():
        score, threshold = row["deepkoala_score"], row["deepkoala_threshold"]
        annotation = row.get("deepkoala_annotate", "")
        annotation = str(annotation).strip() if pd.notna(annotation) else ""
        if simple:
            status, score_type = "threshold_passing", "upstream_filtered"
        elif "*" in annotation:
            status, score_type = "threshold_passing", "probability"
        elif _finite(score) and _finite(threshold):
            status = "threshold_passing" if score >= threshold else "below_threshold"
            score_type = "probability"
        else:
            status, score_type = "unknown", "-"
        for ko in parse_kos(row["predict_label"]):
            rows.append({
                "gene_id": row["gene_id"], "ko": ko, "original_status": status,
                "score_type": score_type, "deepkoala_score": score,
                "deepkoala_threshold": threshold,
            })
    if not rows:
        return []
    frame = _deduplicate(pd.DataFrame(rows), "original_status", ["deepkoala_score"], [False])
    records = []
    for _, row in frame.iterrows():
        records.append(make_evidence_record(
            row["gene_id"], row["ko"], "deepkoala", row["original_status"],
            score_type=row["score_type"], deepkoala_probability=row["deepkoala_score"],
            deepkoala_threshold=row["deepkoala_threshold"],
        ))
    return records

def read_eggnog_tsv(tsv_path: Union[str, Path]) -> pd.DataFrame:
    """Read native or normalized headers, skipping metadata, comments and blanks."""
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"eggNOG file not found: '{path}'")
    if path.stat().st_size == 0:
        raise ValueError(f"eggNOG annotations file is empty (0 bytes): '{path}'")

    header_line = None
    skip_indices: list[int] = []
    with _open_text(path) as f:
        for idx, line in enumerate(f):
            stripped = line.strip()
            if not stripped:
                skip_indices.append(idx)
                continue
            if stripped.startswith("##"):
                skip_indices.append(idx)
                continue
            if header_line is None:
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

    skip_rows = sorted(set(skip_indices) | {header_idx})
    with _open_text(path) as f:
        df = pd.read_csv(f, sep="\t", header=None, names=clean_cols, dtype=str, skiprows=skip_rows)

    if df.empty:
        return pd.DataFrame(columns=clean_cols)
    return df

def extract_eggnog_records(
    tsv_path: Union[str, Path], min_bitscore: float = 60.0, max_evalue: float = 1e-5,
) -> list[dict]:
    """Rank duplicate hits but retain all seed sources for group disambiguation."""
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        return []
    frame = read_eggnog_tsv(path)
    if frame.empty:
        return []
    frame = _rename_columns(frame, {
        "gene_id": ("query", "gene_id"), "eggnog_raw_ko": ("kegg_ko", "ko"),
        "eggnog_bit_score": ("score",), "eggnog_evalue": ("evalue",),
        "eggnog_description": ("description",),
    })
    _validate_gene_ids(frame["gene_id"], path, "eggNOG")
    frame["gene_id"] = frame["gene_id"].astype(str).str.strip()
    _numeric_columns(frame, ["eggnog_bit_score", "eggnog_evalue"])
    rows, sources_by_hit = [], {}
    for index, (_, row) in enumerate(frame.iterrows()):
        kos = parse_kos(row.get("eggnog_raw_ko", "-"))
        bit_score, e_value = row["eggnog_bit_score"], row["eggnog_evalue"]
        if not (_finite(bit_score) and _finite(e_value)):
            status = "unknown"
        elif bit_score >= min_bitscore and e_value <= max_evalue:
            status = "threshold_passing"
        else:
            status = "below_threshold"
        group = f"en_row{index}"
        description = row.get("eggnog_description", "-")
        description = str(description).strip() if pd.notna(description) else "-"
        for ko in kos:
            sources_by_hit.setdefault((row["gene_id"], ko), []).append({
                "group": group, "is_multi": len(kos) > 1, "status": status,
            })
            rows.append({
                "gene_id": row["gene_id"], "ko": ko, "original_status": status,
                "bit_score": bit_score, "e_value": e_value,
                "eggnog_seed_group": group, "eggnog_description": description,
            })
    if not rows:
        return []
    frame = _deduplicate(pd.DataFrame(rows), "original_status", ["bit_score", "e_value"], [False, True])
    records = []
    for _, row in frame.iterrows():
        sources = sources_by_hit[(row["gene_id"], row["ko"])]
        shared = any(source["is_multi"] for source in sources)
        records.append(make_evidence_record(
            row["gene_id"], row["ko"], "eggnog", row["original_status"],
            bit_score=row["bit_score"], e_value=row["e_value"], score_type="seed_hit",
            threshold=min_bitscore, bit_score_threshold=min_bitscore, e_value_threshold=max_evalue,
            eggnog_shared_seed_hit=shared, shared_seed_hit=shared,
            eggnog_seed_group=row["eggnog_seed_group"], eggnog_sources=sources,
            definition="-", eggnog_description=row["eggnog_description"],
        ))
    return records

def get_single_ko_definition(
    ko: str,
    ko_definitions: dict[str, str],
    evidence_definitions: Optional[dict[str, str]] = None,
) -> str:
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

def disambiguate_eggnog(
    en_agreed: set[str],
    tool_recs: dict[str, list[dict]],
    active_tools: list[str],
    confident_kos: dict[str, set[str]],
    candidate_kos: dict[str, set[str]],
    rescued_kos: dict[str, set[str]],
    eggnog_filter_multi: str = "disambiguate",
) -> tuple[set[str], set[str], set[Any], set[Any]]:
    """Narrow passing seed groups independently; singleton sources protect KOs."""
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

        grouped_kos = set()
        for kos in groups.values():
            grouped_kos.update(kos)
        retained = (en_agreed - grouped_kos) | singleton_support
        trusted_other, candidate_other = set(), set()
        for tool in active_tools:
            if tool != "eggnog":
                trusted_other.update(confident_kos.get(tool, set()))
                candidate_other.update(rescued_kos.get(tool, set()))
                candidate_other.update(candidate_kos.get(tool, set()))

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
        # Drop a KO globally only when none of its independent sources survives.
        en_dropped = set(confident_kos.get("eggnog", set())) - en_agreed

    return en_agreed, en_dropped, narrowed_groups, unresolved_groups

@dataclass(slots=True)
class GeneAnnotation:
    """One gene's evidence and the results filled in by the integration steps.

    Records are borrowed from the input adapters. Only this gene's indexes and
    decisions live here; output rows are emitted before processing the next gene.
    """
    gene_id: str
    records: list[dict]
    active_tools: list[str]
    tool_records: dict[str, list[dict]] = field(default_factory=dict)
    confident: dict[str, set[str]] = field(default_factory=dict)
    candidates: dict[str, set[str]] = field(default_factory=dict)
    rescued: dict[str, set[str]] = field(default_factory=dict)
    supporting: dict[str, set[str]] = field(default_factory=dict)
    eggnog_retained: set[str] = field(default_factory=set)
    eggnog_dropped: set[str] = field(default_factory=set)
    narrowed_groups: set = field(default_factory=set)
    unresolved_groups: set = field(default_factory=set)
    accepted_ko: str = "-"
    alternatives: set[str] = field(default_factory=set)
    consensus_level: str = "unannotated"
    evidence: str = "-"


def _categorize_gene(gene):
    gene.tool_records = {tool: [] for tool in gene.active_tools}
    categories = {
        "threshold_passing": gene.confident,
        "below_threshold": gene.candidates,
        "heuristic_rescued": gene.rescued,
    }
    for record in gene.records:
        tool = record["method"]
        if tool not in gene.tool_records:
            continue
        gene.tool_records[tool].append(record)
        category = categories.get(record["original_status"])
        if category is not None:
            category.setdefault(tool, set()).add(record["ko"])
            gene.supporting.setdefault(record["ko"], set()).add(tool)


def _disambiguate_gene(gene, strategy):
    retained, dropped, narrowed, unresolved = disambiguate_eggnog(
        set(gene.confident.get("eggnog", set())), gene.tool_records,
        gene.active_tools, gene.confident, gene.candidates, gene.rescued, strategy,
    )
    gene.eggnog_retained = retained
    gene.eggnog_dropped = dropped
    gene.narrowed_groups = narrowed
    gene.unresolved_groups = unresolved
    if retained:
        gene.confident["eggnog"] = retained


def _set_decision(gene, selected, alternatives, level, evidence):
    """Only a unique selected KO is accepted; ambiguity stays in alternatives."""
    gene.accepted_ko = next(iter(selected)) if len(selected) == 1 else "-"
    gene.alternatives = set(alternatives)
    if gene.accepted_ko == "-":
        gene.alternatives.update(selected)
    else:
        gene.alternatives.discard(gene.accepted_ko)
    gene.consensus_level = level
    gene.evidence = evidence


def _resolve_gene(gene, conflict_strategy="multiple", priority_order=None):
    calls = {tool: gene.confident[tool] for tool in gene.active_tools
             if tool in gene.confident}
    dropped = gene.eggnog_dropped
    if not calls:
        kofam = gene.rescued.get("kofam", set()) if "kofam" in gene.active_tools else set()
        deepkoala = gene.candidates.get("deepkoala", set()) if "deepkoala" in gene.active_tools else set()
        eggnog = gene.candidates.get("eggnog", set()) if "eggnog" in gene.active_tools else set()
        # Preserve the scientific policy: KOfam rescue takes precedence even
        # when DeepKOALA and eggNOG agree on a different candidate.
        if kofam:
            choices = [
                (kofam & deepkoala & eggnog, "deepkoala(candidate),eggnog(candidate),kofam(rescued)"),
                (kofam & eggnog, "eggnog(candidate),kofam(rescued)"),
                (kofam & deepkoala, "deepkoala(candidate),kofam(rescued)"),
            ]
            for selected, evidence in choices:
                if selected:
                    _set_decision(gene, selected, dropped, "dual_candidate", evidence)
                    return
            _set_decision(gene, kofam, dropped, "single_tool", "kofam(rescued)")
        elif deepkoala & eggnog:
            _set_decision(gene, deepkoala & eggnog, dropped, "dual_candidate",
                          "deepkoala(candidate),eggnog(candidate)")
        else:
            _set_decision(gene, set(), set(), "unannotated", "-")
        return

    if len(calls) == 1:
        tool, selected = next(iter(calls.items()))
        supporters = []
        for other, kos, label in [
            ("deepkoala", gene.candidates.get("deepkoala", set()), "deepkoala(candidate)"),
            ("eggnog", gene.candidates.get("eggnog", set()), "eggnog(candidate)"),
            ("kofam", gene.rescued.get("kofam", set()), "kofam(rescued)"),
        ]:
            if other != tool and other in gene.active_tools and selected & kos:
                supporters.append(label)
        level = "single_tool_with_candidate" if supporters else "single_tool"
        _set_decision(gene, selected, dropped, level, ",".join([tool] + supporters))
        return

    counts = Counter()
    for kos in calls.values():
        counts.update(kos)
    all_calls = set(counts)
    common = {ko for ko, count in counts.items() if count == len(calls)}
    selected = common or {ko for ko, count in counts.items() if count >= 2}
    if selected:
        level = "unanimous" if common and len(calls) == len(gene.active_tools) else "majority"
        supporters = sorted(tool for tool, kos in calls.items() if kos & selected)
        _set_decision(gene, selected, (all_calls - selected) | dropped,
                      level, ",".join(supporters))
    elif conflict_strategy == "priority":
        if priority_order is None:
            priority_order = [tool for tool in ("kofam", "eggnog", "deepkoala")
                              if tool in gene.active_tools]
        tool = next((tool for tool in priority_order if tool in calls), next(iter(calls)))
        _set_decision(gene, calls[tool], all_calls | dropped, "conflict_priority", tool)
    else:
        _set_decision(gene, set(), all_calls | dropped, "conflict", ",".join(sorted(calls)))

def adjudicate_gene(
    active_tools: list[str], confident_kos: dict[str, set[str]],
    candidate_kos: dict[str, set[str]], rescued_kos: dict[str, set[str]],
    en_dropped: set[str], conflict_strategy: str = "multiple",
    priority_order: Optional[list[str]] = None,
) -> tuple[str, set[str], str, str]:
    """Compatibility entry point for callers supplying categorized KO sets."""
    gene = GeneAnnotation("", [], active_tools)
    gene.confident = confident_kos
    gene.candidates = candidate_kos
    gene.rescued = rescued_kos
    gene.eggnog_dropped = en_dropped
    _resolve_gene(gene, conflict_strategy, priority_order)
    return gene.accepted_ko, gene.alternatives, gene.consensus_level, gene.evidence

def _summary_metric(records, key, numeric=True):
    """A scalar for one hit, otherwise KO-qualified values in record order."""
    if not records:
        return np.nan if numeric else "-"
    if len(records) == 1:
        return records[0][key]
    values = []
    for record in records:
        value = format_metric_value(record[key]) if numeric else str(record[key])
        values.append(f"{record['ko']}:{value}")
    return ",".join(values)


def _build_summary(gene, master_ko_defs, evidence_defs):
    accepted = {gene.accepted_ko} if gene.accepted_ko != "-" else set()
    row = {
        "gene_id": gene.gene_id,
        "accepted_ko": gene.accepted_ko,
        "alternative_kos": format_kos(gene.alternatives),
        "definition": resolve_definition_for_kos(accepted, master_ko_defs, evidence_defs),
        "alternative_definition": resolve_definition_for_kos(gene.alternatives, master_ko_defs, evidence_defs),
        "consensus_level": gene.consensus_level,
        "evidence": gene.evidence,
    }
    kofam = []
    for record in gene.tool_records.get("kofam", []):
        if record["original_status"] not in ("threshold_passing", "heuristic_rescued"):
            continue
        bit_score, e_value, score_type = select_kofam_scores(
            record.get("score_type"), record.get("bit_score"), record.get("e_value"),
            record.get("domain_bit_score"), record.get("domain_e_value"),
        )
        kofam.append({
            "ko": record["ko"], "score_type": score_type,
            "bit_score": bit_score, "e_value": e_value,
            "assignment": str(record["assignment"]), "threshold": record["threshold"],
        })
    row["kofam_ko"] = format_kos({record["ko"] for record in kofam})
    # score_type historically has no KO prefix and follows sorted KO order.
    types = {record["ko"]: record["score_type"] for record in kofam}
    row["kofam_score_type"] = ",".join(types[ko] for ko in sorted(types)) or "-"
    for column, key, numeric in [
        ("kofam_bit_score", "bit_score", True), ("kofam_evalue", "e_value", True),
        ("kofam_assignment", "assignment", False), ("kofam_threshold", "threshold", True),
    ]:
        row[column] = _summary_metric(kofam, key, numeric)
    for tool, metrics in [
        ("deepkoala", [("deepkoala_score", "deepkoala_probability"), ("deepkoala_threshold", "deepkoala_threshold")]),
        ("eggnog", [("eggnog_bit_score", "bit_score"), ("eggnog_evalue", "e_value")]),
    ]:
        records = sorted(
            (record for record in gene.tool_records.get(tool, [])
             if record["original_status"] in ("threshold_passing", "below_threshold")),
            key=lambda record: record["ko"],
        )
        passing = {record["ko"] for record in records if record["original_status"] == "threshold_passing"}
        if tool == "eggnog":
            passing = gene.eggnog_retained
        row[tool + "_ko"] = format_kos(passing)
        row[tool + "_candidate_ko"] = format_kos({record["ko"] for record in records})
        for column, key in metrics:
            row[column] = _summary_metric(records, key)
    return row

def build_gene_summary_row(
    gid: str, accepted_ko: str, alt_set: set[str], consensus_level: str,
    evidence_str: str, tool_recs: dict[str, list[dict]], en_agreed: set[str],
    master_ko_defs: dict[str, str], evidence_defs: dict[str, str],
) -> dict[str, Any]:
    """Compatibility entry point for an already adjudicated gene."""
    gene = GeneAnnotation(gid, [], [])
    gene.tool_records = tool_recs
    gene.accepted_ko = accepted_ko
    gene.alternatives = alt_set
    gene.consensus_level = consensus_level
    gene.evidence = evidence_str
    gene.eggnog_retained = en_agreed
    return _build_summary(gene, master_ko_defs, evidence_defs)

def _cross_status(gene, record, supporting, accepted):
    ko, method, status = record["ko"], record["method"], record["original_status"]
    if status == "unknown":
        return "unknown"
    if method == "eggnog" and ko in gene.eggnog_dropped:
        return "disambiguated_dropped"
    if ko in accepted:
        if status == "heuristic_rescued":
            return "heuristic_rescued"
        if status == "below_threshold":
            return "rescued"
        if status == "threshold_passing":
            if method == "eggnog" and _has_source(record, gene.narrowed_groups):
                return "disambiguated_retained"
            return "accepted"
        return "unknown"
    if ko in gene.alternatives:
        if status == "heuristic_rescued":
            return "heuristic_rescued"
        if status == "below_threshold":
            return "rescued" if supporting else "below_threshold_unrescued"
        if status == "threshold_passing":
            if gene.consensus_level.startswith("conflict"):
                return "conflict"
            if method == "eggnog" and _has_source(record, gene.unresolved_groups):
                return "unresolved_multi_ko"
            return "alternative"
        return "unknown"
    return "below_threshold_unrescued" if status == "below_threshold" else "unselected"


def _has_source(record, groups):
    for source in record.get("eggnog_sources", []):
        if source["group"] in groups and source["status"] == "threshold_passing":
            return True
    return False


def _build_audit(gene, accepted=None):
    if accepted is None:
        accepted = {gene.accepted_ko} if gene.accepted_ko != "-" else set()
    rows = []
    for record in gene.records:
        supporting = gene.supporting.get(record["ko"], set()) - {record["method"]}
        row = {column: record.get(column, np.nan) for column in EVIDENCE_COLUMNS}
        row.update({
            "gene_id": gene.gene_id, "consensus_level": gene.consensus_level,
            "score_type": record.get("score_type", "-"),
            "bit_score_threshold": record.get("bit_score_threshold", record.get("threshold", np.nan)),
            "eggnog_shared_seed_hit": bool(record.get("eggnog_shared_seed_hit", record.get("shared_seed_hit", False))),
            "cross_method_status": _cross_status(gene, record, supporting, accepted),
            "supporting_methods": ",".join(sorted(supporting)) or "-",
        })
        rows.append(row)
    return rows

def annotate_evidence_records(
    gid: str, recs: list[dict], tool_recs: dict[str, list[dict]],
    active_tools: list[str], accepted_set: set[str], alt_set: set[str],
    accepted_ko: str, consensus_level: str, en_dropped: set[str],
    narrowed_eggnog_groups: set[Any], unresolved_eggnog_groups: set[Any],
) -> list[dict]:
    """Compatibility entry point for callers with an existing decision."""
    gene = GeneAnnotation(gid, recs, active_tools)
    # Supporting records may differ from recs in this legacy helper interface.
    for tool in active_tools:
        for record in tool_recs.get(tool, []):
            if record["original_status"] in ("threshold_passing", "heuristic_rescued", "below_threshold"):
                gene.supporting.setdefault(record["ko"], set()).add(tool)
    gene.accepted_ko = accepted_ko
    gene.alternatives = alt_set
    gene.consensus_level = consensus_level
    gene.eggnog_dropped = en_dropped
    gene.narrowed_groups = narrowed_eggnog_groups
    gene.unresolved_groups = unresolved_eggnog_groups
    return _build_audit(gene, accepted_set)

def export_tables(
    result_df: pd.DataFrame, evidence_df: pd.DataFrame,
    output_tsv: Optional[Union[str, Path]] = None,
    evidence_tsv: Optional[Union[str, Path]] = None,
) -> None:
    """Write both tables, preserving scalar float and boolean formatting."""
    if output_tsv and evidence_tsv is None:
        evidence_tsv = Path(output_tsv).expanduser().resolve().parent / "kolach_evidence.tsv"
    for frame, destination in [(result_df, output_tsv), (evidence_df, evidence_tsv)]:
        if not destination:
            continue
        path = Path(destination).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        exported = frame.copy()
        for column in exported.columns:
            if pd.api.types.is_float_dtype(exported[column]):
                exported[column] = exported[column].apply(lambda value: "-" if pd.isna(value) else str(value))
            else:
                exported[column] = exported[column].fillna("-")
        exported.to_csv(path, sep="\t", index=False)

def integrate_annotations(
    protein_fasta: Optional[Union[str, Path]] = None,
    kofam_tsv: Optional[Union[str, Path]] = None,
    deepkoala_tsv: Optional[Union[str, Path]] = None,
    eggnog_tsv: Optional[Union[str, Path]] = None,
    output_tsv: Optional[Union[str, Path]] = None,
    evidence_tsv: Optional[Union[str, Path]] = None,
    database_dir: Optional[Union[str, Path]] = None,
    eggnog_min_bitscore: float = 60.0, eggnog_max_evalue: float = 1e-5,
    eggnog_filter_multi: str = "disambiguate", conflict_strategy: str = "multiple",
    heuristic_bitscore_fraction: float = 0.75, heuristic_e_value: float = 1e-5,
) -> pd.DataFrame:
    """Return the gene summary and its audit DataFrame in attrs["evidence"]."""
    if conflict_strategy not in {"multiple", "priority"}:
        raise ValueError("conflict_strategy must be 'multiple' or 'priority'")
    if eggnog_filter_multi not in {"disambiguate", "none"}:
        raise ValueError("eggnog_filter_multi must be 'disambiguate' or 'none'")
    parameters = {
        "eggnog_min_bitscore": eggnog_min_bitscore, "eggnog_max_evalue": eggnog_max_evalue,
        "heuristic_bitscore_fraction": heuristic_bitscore_fraction, "heuristic_e_value": heuristic_e_value,
    }
    for name, value in parameters.items():
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
    for name, path in [("protein_fasta", protein_fasta), ("kofam_tsv", kofam_tsv),
                       ("deepkoala_tsv", deepkoala_tsv), ("eggnog_tsv", eggnog_tsv)]:
        if path is not None:
            resolved = Path(path).expanduser().resolve()
            if not resolved.is_file():
                raise FileNotFoundError(
                    f"Supplied {name} file does not exist or is not a regular file: '{path}' (resolved: '{resolved}')"
                )
    active_tools, records, evidence_defs = [], [], {}
    if kofam_tsv is not None:
        records.extend(extract_kofam_records(kofam_tsv, heuristic_bitscore_fraction, heuristic_e_value))
        active_tools.append("kofam")
        for record in records:
            if record.get("definition") and record["definition"] != "-":
                evidence_defs[record["ko"]] = record["definition"]
    if deepkoala_tsv is not None:
        records.extend(extract_deepkoala_records(deepkoala_tsv))
        active_tools.append("deepkoala")
    if eggnog_tsv is not None:
        records.extend(extract_eggnog_records(eggnog_tsv, eggnog_min_bitscore, eggnog_max_evalue))
        active_tools.append("eggnog")
    gene_ids = read_fasta_ids(protein_fasta) if protein_fasta is not None else []
    records_by_gene = {gid: [] for gid in gene_ids}
    for record in records:
        gid = record["gene_id"]
        if gid not in records_by_gene:
            records_by_gene[gid] = []
            gene_ids.append(gid)
        records_by_gene[gid].append(record)
    ko_list = None
    if database_dir:
        directory = Path(database_dir).expanduser().resolve()
        for candidate in [directory / "kofam" / "ko_list", directory / "ko_list"]:
            if candidate.is_file():
                ko_list = candidate
                break
    master_defs = load_ko_definitions(ko_list)
    summaries, audit = [], []
    for gid in gene_ids:
        gene = GeneAnnotation(gid, records_by_gene[gid], active_tools)
        _categorize_gene(gene)
        _disambiguate_gene(gene, eggnog_filter_multi)
        _resolve_gene(gene, conflict_strategy)
        summaries.append(_build_summary(gene, master_defs, evidence_defs))
        audit.extend(_build_audit(gene))
        del gene
    evidence_df = pd.DataFrame(audit, columns=EVIDENCE_COLUMNS)
    if audit:
        evidence_df = evidence_df.sort_values(["gene_id", "ko", "method"]).reset_index(drop=True)
    columns = [
        "gene_id", "accepted_ko", "alternative_kos", "definition", "alternative_definition",
        "consensus_level", "evidence", "kofam_ko", "kofam_score_type",
        "kofam_bit_score", "kofam_evalue", "kofam_assignment", "kofam_threshold",
        "deepkoala_ko", "deepkoala_candidate_ko", "deepkoala_score", "deepkoala_threshold",
        "eggnog_ko", "eggnog_candidate_ko", "eggnog_bit_score", "eggnog_evalue",
    ]
    result = pd.DataFrame(summaries, columns=columns)
    result.attrs["evidence"] = evidence_df
    export_tables(result, evidence_df, output_tsv, evidence_tsv)
    return result
