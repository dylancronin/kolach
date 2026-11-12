import json
import re
from pathlib import Path
from typing import Union
from urllib.request import Request
import pandas as pd

from kolach import package_version
from kolach.download import download_atomic


KEGG_HIERARCHY_URL = (
    "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir="
)


KO_REGEX = re.compile(r"^K\d{5}$")
KO_TOKEN_REGEX = re.compile(r"\bK\d{5}\b")


def is_valid_hierarchy_json(file_path: Path) -> bool:
    """Validate that the cached file exists, is non-empty, and parses as a valid KEGG hierarchy."""
    if not file_path.is_file() or file_path.stat().st_size == 0:
        return False
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return isinstance(data, dict) and "children" in data
    except Exception:
        return False


def download_kegg_hierarchy(
    database_dir: Union[str, Path],
    force_download: bool = False,
    verbose: bool = True,
) -> Path:
    """Download master KEGG Orthology hierarchy (ko00001.json) into database_dir.

    Args:
        database_dir: Directory where kolach databases are stored.
        force_download: If True, re-download even if a valid cached file exists.
        verbose: If True, print progress messages.

    Returns:
        Path: Resolved absolute path to the downloaded ko00001.json file.
    """
    db_path = Path(database_dir).expanduser().resolve()
    db_path.mkdir(parents=True, exist_ok=True)
    dest_path = db_path / "ko00001.json"

    # If already downloaded and valid, skip re-downloading unless forced
    if not force_download and is_valid_hierarchy_json(dest_path):
        return dest_path

    if verbose:
        print(f"[kolach] Downloading KEGG Orthology hierarchy (ko00001.json) to {dest_path}...")

    def validate(path: Path) -> None:
        if not is_valid_hierarchy_json(path):
            raise ValueError(
                "Downloaded file is not a valid KEGG hierarchy JSON (check network/permissions)."
            )

    req = Request(KEGG_HIERARCHY_URL, headers={"User-Agent": f"kolach/{package_version()}"})
    try:
        download_atomic(req, dest_path, validate=validate)
    except Exception as e:
        raise RuntimeError(
            f"Failed to download KEGG hierarchy from {KEGG_HIERARCHY_URL}: {e}"
        ) from e
    if verbose:
        print(f"[kolach] Downloaded ko00001.json ({dest_path.stat().st_size / 1024 / 1024:.1f} MB).")

    return dest_path


def _collect_kegg_hierarchy(json_path: Union[str, Path]) -> dict[str, dict[str, set[str]]]:
    """Map each KO to the sets of level A/B/C names it appears under in ko00001.json."""
    json_path = Path(json_path).expanduser().resolve()
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # Sets deduplicate KOs appearing in multiple places within the hierarchy
    raw_mapping: dict[str, dict[str, set[str]]] = {}

    # Tier 1 (Level A): Top-level categories (e.g., "09100 Metabolism")
    for cat_a_node in data.get("children", []):
        cat_a_name = cat_a_node.get("name", "").strip()

        # Tier 2 (Level B): Subcategories (e.g., "09101 Carbohydrate metabolism")
        for cat_b_node in cat_a_node.get("children", []):
            cat_b_name = cat_b_node.get("name", "").strip()

            # Tier 3 (Level C): Pathways / Brite groups (e.g., "00010 Glycolysis ... [PATH:ko00010]")
            for cat_c_node in cat_b_node.get("children", []):
                cat_c_name = cat_c_node.get("name", "").strip()
                # Remove brackets suffix like [PATH:ko00010] or [BR:ko01000]
                clean_pathway = re.sub(r"\s*\[(PATH|BR):[^\]]+\]", "", cat_c_name).strip()

                # Tier 4 (Level D): KO nodes (e.g., "K00844  HK, glk; hexokinase [EC:2.7.1.1]")
                for ko_node in cat_c_node.get("children", []):
                    parts = ko_node.get("name", "").strip().split(None, 1)
                    if not (parts and KO_REGEX.match(parts[0])):
                        continue
                    info = raw_mapping.setdefault(parts[0], {
                        "category": set(), "subcategory": set(), "pathway": set(),
                    })
                    if cat_a_name:
                        info["category"].add(cat_a_name)
                    if cat_b_name:
                        info["subcategory"].add(cat_b_name)
                    if clean_pathway:
                        info["pathway"].add(clean_pathway)

    return raw_mapping


def _join(values: set[str]) -> str:
    """Sorted, semicolon-separated values, or '-' when empty."""
    return "; ".join(sorted(values)) if values else "-"


def parse_kegg_hierarchy(json_path: Union[str, Path]) -> dict[str, dict[str, str]]:
    """Parse ko00001.json into a mapping of KO to functional categories and pathways.

    Preserves official KEGG numerical identifiers on categories and subcategories
    to ensure canonical sorting.

    Args:
        json_path: Path to the ko00001.json file.

    Returns:
        dict[str, dict[str, str]]: Dictionary mapping ko_id to:
            {
                "category": "09100 Metabolism",
                "subcategory": "09101 Carbohydrate metabolism",
                "pathway": "00010 Glycolysis / Gluconeogenesis; ..."
            }
    """
    return {
        ko_id: {level: _join(names) for level, names in info.items()}
        for ko_id, info in _collect_kegg_hierarchy(json_path).items()
    }


def annotate_pathways(
    annotation_table: Union[str, Path],
    output_tsv: Union[str, Path],
    database_dir: Union[str, Path],
    ko_col: str = "accepted_ko",
    verbose: bool = True,
) -> pd.DataFrame:
    """Annotate a table containing KO identifiers with KEGG functional categories and pathways.

    Args:
        annotation_table: Path to the input tabular annotations file (TSV).
        output_tsv: Destination file path for the annotated output TSV.
        database_dir: Directory where kolach databases are stored (contains or will store ko00001.json).
        ko_col: Name of the column containing KO identifiers (default: 'accepted_ko').
        verbose: If True, print progress and status messages.

    Returns:
        pd.DataFrame: Augmented DataFrame containing category, subcategory, and pathway columns.
    """
    in_path = Path(annotation_table).expanduser().resolve()
    if not in_path.is_file():
        raise FileNotFoundError(f"Annotation table does not exist: {annotation_table} (resolved: {in_path})")

    db_path = Path(database_dir).expanduser().resolve()
    # Resolve existing ko00001.json or download it directly to database_dir
    if is_valid_hierarchy_json(db_path / "ko00001.json"):
        json_path = db_path / "ko00001.json"
    elif is_valid_hierarchy_json(db_path / "kofam" / "ko00001.json"):
        json_path = db_path / "kofam" / "ko00001.json"
    else:
        json_path = download_kegg_hierarchy(db_path, verbose=verbose)

    if verbose:
        print(f"[kolach] Parsing KEGG Orthology hierarchy from {json_path}...")
    kegg_map = _collect_kegg_hierarchy(json_path)
    if verbose:
        print(f"[kolach] Loaded {len(kegg_map):,} KO hierarchy mappings.")

    # keep_default_na=False keeps literal cell values such as "NA" or "None" (for
    # example gene IDs) intact instead of turning them into blanks on write-back.
    df = pd.read_csv(in_path, sep="\t", dtype=str, keep_default_na=False)
    if ko_col not in df.columns:
        raise ValueError(
            f"Column '{ko_col}' not found in {in_path}. Available columns: {list(df.columns)}"
        )

    if verbose:
        print(f"[kolach] Annotating '{ko_col}' across {len(df):,} genes in {in_path.name}...")

    # Function to lookup a single or multi-KO string
    def lookup_kos(val: str) -> tuple[str, str, str]:
        if pd.isna(val) or not val or str(val).strip() in ("", "-", "None", "nan"):
            return "-", "-", "-"
        # Parse all KO identifiers matching standard K\d{5} regex (same as integrate.py)
        kos = KO_TOKEN_REGEX.findall(str(val))
        if not kos:
            return "-", "-", "-"

        cats, subcats, paths = set(), set(), set()
        for k in kos:
            info = kegg_map.get(k)
            if info is not None:
                cats |= info["category"]
                subcats |= info["subcategory"]
                paths |= info["pathway"]
        return _join(cats), _join(subcats), _join(paths)

    results = df[ko_col].apply(lookup_kos)
    cat_series = results.apply(lambda x: x[0])
    subcat_series = results.apply(lambda x: x[1])
    path_series = results.apply(lambda x: x[2])

    # If columns already exist (e.g. re-running or updating), remove them before inserting
    existing_to_drop = [
        c for c in [
            "kegg_brite_category", "kegg_brite_subcategory", "kegg_brite_pathway",
            "category", "subcategory", "pathway",
        ] if c in df.columns
    ]
    if existing_to_drop:
        df = df.drop(columns=existing_to_drop)

    # Insert columns directly adjacent to definition/alternative_definition or ko_col
    insert_after = None
    for candidate in ["alternative_definition", "definition", ko_col]:
        if candidate in df.columns:
            insert_after = candidate
            break

    if insert_after is not None:
        idx = df.columns.get_loc(insert_after) + 1
        df.insert(idx, "kegg_brite_pathway", path_series)
        df.insert(idx, "kegg_brite_subcategory", subcat_series)
        df.insert(idx, "kegg_brite_category", cat_series)
    else:
        df["kegg_brite_category"] = cat_series
        df["kegg_brite_subcategory"] = subcat_series
        df["kegg_brite_pathway"] = path_series

    out_path = Path(output_tsv).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

    if verbose:
        annotated_count = (cat_series != "-").sum()
        print(
            f"[kolach] Successfully annotated {annotated_count:,}/{len(df):,} genes with KEGG BRITE hierarchies.\n"
            f"[kolach] Wrote {out_path} ({out_path.stat().st_size / 1024 / 1024:.2f} MB)."
        )

    return df
