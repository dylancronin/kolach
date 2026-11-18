"""Genome-level CAZy substrate calls from Kolach's dbCAN family/EC evidence.

Substrate definitions are DRAM2's current CAZy rules (dram-rules
`rules_common.tsv`, pinned in `data/dram_rules/manifest.json`), not AnnoGuild's
legacy `DRAM_CAZy_definitions.tsv`. Each DRAM2 substrate rule is an OR list of
`FAMILY;EC:x.x.x.x` pairs. A genome has a substrate when any of its proteins
carries one of those pairs in the `dbcan_family_ec` annotation column.

Differences from DRAM2's own evaluator, by design:
- Pairs come from the same dbCAN-sub domain (Kolach's `dbcan_family_ec`). DRAM2
  pairs every family on a protein with every EC on it.
- Subfamilies are reduced to families (GH5_2 -> GH5), as in Kolach's dbCAN
  parser, so rules written at family level can match subfamily hits.
Matching is otherwise exact, as in DRAM2: case-insensitive, and a partial EC
such as 3.2.1.- matches only the same partial EC.

Only substrate groups that are polymers count toward AnnoGuild's generalist and
macromolecule-degrader labels. The grouping is provisional (see POLYMER_GROUPS).
"""
import csv
from pathlib import Path
import re

from .evaluate import DATA_DIR

RULES_FILE = DATA_DIR / "dram_rules" / "rules_common.tsv"

# Product column -> DRAM2 rule aliases. Subcategories of one polymer are merged
# so they cannot be counted twice toward the generalist threshold.
# Provisional: chosen for this release and expected to be revisited.
POLYMER_GROUPS = {
    "CAZy-Cellulose": ("Cellulose",),
    "CAZy-Hemicellulose": ("Hemicellulose-Xylan", "Hemicellulose-Other"),
    "CAZy-Pectin": ("Pectin",),
    "CAZy-Arabinogalactan": ("Arabinogalactan",),
    "CAZy-Cutin": ("Cutin",),
    "CAZy-Fructan": ("Fructan",),
    "CAZy-Starch_glycogen": ("StarchSpecific", "StarchGlycogen"),
    "CAZy-Chitin": ("Chitin",),
    "CAZy-Microbial_cell_wall": ("Microbially-derived",),
}
# Families AnnoGuild judged to match poorly; used only by the signature-only
# comparison column, never by the guild calls.
LOW_MATCHING_FAMILIES = {"GH63", "GH49", "GH127", "CE2", "GH23", "PL30", "PL8", "PL33", "GH88", "GH67", "GH98", "CE4"}
# Matched after upper-casing (as DRAM2 does), so provisional ECs read 2.4.1.N2.
PAIR = re.compile(r"^([A-Z]+\d+);EC:(\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:N?\d+|-))$")
DETAIL_COLUMNS = ["genome", "substrate", "dram2_rules", "counts_toward_guild", "called",
                  "called_signature_only", "matched_pairs", "genes"]


def normalise(pair: str) -> str:
    return pair.strip().upper()


def load_cazy_rules(path=RULES_FILE) -> dict[str, frozenset]:
    """Return every DRAM2 rule made only of family/EC pairs: alias -> pairs."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0][:3]
    if header != ["name", "alias", "rule"]:
        raise ValueError(f"Unexpected DRAM2 rule header in {path}: {rows[0]}")
    rules = {}
    for row in rows[1:]:
        if len(row) < 3 or ";EC:" not in row[2]:
            continue
        alias = row[1].strip()
        terms = [term.strip() for term in row[2].split("|")]
        if all(PAIR.match(normalise(term)) for term in terms):
            if alias in rules:
                raise ValueError(f"Duplicate DRAM2 CAZy rule {alias!r}")
            rules[alias] = frozenset(normalise(term) for term in terms)
    for aliases in POLYMER_GROUPS.values():
        for alias in aliases:
            if alias not in rules:
                raise ValueError(f"DRAM2 rule {alias!r} is missing or not a family/EC list")
    return rules


def substrate_groups(rules) -> dict[str, tuple]:
    """Polymer groups first, then every other family/EC rule on its own."""
    grouped = {alias for aliases in POLYMER_GROUPS.values() for alias in aliases}
    groups = dict(POLYMER_GROUPS)
    for alias in rules:
        if alias not in grouped:
            groups[f"CAZy_other-{alias}"] = (alias,)
    return groups


def signature_pairs(rules) -> set[str]:
    """Pairs that belong to exactly one polymer group and a well-matching family."""
    owners = {}
    for group, aliases in POLYMER_GROUPS.items():
        for alias in aliases:
            for pair in rules[alias]:
                owners.setdefault(pair, set()).add(group)
    return {pair for pair, groups in owners.items()
            if len(groups) == 1 and pair.split(";")[0] not in LOW_MATCHING_FAMILIES}


def gene_pairs(value: str) -> set[str]:
    """Parse a `dbcan_family_ec` cell such as 'GH13;EC:3.2.1.1|GH5;EC:3.2.1.4'."""
    if value.strip() in ("", "-"):
        return set()
    pairs = {normalise(pair) for pair in value.split("|")}
    bad = [pair for pair in pairs if not PAIR.match(pair)]
    if bad:
        raise ValueError(f"Malformed dbcan_family_ec value(s): {bad}")
    return pairs


def call_substrates(genome_gene_pairs: dict, rules=None):
    """Evaluate substrates per genome.

    `genome_gene_pairs` is {genome: {gene: set(pairs)}}. Returns the polymer
    product columns {genome: {column: bool}} and detail rows for every substrate.
    """
    rules = rules or load_cazy_rules()
    groups = substrate_groups(rules)
    signatures = signature_pairs(rules)
    products, details = {}, []
    for genome in sorted(genome_gene_pairs):
        genes = genome_gene_pairs[genome]
        products[genome] = {}
        for group, aliases in groups.items():
            wanted = set().union(*(rules[alias] for alias in aliases))
            matched = {gene: genes[gene] & wanted for gene in sorted(genes) if genes[gene] & wanted}
            pairs = set().union(*matched.values()) if matched else set()
            called = bool(pairs)
            polymer = group in POLYMER_GROUPS
            if polymer:
                products[genome][group] = called
            details.append({
                "genome": genome, "substrate": group, "dram2_rules": ";".join(aliases),
                "counts_toward_guild": polymer, "called": called,
                "called_signature_only": bool(pairs & signatures) if polymer else "",
                "matched_pairs": ";".join(sorted(pairs)), "genes": ";".join(matched),
            })
    return products, details


def read_gene_pairs(annotation_table, gene_genomes: dict, genomes=()) -> dict:
    """Collect dbCAN family/EC pairs per genome from a Kolach annotation table.

    `genomes` lists every assessed genome, including any without genes.
    """
    with open(annotation_table, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "dbcan_family_ec" not in (reader.fieldnames or []):
            raise ValueError("Annotation table has no dbcan_family_ec column; run annotation with dbCAN")
        result = {genome: {} for genome in {*genomes, *gene_genomes.values()}}
        for row in reader:
            genome = gene_genomes.get(row["gene_id"])
            if genome is None:
                raise ValueError(f"Gene {row['gene_id']!r} has no genome")
            result[genome][row["gene_id"]] = gene_pairs(row["dbcan_family_ec"])
    return result


def write_details(path, details):
    with open(Path(path), "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DETAIL_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(details)
