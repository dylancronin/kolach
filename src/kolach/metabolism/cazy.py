"""Genome-level CAZy substrate calls from Kolach's dbCAN family/EC evidence.

Substrate definitions are DRAM2's current CAZy rules (dram-rules, pinned in
RULES_SOURCE), not AnnoGuild's legacy `DRAM_CAZy_definitions.tsv` (which DRAM2
keeps only as its opt-in "legacy" sheet). Rules come from the common alias
sheet and from the ecosystem sheets that add CAZy substrates of their own (gut:
host glycans and beta-glucans; aquatic: marine algal polysaccharides; eng_sys:
galactomannan). Each rule used is an OR list of `FAMILY;EC:x.x.x.x` pairs; rules
mixing KOs are not evaluated. A genome has a substrate when any of its proteins
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

# Vendored unchanged from dram-rules src/dram_rules/data/ (MIT).
RULE_SHEETS = {
    "common": DATA_DIR / "dram_rules" / "rules_common.tsv",
    "gut": DATA_DIR / "dram_rules" / "rules_gut.tsv",
    "aquatic": DATA_DIR / "dram_rules" / "rules_aquatic.tsv",
    "eng_sys": DATA_DIR / "dram_rules" / "rules_eng_sys.tsv",
}
RULES_SOURCE = {
    "repository": "https://github.com/BortonWrightonLabs/dram-rules",
    "commit": "3257b8604ab6e40c00e811d6ba0211b44dcf5648",
    "license": "MIT",
    "files": {
        "src/dram_rules/data/rules_common.tsv": "9b0530148470801987460e46d6508ed5197c4d0177d334584eaedef50a626c41",
        "src/dram_rules/data/rules_gut.tsv": "96fff41adebabcfa2ce9e1eeaed29b8e8b8ddf270041c76822ca127cf51d24af",
        "src/dram_rules/data/rules_aquatic.tsv": "cbeebf9852e9e8b2d763d6a3cb3c3791be79e24a0ab985a07735c252fd70ffd7",
        "src/dram_rules/data/rules_eng_sys.tsv": "bc18f6e28c4f586b713f3cc82d34b4e41d0127c36c2bede8edde178a3a2d1f6c",
    },
}

# Product column -> DRAM2 rule aliases. Subcategories of one polymer are merged
# so they cannot be counted twice toward the generalist threshold.
# Provisional: chosen for this release and expected to be revisited.
# DRAM2's "Microbially-derived" (peptidoglycan) rule is deliberately not a
# group: it includes GH23/GH102-104 lytic transglycosylases that nearly every
# bacterium uses for its own cell wall, which would make almost every genome a
# macromolecule degrader. It is still reported as CAZy_other-Microbially-derived.
# Ecosystem rules that DRAM2 classes as oligomers, or that act on oligosaccharides,
# glycolipids or terminal residues (milk glycans, HMOs, sialic-acid release,
# host glycosphingolipids, alginate/agarose "Associated"), are reported only.
POLYMER_GROUPS = {
    "CAZy-Cellulose": ("Cellulose",),
    "CAZy-Hemicellulose": ("Hemicellulose-Xylan", "Hemicellulose-Other", "Galactomannan (e.g., Guar-Gum)"),
    "CAZy-Pectin": ("Pectin",),
    "CAZy-Arabinogalactan": ("Arabinogalactan",),
    "CAZy-Cutin": ("Cutin",),
    "CAZy-Fructan": ("Fructan",),
    "CAZy-Starch_glycogen": ("StarchSpecific", "StarchGlycogen"),
    "CAZy-Chitin": ("Chitin",),
    "CAZy-Beta_glucan": ("Beta-Glucan Other", "Laminarin"),
    "CAZy-Mucin": ("Mucin-O-glycan", "Blood-group-glycan"),
    "CAZy-Glycosaminoglycans": ("Hyaluronan", "Chondroitin-sulfate", "Dermatan-sulfate", "Keratan-sulfate",
                                "Heparin-heparan-sulfate", "Mixed-heparan-sulfate/dermatan-sulfate"),
    "CAZy-Host_N_glycan": ("Host-N-glycan",),
    "CAZy-Alginate": ("Alginate",),
    "CAZy-Fucoidan": ("Fucoidan",),
    "CAZy-Ulvan": ("Ulvan",),
    "CAZy-Red_algal_galactans": ("Agarose", "Carrageenan", "Porphyran", "Funoran"),
}
# Families AnnoGuild judged to match poorly; used only by the signature-only
# comparison column, never by the guild calls.
LOW_MATCHING_FAMILIES = {"GH63", "GH49", "GH127", "CE2", "GH23", "PL30", "PL8", "PL33", "GH88", "GH67", "GH98", "CE4"}
# Matched after upper-casing (as DRAM2 does), so provisional ECs read 2.4.1.N2.
PAIR = re.compile(r"^([A-Z]+\d+);EC:(\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:N?\d+|-))$")
DETAIL_COLUMNS = ["genome", "substrate", "dram2_rules", "counts_toward_guild", "called",
                  "called_signature_only", "matched_pairs", "genes"]


def normalise(pair: str) -> str:
    """Upper-case, and add the EC: prefix some sheets omit (GH16;3.2.1.81),
    as DRAM2's matcher treats both spellings alike."""
    pair = pair.strip().upper()
    family, separator, ec = pair.partition(";")
    if separator and not ec.startswith("EC:"):
        pair = f"{family};EC:{ec}"
    return pair


def load_cazy_rules(sheets=None) -> dict[str, frozenset]:
    """Return every DRAM2 rule made only of family/EC pairs: name -> pairs.

    The common sheet names rules by `alias`; ecosystem sheets by `name`.
    """
    rules = {}
    for path in (sheets or RULE_SHEETS).values():
        with open(path, newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        header = rows[0]
        if "rule" not in header or "name" not in header:
            raise ValueError(f"Unexpected DRAM2 rule header in {path}: {header}")
        name_column = header.index("alias" if "alias" in header else "name")
        rule_column = header.index("rule")
        for row in rows[1:]:
            if len(row) <= max(name_column, rule_column) or ";" not in row[rule_column]:
                continue
            name = row[name_column].strip()
            terms = [term.strip() for term in row[rule_column].split("|")]
            if name and all(PAIR.match(normalise(term)) for term in terms):
                if name in rules:
                    raise ValueError(f"Duplicate DRAM2 CAZy rule {name!r} ({path})")
                rules[name] = frozenset(normalise(term) for term in terms)
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
    """Pairs unique to one substrate (over all family/EC rules), as in AnnoGuild,
    excluding its low-matching families."""
    owners = {}
    for group, aliases in substrate_groups(rules).items():
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
    product columns {genome: {column: bool}}, the same columns using signature
    pairs only (for comparison), and detail rows for every substrate.
    """
    rules = rules or load_cazy_rules()
    groups = substrate_groups(rules)
    signatures = signature_pairs(rules)
    products, signature_products, details = {}, {}, []
    for genome in sorted(genome_gene_pairs):
        genes = genome_gene_pairs[genome]
        products[genome], signature_products[genome] = {}, {}
        for group, aliases in groups.items():
            wanted = set().union(*(rules[alias] for alias in aliases))
            matched = {gene: genes[gene] & wanted for gene in sorted(genes) if genes[gene] & wanted}
            pairs = set().union(*matched.values()) if matched else set()
            called = bool(pairs)
            polymer = group in POLYMER_GROUPS
            if polymer:
                products[genome][group] = called
                signature_products[genome][group] = bool(pairs & signatures)
            details.append({
                "genome": genome, "substrate": group, "dram2_rules": ";".join(aliases),
                "counts_toward_guild": polymer, "called": called,
                "called_signature_only": bool(pairs & signatures) if polymer else "",
                "matched_pairs": ";".join(sorted(pairs)), "genes": ";".join(matched),
            })
    return products, signature_products, details


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
