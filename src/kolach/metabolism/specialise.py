"""Aggregate accepted KOs, evaluate pathways, and write legacy-format guilds.

The four-column specialisation table and wide Boolean product match AnnoGuild.
Unassessed evidence and pathway explanations live in separate reports so existing
analysis scripts can continue to read the original schemas.
"""

import csv
import hashlib
import json
import math
from pathlib import Path
import re

from .evaluate import DATA_DIR, evaluate_pathway, load_rules, read_tsv


# KO signatures retain the old thresholds; only the missing comma is corrected.
KO_MCR = {"K00399", "K00401", "K00402"}
KO_HDR = {"K22480", "K22481", "K22482", "K03388", "K03389", "K03390", "K08264", "K08265"}
KO_FWD_FMD = {"K00200", "K00201", "K00202", "K00203", "K00204", "K00205", "K11260", "K11261"}
KO_RAMX = {"K16183"}
KO_MTXBCA = {"K14080", "K14081", "K04480", "K14082", "K16176", "K16177",
             "K16178", "K16179", "K14083", "K14084", "K16954"}
KO_ACETO = {"K01895", "K00192", "K00193", "K00194", "K00195", "K00197"}
KO_METHANOTROPHY = {"K10944", "K10945", "K10946", "K109440", "K109450", "K109460",
                    "K16157", "K16158", "K16159", "K16160", "K16161", "K16162"}

# Custom IDs are AnnoGuild's refined features, not official six-digit KOs.
# Keep them internally for rule compatibility, while allowing named marker input.
MARKER_IDS = {
    "dsrA": "K111800", "rdsrA": "K111801", "dsrB": "K111810", "rdsrB": "K111811",
    "narG": "K003700", "nxrA": "K003701", "narH": "K003710", "nxrB": "K003711",
    "pmoA": "K109440", "amoA": "K109441", "pmoB": "K109450", "amoB": "K109451",
    "pmoC": "K109460", "amoC": "K109461", "norB": "K045610", "norZ": "K045611",
    "mmoX": "K16157",
}
NITROGEN_ROLES = [
    ("nitrogen_redox-comammox", "comammox"),
    ("nitrogen_redox-ammonia_oxidation", "ammonia_oxidiser"),
    ("nitrogen_redox-denitrification", "denitrifier"),
    ("nitrogen_redox-nitrate_reduction", "nitrate_reducer"),
    ("nitrogen_redox-nitrite_oxidation", "nitrite_oxidiser"),
    ("nitrogen_redox-nitrogen_fixation", "nitrogen_fixer"),
]
SULFUR_ROLES = [
    ("sulfur_redox-dissimilatory_sulfate_reduction", "sulfate_reducer"),
    ("sulfur_redox-dissimilatory_sulfur_oxidation", "sulfur_oxidiser"),
]
SPECIALISATION_COLUMNS = ["genome", "specialisation", "nitrogen_specialisation", "sulfur_specialisation"]
DETAIL_COLUMNS = ["genome", "pathway", "called", "reaction_coverage", "ko_coverage",
                  "signature_present", "ko_coverage_route", "supporting_features", "missing_features",
                  "missing_signature_features", "required_refinements", "assessment"]


def ko_predict_methanogen(kos: set[str]) -> bool:
    """Apply the original core-plus-route heuristic, including its exclusions."""
    if kos & KO_METHANOTROPHY or not kos & (KO_HDR | KO_MCR):
        return False
    return (len(kos & KO_FWD_FMD) >= 4
            or bool(kos & KO_RAMX) and bool(kos & KO_MTXBCA)
            or len(kos & KO_ACETO) >= 4)


def has_methanotrophy_marker(kos: set[str]) -> bool:
    """Original 'ko_fdh' check, renamed to describe what it actually tests."""
    return bool(kos & KO_METHANOTROPHY)


def parse_bool(value) -> bool:
    """Accept actual booleans and explicit TSV booleans; never bool('False')."""
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1"}:
        return True
    if text in {"false", "0", "", "none", "nan"}:
        return False
    raise ValueError(f"Expected a Boolean pathway call, got {value!r}")


def classify_specialisation(pathways: dict, kos: set[str]) -> dict:
    """Return the original three labels, preserving their priority ordering."""
    called = {pathway for pathway, value in pathways.items() if parse_bool(value)}
    methanogenesis = {
        "hydrogenotrophic_methanogenesis-all", "acetoclastic_methanogenesis-all",
        "methylotrophic_methanogenesis-h2_dependent", "methylotrophic_methanogenesis-h2_independent",
        # Also accept old shortened labels in externally supplied product matrices.
        "methylotrophic_methanogenesis-h2_dep", "methylotrophic_methanogenesis-h2_indep",
    }
    methanotrophy = "carbon_redox-dissimilatory_methanotrophy" in called
    # Parentheses ensure the methanotrophy exclusion applies to BOTH evidence routes.
    if (bool(called & methanogenesis) or ko_predict_methanogen(kos)) and not methanotrophy:
        primary = "methanogen"
    elif methanotrophy and has_methanotrophy_marker(kos):
        primary = "methanotroph"
    elif "Wood_Ljungdahl-acetogen" in called:
        primary = "homoacetogen"
    elif sum(p.startswith("CAZy") for p in called) >= 3 and sum("degradation" in p for p in called) >= 3:
        primary = "generalist"
    elif sum("fermentation" in p for p in called) >= 2:
        primary = "fermenter"
    elif any(p.startswith("CAZy") for p in called):
        primary = "macromolecule_degrader"
    elif any("degradation" in p for p in called):
        primary = "monomer_degrader"
    else:
        primary = ""
    return {
        "specialisation": primary,
        "nitrogen_specialisation": next((role for pathway, role in NITROGEN_ROLES if pathway in called), "None detected"),
        "sulfur_specialisation": next((role for pathway, role in SULFUR_ROLES if pathway in called), ""),
    }


def parse_kos(value: str) -> set[str]:
    """Split whole KO tokens, retaining legacy six-digit refined identifiers."""
    if value.strip() in {"", "-", "None", "nan"}:
        return set()
    tokens = set(re.split(r"[,;\s]+", value.strip()))
    if any(not re.fullmatch(r"K\d{5,6}", token) for token in tokens):
        raise ValueError(f"Invalid KO identifier(s): {value!r}")
    return tokens


def read_annotations(annotation_table, ko_column="accepted_ko", genome_id=None,
                     genome_column=None, gene_genome_map=None) -> dict:
    """Build genome -> gene -> feature sets, retaining unannotated genomes."""
    columns, rows = read_tsv(annotation_table)
    if ko_column not in columns:
        raise ValueError(f"Missing KO column {ko_column!r}")
    if sum(value is not None for value in (genome_id, genome_column, gene_genome_map)) > 1:
        raise ValueError("Choose one of genome_id, genome_column or gene_genome_map")
    mapping = {}
    genomes = {}
    if gene_genome_map:
        map_columns, map_rows = read_tsv(gene_genome_map)
        if not {"gene_id", "genome"} <= set(map_columns):
            raise ValueError("Gene/genome map requires gene_id and genome columns")
        for row in map_rows:
            if not row["gene_id"] or not row["genome"] or row["gene_id"] in mapping:
                raise ValueError("Empty or duplicate gene ID in gene/genome map")
            mapping[row["gene_id"]] = row["genome"]
            genomes.setdefault(row["genome"], {})
    elif genome_id is not None:
        if not genome_id.strip():
            raise ValueError("Genome ID cannot be empty")
        genomes[genome_id] = {}
    else:
        genome_column = genome_column or "genome"
        if genome_column not in columns:
            raise ValueError("Supply --genome-id, --genome-column or --gene-genome-map")
    if "gene_id" not in columns:
        raise ValueError("Annotation table requires gene_id (rename the DRAM index column when importing)")
    for row in rows:
        gene = row["gene_id"]
        if genome_id is not None:
            genome = genome_id
        elif gene_genome_map:
            genome = mapping.get(gene)
        else:
            genome = row[genome_column]
        if not gene or not genome:
            raise ValueError(f"Missing gene or genome mapping for {gene!r}")
        genes = genomes.setdefault(genome, {})
        if gene in genes:
            raise ValueError(f"Duplicate annotation for {genome}/{gene}")
        kos = parse_kos(row[ko_column])
        if ko_column == "accepted_ko" and len(kos) > 1:
            raise ValueError(f"accepted_ko must be a single assignment for {gene}; alternatives are not confirmed evidence")
        genes[gene] = kos
    if not genomes:
        raise ValueError("No genomes to assess")
    return genomes


def add_markers(genomes: dict, marker_table) -> None:
    """Import gene-level marker calls, replacing the corresponding ambiguous KO.

    TSV columns are gene_id, feature, and optionally genome. A missing genome is
    allowed only for globally unique gene IDs. feature accepts marker:pmoA, pmoA,
    or the legacy K109440. The source annotation table is never modified.
    """
    columns, rows = read_tsv(marker_table)
    if not {"gene_id", "feature"} <= set(columns):
        raise ValueError("Marker table requires gene_id and feature columns")
    locations = {}
    for genome, genes in genomes.items():
        for gene in genes:
            locations.setdefault(gene, []).append(genome)
    updates = {}
    for row in rows:
        gene = row["gene_id"]
        candidates = [row["genome"]] if "genome" in columns else locations.get(gene, [])
        if len(candidates) != 1 or gene not in genomes.get(candidates[0], {}):
            raise ValueError(f"Unknown or ambiguous marker gene: {gene!r}; supply its genome")
        name = row["feature"].removeprefix("marker:")
        feature = MARKER_IDS.get(name, name)
        if feature not in MARKER_IDS.values():
            raise ValueError(f"Unknown refined marker: {name!r}")
        key = (candidates[0], gene)
        previous = updates.setdefault(key, set())
        if any(len(old) == len(feature) == 7 and old[:6] == feature[:6] and old != feature for old in previous):
            raise ValueError(f"Conflicting refinements for {gene}: {previous} and {feature}")
        previous.add(feature)
    for (genome, gene), features in updates.items():
        # Mirrors AnnoGuild's replacement, rather than keeping both AMO and pMMO
        # evidence via the ambiguous parent KO after a specific refinement.
        parents = {feature[:6] for feature in features if len(feature) == 7}
        genomes[genome][gene] = (genomes[genome][gene] - parents) | features


def read_product(path) -> tuple[list[str], dict]:
    """Import an existing refined Boolean matrix, preserving its genome order."""
    columns, rows = read_tsv(path)
    genome_column = "genome" if "genome" in columns else "#genome"
    if genome_column not in columns:
        raise ValueError("Product requires genome or #genome column")
    names = [column for column in columns if column != genome_column]
    products = {}
    for row in rows:
        genome = row[genome_column]
        if not genome or genome in products:
            raise ValueError("Empty or duplicate genome in product")
        products[genome] = {name: parse_bool(row[name]) for name in names}
    return names, products


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def specialise(annotation_table, output_dir, *, genome_id=None, genome_column=None,
               gene_genome_map=None, ko_column="accepted_ko", markers=None,
               product_table=None, cazy_product=None, reaction_threshold=0.7,
               ko_threshold=0.6, reactions_file=None, pathways_file=None) -> dict[str, Path]:
    """Run the complete genome-level analysis; no databases or network required."""
    for threshold in (reaction_threshold, ko_threshold):
        if not math.isfinite(threshold) or not 0 <= threshold <= 1:
            raise ValueError("Coverage thresholds must be finite and within [0, 1]")

    # 1. Aggregate the selected accepted annotations and optional refinements.
    genomes = read_annotations(annotation_table, ko_column, genome_id, genome_column, gene_genome_map)
    if markers:
        add_markers(genomes, markers)
    features = {genome: set().union(*genes.values()) for genome, genes in genomes.items()}

    # 2. Evaluate bundled rules, or import an already evaluated legacy product.
    details = []
    if product_table:
        pathway_names, products = read_product(product_table)
        unknown = set(products) - set(genomes)
        if unknown:
            raise ValueError(f"Product genomes lack annotation inventory: {sorted(unknown)}")
    else:
        rules = load_rules(reactions_file, pathways_file)
        pathway_names = [rule.pathway_id for rule in rules]
        products = {}
        for genome in sorted(genomes):
            products[genome] = {}
            for rule in rules:
                result = evaluate_pathway(rule, features[genome], reaction_threshold, ko_threshold)
                # A successful rule is supported even if some optional refinements
                # are unavailable. An unsuccessful rule needing an unassessed
                # refinement is unresolved, not a claim that the pathway is absent.
                required = set(filter(None, result["required_refinements"].split(";")))
                assessed = markers is not None or required <= features[genome]
                result["assessment"] = "supported" if result["called"] else "not_detected" if assessed else "unresolved_refinement"
                products[genome][rule.pathway_id] = result["called"]
                details.append({"genome": genome, **result})

    # 3. CAZy calls must come from the legacy verified-signature matrix. They
    # cannot be inferred reliably from KO-only evidence or raw family presence.
    if cazy_product:
        names, cazy = read_product(cazy_product)
        names = [name for name in names if name.startswith("CAZy")]
        if not names or set(products) != set(cazy):
            raise ValueError("CAZy product must contain CAZy columns and exactly the assessed genomes")
        for name in names:
            if name not in pathway_names:
                pathway_names.append(name)
        for genome in products:
            products[genome].update({name: cazy[genome][name] for name in names})

    # 4. Apply the original priorities; preserve empty cells and 'None detected'.
    specialisations = []
    assessments = []
    cazy_assessed = any(name.startswith("CAZy") for name in pathway_names)
    for genome, calls in products.items():
        labels = classify_specialisation(calls, features[genome])
        specialisations.append({"genome": genome, **labels})
        assessments.append({
            "genome": genome,
            "cazy_assessment": "assessed" if cazy_assessed else "unassessed",
            "refinement_assessment": "imported_product" if product_table else "assessed" if markers is not None else "partial_or_unassessed",
            "supported_pathways": ";".join(name for name, call in calls.items() if call),
            "nitrogen_roles": ";".join(role for name, role in NITROGEN_ROLES if calls.get(name)),
            "sulfur_roles": ";".join(role for name, role in SULFUR_ROLES if calls.get(name)),
            "ko_methanogen": ko_predict_methanogen(features[genome]),
            "methanotrophy_marker": has_methanotrophy_marker(features[genome]),
        })

    # 5. Write legacy schemas separately from explanations and reproducibility data.
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {name: out / filename for name, filename in {
        "product": "kolach_product_refined.tsv",
        "specialisations": "kolach_specialisations.tsv",
        "pathways": "kolach_metabolic_pathways.tsv",
        "assessment": "kolach_specialisation_assessment.tsv",
        "metadata": "kolach_specialisation_metadata.json",
    }.items()}
    write_tsv(paths["product"], ["genome", *pathway_names],
              [{"genome": genome, **calls} for genome, calls in products.items()])
    write_tsv(paths["specialisations"], SPECIALISATION_COLUMNS, specialisations)
    write_tsv(paths["pathways"], DETAIL_COLUMNS, details)
    assessment_columns = list(assessments[0]) if assessments else ["genome"]
    write_tsv(paths["assessment"], assessment_columns, assessments)
    provenance = json.loads((DATA_DIR / "manifest.json").read_text(encoding="utf-8"))
    inputs = {"annotations": annotation_table, "genome_map": gene_genome_map, "markers": markers,
              "product": product_table, "cazy_product": cazy_product}
    rule_inputs = {} if product_table else {"reactions": reactions_file or DATA_DIR / "reactions.tsv",
                                           "pathways": pathways_file or DATA_DIR / "pathways.tsv"}
    provenance.update({
        "reaction_threshold": reaction_threshold, "ko_threshold": ko_threshold,
        "ko_column": ko_column, "genome_id": genome_id, "genome_column": genome_column,
        "evidence_policy": "selected KO column only; alternative_kos are never added",
        "pathway_evaluation": "imported_product" if product_table else "legacy_compatible",
        "inputs": {name: {"path": str(Path(path).resolve()), "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
                   for name, path in {**inputs, **rule_inputs}.items() if path is not None},
    })
    paths["metadata"].write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return paths
