#!/usr/bin/env python3
"""
Classify EMERGE genomes by metabolic specialisation.

The script merges *pathway-level evidence* (from DRAM's
`product_refined` table) with a *KO-pattern heuristic* that detects
methanogens while excluding methanotrophs.
"""
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate genome specialisations.")
    parser.add_argument("--product", required=True, help="Refined product TSV file with CAZy annotations")
    parser.add_argument("--annotations", required=True, help="DRAM annotations file (filled/updated)")
    parser.add_argument("--output", required=True, help="Output destination for specialisations TSV")
    return parser.parse_args()

# ------------------------------------------------------------------
# 1. ------------------------ KO SIGNATURES ------------------------
# ------------------------------------------------------------------
# Core machinery (present in all methanogens)
KO_MCR = {"K00399", "K00401", "K00402"}                    # mcrA/B/G
KO_HDR = {
    "K22480","K22481","K22482",
    "K03388","K03389","K03390",
    "K08264","K08265",
}

# CO2-reducing / hydrogenotrophic route
KO_FWD_FMD = {
    "K00200","K00201","K00202","K00203","K00204","K00205",
    "K11260","K11261",
}

# Methylotrophic route (methyl-amine / methanol)
KO_RAMX   = {"K16183",}#"K16184","K16185"}                  # activases
KO_MTXBCA = {                                             # methyltransferases
    "K14080","K14081","K04480"                            # mtaA-C - https://www.kegg.jp/module/M00356+K14080
    "K14082",                                             # mtbA
    "K16176","K16177","K16178","K16179","K14083","K14084", # mtt, mtb, mtm https://www.kegg.jp/module/M00563+K14082
    "K16954",                                             # mtsA
}

# Acetoclastic route
KO_ACETO = {"K01895","K00192","K00193","K00194","K00195","K00197"}
ACETO_MIN_HITS = 4            # >= 4 of the 6 above, no extra Mcr check

# Methanotrophy exclusion (pMMO / sMMO)
KO_METHANOTROPHY = {
    "K10944","K10945","K10946","K109440","K109450","K109460",   # pMMO / AMO
    "K16157","K16158","K16159","K16160","K16161","K16162" # sMMO
}

def ko_predict_methanogen(kos: set) -> bool:
    """Return True if KO profile matches *any* methanogenesis route
    and lacks methanotrophy mono-oxygenases."""
    if not kos:
        return False

    # Exclude clear methanotrophs
    if kos & KO_METHANOTROPHY:
        return False

    # Must have some core machinery (Hdr or Mcr)
    if not kos & (KO_HDR | KO_MCR):
        return False

    # Route-specific evidence
    has_h2_route   = len(kos & KO_FWD_FMD) >= 4
    has_methyl     = bool(kos & KO_RAMX) and bool(kos & KO_MTXBCA)
    has_acetate    = len(kos & KO_ACETO) >= ACETO_MIN_HITS

    return has_h2_route or has_methyl or has_acetate

def ko_predict_fdh(kos: set) -> bool:
    """Return TRUE if KO profile matches *any* FDH route."""
    if not kos:
        return False

    if (kos & KO_METHANOTROPHY):
        return True
    else:
        return False

def classify_specialisation(row: pd.Series) -> pd.Series:
    """Return primary + N/S specialisation calls for one genome."""

    called = {p for p, v in row.items() if v is True and p != "ko_methanogen"}

    # ----- Primary role (mutually exclusive) -----------------------
    if (
        any(p in called for p in [
            "hydrogenotrophic_methanogenesis-all",
            "acetoclastic_methanogenesis-all",
            "methylotrophic_methanogenesis-h2_dep",
            "methylotrophic_methanogenesis-h2_indep",
        ])
        or row["ko_methanogen"]
        and "carbon_redox-dissimilatory_methanotrophy" not in called
    ):
        primary = "methanogen"

    elif "carbon_redox-dissimilatory_methanotrophy" in called and row["ko_fdh"]:
        primary = "methanotroph"

    elif "Wood_Ljungdahl-acetogen" in called:
        primary = "homoacetogen"

    elif (
        sum(p.startswith("CAZy") for p in called) >= 3
        and sum("degradation" in p for p in called) >= 3
    ):
        primary = "generalist"

    elif sum("fermentation" in p for p in called) >= 2:
        primary = "fermenter"

    elif any(p.startswith("CAZy") for p in called):
        primary = "macromolecule_degrader"

    elif any("degradation" in p for p in called):
        primary = "monomer_degrader"

    else:
        primary = None

    # ----- Nitrogen / sulfur tags ---------------------------------
    n_tag = (
        "comammox"         if "nitrogen_redox-comammox"           in called else
        "ammonia_oxidiser" if "nitrogen_redox-ammonia_oxidation"  in called else
        "denitrifier"      if "nitrogen_redox-denitrification"    in called else
        "nitrate_reducer"  if "nitrogen_redox-nitrate_reduction"  in called else
        "nitrite_oxidiser" if "nitrogen_redox-nitrite_oxidation"  in called else
        "nitrogen_fixer"   if "nitrogen_redox-nitrogen_fixation"  in called else
        "None detected"
    )

    s_tag = (
        "sulfate_reducer"
        if "sulfur_redox-dissimilatory_sulfate_reduction" in called
        else "sulfur_oxidiser"
        if "sulfur_redox-dissimilatory_sulfur_oxidation" in called
        else None
    )

    return pd.Series({
        "specialisation": primary,
        "nitrogen_specialisation": n_tag,
        "sulfur_specialisation": s_tag,
    })

def main():
    args = parse_args()

    # 1. Product Matrix
    product = (
        pd.read_csv(args.product, sep="\t")
          .rename(columns=lambda c: c.strip())
          .rename(columns={"#genome": "genome"}, errors="ignore")
          .set_index("genome")
    )

    # 2. Annotations Matrix
    ann = pd.read_csv(args.annotations, sep="\t", dtype=str)
    ann["genome"] = ann["fasta"].str.replace(r"_annotate$", "", regex=True)

    ko_sets = ann.groupby("genome")["ko_id"].apply(lambda col: set(col.dropna()))

    # Apply predictors
    ko_meth_flag = ko_sets.apply(ko_predict_methanogen).rename("ko_methanogen")
    product = product.merge(ko_meth_flag.to_frame(),
                            how="left", left_index=True, right_index=True
                           ).infer_objects(copy=False).fillna({"ko_methanogen": False})

    ko_fdh_flag = ko_sets.apply(ko_predict_fdh).rename("ko_fdh")
    product = product.merge(ko_fdh_flag.to_frame(),
                            how="left", left_index=True, right_index=True
                           ).infer_objects(copy=False).fillna({"ko_fdh": False})

    # Classify
    specialisations = product.apply(classify_specialisation, axis=1)

    # Write output
    specialisations.reset_index().to_csv(args.output, sep="\t", index=False)
    print(f"Finished. Written specialisations to: {args.output}")

if __name__ == "__main__":
    main()
