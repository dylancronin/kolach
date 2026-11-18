# Metabolic specialisation

`kolach annotate --add-specialisations` predicts genome-level metabolism and
AnnoGuild's guild labels from Kolach's own outputs. It needs no manual inputs:

1. KO annotation and consensus (KOfam/eggNOG/DeepKOALA, as selected).
2. dbCAN, added automatically, for CAZy substrate evidence.
3. GraftM, which refines KOs shared by several functions (pmoA/amoA,
   dsrA/rdsrA, narG/nxrA, ...).
4. EMERGE pathway rules and AnnoGuild's guild priorities, per genome.

`kolach pathway` remains a separate gene-level KEGG BRITE lookup.

Code: `methods/graftm.py` (marker refinement), `metabolism/cazy.py` (CAZy
substrates), `metabolism/evaluate.py` (rule parsing and scoring) and
`metabolism/specialise.py` (inputs, classification, outputs). Bundled pathway
rules contain 225 reactions and 55 pathways, pinned to AnnoGuild commit
`291b232907e6568597178e21384ff0f964993a82`.

## 1. Run it

```bash
kolach download --database-dir kolach-db --databases kofam dbcan graftm
kolach annotate --protein-fasta MAG001.faa --database-dir kolach-db \
  --output-dir results --add-specialisations --threads 8
```

The genome is named after the FASTA file (`MAG001` here; `.gz` and one FASTA
suffix are removed). Use `--genome-id` to choose another name. For a FASTA
holding several genomes, give `--gene-genome-map` (columns `gene_id`, `genome`);
every gene must be mapped. Names are never guessed from gene IDs. Running a
directory of genome FASTAs in one call is planned.

`--add-specialisations` adds `dbcan` to `--databases` if it is missing, and a
gene/genome map is checked against the FASTA before any search starts.
`environment.yaml` installs the GraftM tools (`graftM`, `pplacer`, `guppy`,
`hmmsearch`, `hmmalign`); update an existing environment with
`conda env update -f environment.yaml`. GraftM is pip-installed so HMMER stays
at 3.4; its bioconda recipe pins 3.2.1. HMMER 3.4 reproduced all recorded
HMMER 3.2.1 placements (see `tests/annoguild_validation/`).

To re-evaluate a finished run, e.g. with other thresholds, without repeating
any search:

```bash
kolach specialise --annotation-dir results --output-dir results_ko05 --ko-threshold 0.5
```

It reads `kolach_annotations.tsv` (with its dbCAN columns) and
`kolach_graftm_markers.tsv`, reuses the genome grouping recorded by the annotate
run (or `--genome-id` / `--gene-genome-map`), and writes to a separate
`--output-dir` so the workflow's own reports are never overwritten.

## 2. Which KOs count

Each gene contributes its `accepted_ko`. `alternative_kos` are not confirmed and
are not used directly; an alternative KO counts only after GraftM places the
protein in a functional clade (section 5).

```tsv
gene_id	accepted_ko	alternative_kos
gene_1	K00001	-
gene_2	K00003	-
gene_3	-	K00002
```

gives genome evidence `{K00001, K00003}`.

## 3. Reaction and pathway rules

The small included reaction example (`examples/specialisation/reactions.tsv`):

```tsv
pathway	reaction_id	definition
example_degradation	1	K00001+K00002
example_degradation	2	K00003
example_degradation	3	K00004
```

The pathway example (`examples/specialisation/pathways.tsv`):

```tsv
pathway	reaction	subpathway	signature_definition
example_degradation	1+2+3	all	K00001
```

The grammar supports `+`, comma-separated alternatives, and parentheses:

```python
parse_expression("(K00001+K00002),K00003")
# OR(AND(K00001, K00002), K00003)
```

One representation serves three different scoring policies:

| Context | Tokens | Interpretation |
|---|---|---|
| Reaction | KOs/refined IDs | Any evidence makes a reaction partially present; missing AND subunits affect KO coverage |
| Pathway | Reaction numbers | Each alternative route receives reaction coverage |
| Signature | KOs/refined IDs | Strict AND/OR requirements |

Rules are parsed once per run. Invalid syntax, duplicate IDs, and missing
reaction references are rejected. The implementation uses expressions and count
combinations, with no graph library, graph rendering, process pool, or `eval`.

## 4. Coverage, using example evidence

Feed `{K00001, K00003}` into the three example reactions:

| Reaction | Required | Observed | Partial presence |
|---|---|---|---|
| 1 | K00001 + K00002 | K00001 | True |
| 2 | K00003 | K00003 | True |
| 3 | K00004 | None | False |

Reaction coverage is `2/3 = 0.6667`. KO coverage is also `2/3`: it counts two
observed KOs out of the three required KOs in reactions having some evidence.
The entirely absent reaction 3 is excluded from the KO denominator.

Default requirements are 0.7 reaction coverage and 0.6 KO coverage, so the call
is **False**. Run the included small example with the Python API:

```python
from kolach.metabolism.specialise import specialise
specialise("examples/specialisation/small_annotations.tsv", "small_output",
           reactions_file="examples/specialisation/reactions.tsv",
           pathways_file="examples/specialisation/pathways.tsv")
```

Expected product:

```tsv
genome	example_degradation-all
SMALL	False
```

With `reaction_threshold=0.66` the call becomes True. This demonstrates a
threshold change; it is not a recommendation to alter the scientific method.

Legacy choices deliberately preserved:

- A multi-subunit reaction is partially present when any subunit is found.
- Entirely absent OR alternatives are excluded from KO-coverage candidates.
- Repeated KOs in different reactions count separately.
- Reaction and KO coverage use independent maxima over alternative routes. They
  may come from different routes; the report explicitly labels `ko_coverage_route`.
- A required signature must pass regardless of coverage.

The `missing_features` report lists missing KOs in the selected KO-coverage
route. For an entirely absent reaction it lists all alternative candidates,
without adding them to the denominator. With no observed route it lists all
pathway candidates. A failed pathway that needs a refined feature is labelled
`unresolved_refinement` when a protein carrying that feature's gate KO could not
be refined (GraftM unresolved or no hit), or when no marker table exists;
otherwise it is `not_detected`.

## 5. Marker refinement with GraftM

Some KOs cover several functions. EMERGE GraftM packages place a protein in a
curated tree, and the clade names its function. Refined features use AnnoGuild's
encoding (shared KO + 0/1); they are not official KOs.

| Package | Gate KO | Functions (feature IDs) |
|---|---|---|
| amoA_pmoA | K10944 | pmoA K109440 / amoA K109441 |
| amoB_pmoB | K10945 | pmoB K109450 / amoB K109451 |
| amoC_pmoC | K10946 | pmoC K109460 / amoC K109461 |
| narG_nxrA | K00370 | narG K003700 / nxrA K003701 |
| narH_nxrB | K00371 | narH K003710 / nxrB K003711 |
| dsrA | K11180 | dsrA K111800 / rdsrA K111801 |
| dsrB | K11181 | dsrB K111810 / rdsrB K111811 |
| norZ_norB | K04561 | cNor: norB K045610 / qNor: norZ K045611 |
| mmoX_sMMO | K16157 | mmoX K16157 |

The method follows the temporal-paper workflow
(`GraftM-analysis/update_annotations.R`), not AnnoGuild's `process_graftm.py`:

- **Gating.** A protein is placed with a package only if the gate KO is its
  accepted KO or one of its alternative KOs. AnnoGuild placed every protein with
  every package, letting a package's HMM classify unrelated enzymes.
- **Exact lookup.** Every clade of every package is listed in
  `metabolism/data/graftm_placements.tsv` with its function, or a reason it has
  none. Placements are matched exactly by (package, placement). Installing the
  packages checks that the packages and the table list the same clades.
- **One function per protein.** Paralogous packages (dsrA/dsrB, narG/narH,
  amoA/B/C) can all place one protein when their KOs are among its
  alternatives. The package gated by the accepted KO wins and the others are
  `superseded`; if no accepted-gated package refines it, competing refinements
  are all `unresolved`.
- **Effect on evidence.** A refined accepted KO is replaced by its feature
  (K10944 becomes K109440). A protein confirmed from an alternative KO gains the
  feature beside its accepted KO. A protein placed among non-target homologues
  (`excluded`: the mmoX outgroup `Homologues_mmoX`, i.e. sHMO homologues, and
  `Nod`, nitric oxide dismutase) loses its accepted gate KO, so it no longer
  supports e.g. methanotrophy. Unresolved and unplaced proteins change nothing:
  the shared KO stays, and rules needing a refined feature ignore it.

Unresolved placements: the root and `unknown*`/`Unknown` clades.

Corrections relative to earlier code:

| Source | Problem | Now |
|---|---|---|
| AnnoGuild substring matching | `"dsrA" in "Root; rdsrA; …"`, so every rdsrA/rdsrB became dsrA/dsrB (249 of the 1,414 dsrA/dsrB placements in the temporal paper's `emerge_all_read_tax.tsv`) | rdsrA/rdsrB |
| AnnoGuild substring matching | `Root; Homologues_mmoX` (outgroup) became mmoX | excluded |
| AnnoGuild (all proteins placed) | One protein could be refined by several paralogous packages | one function per protein |
| AnnoGuild substring matching | Clades absent from the EMERGE set (e.g. `NarG_Archaea`) worked only if a substring happened to match | every clade listed |
| Temporal-paper lookup | Joined on placement alone, so a bare `Root` in any package became norB | unresolved |
| Temporal-paper lookup | qNor clades mapped to norB (K045610). The norB reaction requires norC (K02305), which single-subunit qNor lacks; the qNor reaction is K045611 alone | qNor → norZ (K045611) |

The qNor mapping is a scientific decision recorded here for review.

`kolach_graftm_markers.tsv` has one row per gated protein and package:
`gene_id, package, gate_ko, gate_source (accepted|alternative), placement,
marker, feature, status (refined|excluded|unresolved|superseded|not_placed),
unresolved_reason`. GraftM's logs are kept in `kolach_graftm_markers.tsv.log`.
The assessment report counts each status per genome (rows, i.e. protein and
package pairs).

Packages: eight EMERGE packages from the temporal paper and the public GraftM
`7.60.mmoX_sMMO` package, distributed as one reproducible bundle
(`scripts/build_graftm_bundle.py`; SHA-256 pinned in `methods/graftm.py`).

## 6. CAZy substrates from dbCAN

Substrate definitions are DRAM2's current CAZy rules (`rules_common.tsv` from
[dram-rules](https://github.com/BortonWrightonLabs/dram-rules) at
`3257b8604ab6e40c00e811d6ba0211b44dcf5648`, MIT licence), not AnnoGuild's
legacy `DRAM_CAZy_definitions.tsv`. Each rule is an OR list of
`FAMILY;EC:x.x.x.x` pairs. A genome has a substrate when any protein carries
one of those pairs in its `dbcan_family_ec` column.

- Pairs come from one dbCAN-sub domain. DRAM2 instead pairs every family on a
  protein with every EC on it.
- Subfamilies are reduced to families (GH5_2 → GH5), as in Kolach's dbCAN parser.
- Otherwise matching is exact and case-insensitive, as in DRAM2. `3.2.1.-`
  matches only `3.2.1.-`; provisional ECs such as `2.4.1.n2` are accepted.

Only polymer substrates count toward the generalist (three or more) and
macromolecule-degrader labels. Subcategories of one polymer are merged so they
cannot count twice. **This grouping is provisional.**

| Product column | DRAM2 rules |
|---|---|
| CAZy-Cellulose | Cellulose |
| CAZy-Hemicellulose | Hemicellulose-Xylan, Hemicellulose-Other |
| CAZy-Pectin | Pectin |
| CAZy-Arabinogalactan | Arabinogalactan |
| CAZy-Cutin | Cutin |
| CAZy-Fructan | Fructan |
| CAZy-Starch_glycogen | StarchSpecific, StarchGlycogen |
| CAZy-Chitin | Chitin |

Not counted: `Microbially-derived` (peptidoglycan; it includes GH23/GH102–104
lytic transglycosylases that nearly every bacterium uses on its own cell wall,
so it would make almost every genome a macromolecule degrader),
`PlantPolymerMixed` (not substrate-specific), oligosaccharides,
disaccharides, β-glucosides and oligosaccharide oxidation. They still appear in
`kolach_cazy_substrates.tsv` as `CAZy_other-*` rows. DRAM2 rules that mix KOs
with family/EC pairs (cellobiose, sucrose, trehalose, lactose, other
disaccharides) are not evaluated. AnnoGuild's generalist threshold was three of
18 legacy categories; three of eight is stricter.

`kolach_cazy_substrates.tsv` lists every substrate per genome with the matched
pairs and genes, plus `called_signature_only`: the call using only pairs unique
to one substrate (over all DRAM2 family/EC rules, as AnnoGuild did over its
categories) and outside AnnoGuild's low-matching families (GH63, GH49, GH127,
CE2, GH23, PL30, PL8, PL33, GH88, GH67, GH98, CE4). The assessment report's
`specialisation_cazy_signature_only` gives the primary label under those
signature-only calls. Both are for comparison and never affect the guild labels.

## 7. Convert supported pathways to guild labels

Primary priority remains:

1. Methanogen: a supported methanogenesis route or the original KO heuristic,
   excluding a supported dissimilatory methanotrophy pathway.
2. Methanotroph: dissimilatory methanotrophy plus the original monooxygenase check.
3. Homoacetogen: Wood–Ljungdahl acetogen pathway.
4. Generalist: three CAZy substrate calls and three degradation calls.
5. Fermenter: two fermentation pathway calls.
6. Macromolecule degrader: any CAZy substrate call.
7. Monomer degrader: any degradation call.
8. Otherwise an empty primary cell.

The methanogen heuristic retains some Mcr/Hdr evidence plus four Fwd/Fmd KOs,
or RamA with a methyltransferase, or four acetoclastic-route KOs. Its original
monooxygenase exclusion is retained. These are historical metabolic-potential
heuristics, not a newly validated biological classifier.

```python
classify_specialisation(
    {"ethanol_fermentation-all": True, "lactate_fermentation-all": True},
    {"K00656", "K00132", "K00001", "K00016"},
)
# {'specialisation': 'fermenter',
#  'nitrogen_specialisation': 'None detected',
#  'sulfur_specialisation': ''}
```

Nitrogen priority remains comammox, ammonia oxidiser, denitrifier, nitrate
reducer, nitrite oxidiser, then nitrogen fixer. Sulfate reduction retains
priority over sulfur oxidation. All supported N/S roles appear separately
in the assessment report.

## 8. Outputs

| File | Contents |
|---|---|
| `kolach_product_refined.tsv` | Legacy Boolean matrix: genome, 55 pathways, eight CAZy polymer columns |
| `kolach_specialisations.tsv` | Exactly genome, specialisation, nitrogen_specialisation, sulfur_specialisation |
| `kolach_metabolic_pathways.tsv` | Coverage, signatures, route, supporting/missing evidence and assessment |
| `kolach_specialisation_assessment.tsv` | Evidence sources, GraftM status counts, signature-only CAZy label, supported pathways, all N/S roles, heuristic flags |
| `kolach_cazy_substrates.tsv` | Every DRAM2 substrate per genome, with matched pairs and genes |
| `kolach_graftm_markers.tsv` (+ `.log`) | GraftM outcome per gated protein and package; GraftM logs |
| `kolach_specialisation_metadata.json` | Sources (AnnoGuild and dram-rules commits, GraftM bundle manifest), thresholds, policies and SHA-256 hashes of inputs and rules |

Booleans remain `True`/`False`. Missing primary and sulfur labels are empty;
missing nitrogen retains `None detected`. Legacy schemas are preserved; the
additional reports are separate files.

## 9. Bundled examples and legacy imports (Python API)

The synthetic profiles illustrate the rules; they are not reference genomes.
Their marker and CAZy inputs are supplied by hand. The Python API still allows
that for comparison with AnnoGuild; the CLI always computes both.

```python
from kolach.metabolism.specialise import specialise
specialise("examples/specialisation/annotations.tsv", "example_output",
           markers="examples/specialisation/markers.tsv",
           cazy_product="examples/specialisation/cazy_product.tsv")
```

| Genome | Primary | Nitrogen | Sulfur |
|---|---|---|---|
| EMPTY | Empty | None detected | Empty |
| FERMENTER | fermenter | None detected | Empty |
| GENERALIST | generalist | None detected | Empty |
| HOMOACETOGEN | homoacetogen | None detected | Empty |
| MACROMOLECULE | macromolecule_degrader | None detected | Empty |
| METHANOGEN | methanogen | None detected | Empty |
| METHANOTROPH | methanotroph | None detected | Empty |
| MONOMER | monomer_degrader | None detected | Empty |
| NITROGEN | Empty | comammox | Empty |
| SULFATE | Empty | None detected | sulfate_reducer |

The expected tables under `examples/specialisation/expected/` are checked
byte-for-byte. `product_table=` classifies an existing AnnoGuild product matrix
and `cazy_product=` imports AnnoGuild's verified CAZy columns; both exist for
parity checks. A `markers=` table needs `gene_id` and `feature` (`pmoA`,
`marker:pmoA` or `K109440`); with a `status` column only `refined` rows count.
Other options: `genome_column=` for a table with a genome column,
`ko_column="ko_id"` for legacy DRAM tables (comma/semicolon-separated KOs).

## Corrections and validation

[CHANGELOG.md](../CHANGELOG.md) records the classifier corrections (comma,
methylotrophic naming, methanotrophy exclusion, `ko_fdh` name, signature node
collisions, explicit Boolean parsing) and the GraftM corrections above. No
pathway thresholds, guild priorities, denominators or route policies changed.

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The pathway tests cover 27 profiles across 55 original pathways (1,485 Boolean
comparisons) captured from AnnoGuild; `scripts/verify_annoguild_parity.py`
re-checks them against a pinned AnnoGuild checkout with its real libraries.

The offline workflow test runs the real Snakemake DAG with stand-in `run_dbcan`
and `graftM` executables:

```bash
python -m unittest tests.test_metabolism.TestSpecialisationWorkflow -v
```

A captured KOfam call (K10944) on protein `001` in `001.faa` is placed as pmoA,
dbCAN evidence gives a starch substrate, and the genome is named `001` from the
FASTA. An unchanged rerun does nothing; `--genome-id 002` regenerates the reports.
GraftM gating, lookup, bundle reproducibility and archive safety are covered by
`tests/test_graftm.py`; CAZy matching by `tests/test_metabolism.py`.

Real-data checks are under `tests/annoguild_validation/` (see its README).
