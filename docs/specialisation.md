# Metabolic specialisation: walkthrough

`kolach specialise` aggregates accepted KOs by genome, evaluates EMERGE pathway
rules, and applies AnnoGuild's guild priorities. It runs offline without new
dependencies or rerunning annotation. `kolach pathway` remains a separate
gene-level KEGG BRITE lookup.

The implementation has two commented modules: `metabolism/evaluate.py` parses
and scores rules; `metabolism/specialise.py` handles inputs, classification, and
outputs. Bundled rules contain 225 reactions and 55 pathways, pinned to AnnoGuild
commit `291b232907e6568597178e21384ff0f964993a82`.

## 1. Feed in annotations and assign genomes

For one MAG:

```bash
kolach specialise \
  --annotation-table results/kolach_annotations.tsv \
  --genome-id MAG001 --output-dir results/specialisation
```

For a combined table, use its `genome` column automatically, or select another
column with `--genome-column bin_id`. Example input:

```tsv
genome	gene_id	accepted_ko	alternative_kos
MAG001	gene_1	K00001	-
MAG001	gene_2	K00003	-
MAG001	gene_3	-	K00002
MAG002	gene_4	-	-
```

This becomes:

```python
{"MAG001": {"K00001", "K00003"}, "MAG002": set()}
```

`alternative_kos` are not confirmed functions and are never added to evidence.
The default `accepted_ko` column must contain at most one KO per gene. For an
intentional legacy import, `--ko-column ko_id` supports comma/semicolon-separated
calls but does not independently resolve their ambiguity.

Alternatively, provide `--gene-genome-map genes_to_genomes.tsv`:

```tsv
gene_id	genome
gene_1	MAG001
gene_2	MAG001
gene_3	MAG002
```

Choose one grouping method. Names are never guessed by trimming gene identifiers.
All annotation genes must have a mapping. Duplicate map IDs are rejected;
genomes present only in the map remain in the output with no evidence. Gene IDs
may repeat across genomes when the annotation itself has a genome column.

## 2. Feed in reaction and pathway rules

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

## 3. Calculate coverage from example evidence

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
is **False**. Run the included small example:

```bash
kolach specialise \
  --annotation-table examples/specialisation/small_annotations.tsv \
  --reactions-file examples/specialisation/reactions.tsv \
  --pathways-file examples/specialisation/pathways.tsv \
  --output-dir small_output
```

Expected product:

```tsv
genome	example_degradation-all
SMALL	False
```

With `--reaction-threshold 0.66`, the call becomes True. This demonstrates a
threshold change; it is not a recommendation to alter the scientific method.

Legacy choices deliberately preserved:

- A multi-subunit reaction is partially present when any subunit is found.
- Entirely absent OR alternatives are excluded from KO-coverage candidates.
- Repeated KOs in different reactions count separately.
- Reaction and KO coverage use independent maxima over alternative routes. They
  may come from different routes; the report explicitly labels `ko_coverage_route`.
- A required signature must pass regardless of coverage.

The evaluator combines distinct `(observed, required)` counts across reactions,
avoiding enumeration of every full enzyme-alternative combination. The
`missing_features` report lists missing KOs in the selected KO-coverage route.
For an entirely absent reaction it lists all alternative candidates, without
adding them to the denominator. With no observed route it lists all pathway
candidates. Missing signature features likewise include missing alternatives;
an OR signature can pass while some alternatives are absent.

## 4. Feed in refined marker evidence

AnnoGuild used GraftM-derived custom IDs to distinguish functions sharing a KO.
These are preserved internally; six-digit custom IDs are not official KOs.

Before refinement:

```tsv
genome	gene_id	accepted_ko
MAG001	pmo_gene	K10944
MAG002	amo_gene	K10944
```

Supply `--markers markers.tsv`:

```tsv
genome	gene_id	feature
MAG001	pmo_gene	marker:pmoA
MAG002	amo_gene	marker:amoA
```

Internal features become `K109440` and `K109441`. The ambiguous parent `K10944`
is replaced for each gene, mirroring AnnoGuild; the original annotation file
stays unchanged. This prevents an AMO refinement retaining ambiguous methane
monooxygenase evidence via its parent KO.

| Names | Legacy IDs |
|---|---|
| pmoA / amoA | K109440 / K109441 |
| pmoB / amoB | K109450 / K109451 |
| pmoC / amoC | K109460 / K109461 |
| narG / nxrA | K003700 / K003701 |
| narH / nxrB | K003710 / K003711 |
| dsrA / rdsrA | K111800 / K111801 |
| dsrB / rdsrB | K111810 / K111811 |
| norB / norZ | K045610 / K045611 |
| mmoX | K16157 |

`feature` accepts `marker:pmoA`, `pmoA`, or `K109440`. The genome column is
optional only for globally unique gene IDs. Unknown genes/features and
contradictory refinements of one parent KO on one gene are rejected.

This imports already interpreted GraftM or equivalent calls. It does not run
GraftM or interpret raw taxonomy strings. Use the exact Kolach gene IDs and
provide complete results for the assessed genomes. A header-only marker table
explicitly means refinement was assessed with no positive results.

Without refinements, an unsupported pathway needing custom IDs is labelled
`unresolved_refinement` in the detailed report; its legacy product cell remains
False. A pathway with enough available alternative evidence remains supported.

## 5. Feed in verified CAZy calls

The original generalist/macromolecule-degrader decisions require CAZy evidence.
Supply `--cazy-product` with the verified Boolean CAZy columns from AnnoGuild's
`product_refined_with_CAZy.tsv`. The full old matrix is accepted: only CAZy
columns are imported. Example:

```tsv
genome	CAZy-Starch	CAZy-Chitin	CAZy-Pectin
MAG001	True	True	True
MAG002	False	False	False
```

The table must contain exactly the assessed genomes, preventing missing CAZy
assessments being silently interpreted as negatives. Raw dbCAN family presence
is insufficient: AnnoGuild also uses its verified signature-family checks and
DRAM product cutoff. This feature imports the completed calls; running dbCAN
and reconstructing its distillation layer are outside this implementation.

Without CAZy input, the assessment report says `unassessed`. Other primary labels
may still be returned, but CAZy-dependent guild priorities remain incompletely
assessed. Do not interpret the legacy False cells as demonstrated absence.

## 6. Convert supported pathways to guild labels

Primary priority remains:

1. Methanogen: a supported methanogenesis route or the original KO heuristic,
   excluding a supported dissimilatory methanotrophy pathway.
2. Methanotroph: dissimilatory methanotrophy plus the original monooxygenase check.
3. Homoacetogen: Wood–Ljungdahl acetogen pathway.
4. Generalist: three CAZy calls and three degradation calls.
5. Fermenter: two fermentation pathway calls.
6. Macromolecule degrader: any CAZy call.
7. Monomer degrader: any degradation call.
8. Otherwise an empty primary cell.

The methanogen heuristic retains some Mcr/Hdr evidence plus four Fwd/Fmd KOs,
or RamA with a methyltransferase, or four acetoclastic-route KOs. Its original
monooxygenase exclusion is retained. These are historical metabolic-potential
heuristics, not a newly validated biological classifier.

Example function input and output:

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

## 7. Run the full included example and inspect outputs

The synthetic example profiles illustrate the rules; they are not reference
genomes or biological benchmarks. The CAZy calls are illustrative supplied calls.

```bash
kolach specialise \
  --annotation-table examples/specialisation/annotations.tsv \
  --markers examples/specialisation/markers.tsv \
  --cazy-product examples/specialisation/cazy_product.tsv \
  --output-dir example_output
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

| File | Contents |
|---|---|
| `kolach_product_refined.tsv` | Legacy Boolean matrix: genome, 55 original pathways, optional CAZy columns |
| `kolach_specialisations.tsv` | Exactly genome, specialisation, nitrogen_specialisation, sulfur_specialisation |
| `kolach_metabolic_pathways.tsv` | Coverage, signatures, route, supporting/missing evidence and assessment |
| `kolach_specialisation_assessment.tsv` | Supplemental assessment status, supported pathways, all N/S roles, heuristic flags |
| `kolach_specialisation_metadata.json` | Source/version, thresholds, policy and SHA-256 hashes of inputs and rules |

Booleans remain `True`/`False`. Missing primary and sulfur labels are empty;
missing nitrogen retains `None detected`. File names have a Kolach prefix;
legacy schemas are preserved and additional reports are separate.
Expected legacy tables under `examples/specialisation/expected/` are checked
byte-for-byte against the new outputs.

## 8. Import an old product or run after annotation

To compare classifiers independently, prepare legacy annotations with `gene_id`,
`genome`, and `ko_id`, retaining custom six-digit IDs. Rename the DRAM index
column to `gene_id`. If old genome IDs end in `_annotate`, remove that suffix
consistently to match the product; Kolach uses explicit IDs without auto-trimming.

```bash
kolach specialise \
  --annotation-table legacy_annotations_normalized.tsv --ko-column ko_id \
  --product EMERGE_product_refined_with_CAZy.tsv --output-dir comparison
```

`--product` preserves product order, accepts `#genome`, and skips pathway
recalculation. Product genomes must exist in the annotation inventory;
annotation-only genomes are excluded, matching the old left join. The detailed
pathway file has only a header because no coverage was calculated. Metadata
records the import mode.

To run after protein annotation:

```bash
kolach annotate \
  --protein-fasta MAG001.faa --database-dir kolach-db --output-dir results \
  --add-specialisations --genome-id MAG001
```

For a mixed FASTA replace `--genome-id` with `--gene-genome-map`. Optional markers
and CAZy inputs work here too. Snakemake tracks those files, annotations, and
rules, along with genome grouping parameters. Changing the genome ID triggers
regeneration; IDs such as `001` retain their leading zeros. Its `all` rule is
explicitly the default target. Supplemental paths are resolved before execution.

## Corrections and validation

[CHANGELOG.md](../CHANGELOG.md) records the comma, methylotrophic naming,
methanotrophy exclusion, misleading `ko_fdh` name, signature node collisions,
and explicit Boolean parsing corrections. No broader scientific thresholds,
guild priorities, denominators, or default route policies were changed.

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_metabolism.py -v
```

The standard-library tests cover 27 profiles across 55 original pathways (1,485
Boolean comparisons), coverage, the corrected original classifier, errors,
markers, ambiguity, CLI configuration and example outputs. Comammox signature
diagnostics fix the old node-collision bug; all captured default Boolean calls
agree. Changing thresholds can expose a corrected signature call difference.

Snapshots were generated by executing the original RuleParser source unchanged
with minimal adapters for graph, scalar-array, and lookup operations, because
the build environment could not install pandas/numpy/networkx. Fixture provenance
records that limitation. `scripts/verify_annoguild_parity.py` independently verifies
the snapshots with the real upstream libraries. The full existing suite requires
the usual installed Kolach environment: `python -m unittest discover -s tests -v`.

Final validation on 2026-10-05 verified all 1,485 snapshots against the pinned
AnnoGuild checkout with actual dependencies. The full local suite passed 89 tests
(61 tracked and 28 pre-existing local tests); wheel assets and extracted-wheel
example execution also passed. No editable reinstall was needed.

To reproduce the controlled offline workflow test:

```bash
python -m unittest tests.test_metabolism.TestSpecialisationWorkflow -v
```

It supplies one protein and a captured KOfam call (`g1`, `K10944`), imports the
`pmoA` marker, and runs the real integration and specialisation jobs. It checks
that genome `001` appears exactly in the product, an unchanged dry run has no
work, and changing the genome to `002` reruns the reports and changes the output
ID. Annotation tools and databases are unnecessary for this fixture test. The
test skips when Snakemake is unavailable; the core feature tests remain
standard-library-only.
