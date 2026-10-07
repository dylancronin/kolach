# DRAM2 metabolism rules

Run the current DRAM2 trait and product rules on Kolach annotations:

```sh
kolach specialise --ruleset dram2 --annotation-table kolach_annotations.tsv \
  --genome-id MAG001 --output-dir metabolism
```

Combined tables can use `--genome-column` or `--gene-genome-map`, as with the
existing specialisation command. `--dram2-system` selects `default`, `ag`, `bgc`,
`eng_sys`, `gut`, or `marine`. Custom tables use `--dram2-rules` and
`--dram2-common`. The AnnoGuild default and its output files remain available.

`kolach_dram2.tsv` contains one row per genome/rule/group. Boolean rules report
`present`. Modules and complexes report their units, units present, and coverage
as a fraction; these are **not** converted into Boolean calls with an invented
cutoff. Traits are independent labels, rather than AnnoGuild's priority-based
primary specialisation. The metadata JSON records source commits and input hashes.

Only `accepted_ko` contributes Kolach KO evidence. Supplemental evidence must be
explicit gene-level columns: `dbcan_id`, `dbcan_sub_id`, `dbcan_sub_ec`,
`dbcan_sub_substrate`, `fegenie_id`, `sulfur_id`, `camper_id`, `methyl_id`,
`merops_family`, `peptidase_family`, and DRAM-style Pfam/EC columns. A different
confirmed KO column can be selected with `--ko-column`. Conflicting or candidate
KO columns are never unioned into the evidence. This command evaluates supplied
annotations; it does not run these supplementary annotation tools.

DRAM2's taxonomy inclusion options are applied to the input annotation columns.
Without taxonomy, the combined methanogen/ANME trait is emitted instead of the
taxonomy-dependent methanogen subtypes. With taxonomy, its regex rules apply.
Missing annotations provide zero evidence, not proof of biological absence.
In particular, KO-only inputs cannot fully assess specialised iron/sulfur,
taxonomy, CAZy, or marker-dependent rules. Record which supplementary tools were
run alongside the results.

## Small implementation, unchanged data

The independent parser/evaluator is `src/kolach/metabolism/dram2.py` and uses only
the standard library. It parses square brackets, explicit AND/OR, comma steps,
aliases, thresholds, taxonomy, and row filters. It has no graph representation,
dynamic Python execution, or dependency on DRAM2 at runtime. The data under
`data/dram2/` is copied byte-for-byte from upstream Git objects with SHA-256 hashes
and source commits in `manifest.json`. Licenses are included with the assets.

The pinned public DRAM2 checkout is `4536cd93c97c07d783fab9e9ef99e56b271f3782`,
with Rule-Parser submodule `b87133283d4d07924007902ae49cc8632f7ad908`.
Product tables use dram-viz `18ae54a6eaf7f0887cfbc34432a44ab592f27e5f`.
The importer accepts new checkouts so future changes remain a reviewable data
diff and can be tested before inclusion:

```sh
git clone https://github.com/BortonWrightonLabs/DRAM.git upstream-dram
git -C upstream-dram checkout 4536cd93c97c07d783fab9e9ef99e56b271f3782
git -C upstream-dram submodule update --init bin/rule_parser
git clone https://github.com/BortonWrightonLabs/dram-viz.git upstream-viz
git -C upstream-viz checkout 18ae54a6eaf7f0887cfbc34432a44ab592f27e5f
python scripts/import_dram2_rules.py upstream-dram upstream-viz
pip install numpy polars lark networkx
PYTHONPATH=src python scripts/verify_dram2_parity.py upstream-dram --dram-viz upstream-viz
PYTHONPATH=src python -m unittest discover -s tests -p test_dram2.py
```

The independent verifier runs the actual upstream parser, without stubs. It
compares empty, complete, seeded random and every singleton feature profile
across all seven tables. Initial validation passed **530,929 comparisons** of
Boolean and numeric results. The checked-in smaller snapshots are generated
from upstream results with `--write-fixture` and checked offline by unit tests.
These are rule-engine checks, not a biological benchmark or a complete DRAM2
Nextflow annotation run. Input feature conversion is a separate boundary.

## Explicit compatibility choices

- Commas mean steps, not AnnoGuild OR alternatives. Mixed `&`/`|` requires brackets.
- `at_least(..., COUNTS, ...)` currently casts steps to Boolean upstream; Kolach
  preserves this behavior, including repeated genes not increasing the score.
- Shared aliases override local aliases. Duplicate names use the last expression,
  while each group row is retained, matching upstream table/dictionary behavior.
- Some upstream ecosystem sheets have ragged metadata. The loader ignores extra
  trailing cells and fills missing cells, matching DRAM2's CSV loader.
- Upstream accidentally exposes its last common alias as an unnamed output.
  Kolach omits that unnamed output; every named result is compared.
- The rule inclusion gate checks actual annotation columns. The upstream
  `evaluate_cycles` helper checks rule-table columns instead; this is a deliberate
  correction covered by input/output tests, separate from expression parity.

CAZy family-plus-EC rules depend on the upstream annotation encoding, so do not
infer an EC-specific substrate merely from a CAZy family. Preserve gene/domain
evidence and verify any external converter separately.
