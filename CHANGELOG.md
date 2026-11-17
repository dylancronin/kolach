# Changelog

## Unreleased — metabolic specialisation development branch

### Added

- Offline `kolach specialise` and `annotate --add-specialisations`.
- Explicit genome grouping by ID, annotation column, or gene/genome map.
- Pinned AnnoGuild rules: 225 reactions and 55 pathways.
- Commented expression/count evaluator replacing graphs and process pools.
  No new dependencies required.
- Legacy product and four-column guild outputs, plus separate evidence,
  assessment and reproducibility reports.
- Refined marker and verified CAZy product imports; existing-product mode.
- Example data, expected tables, walkthroughs, source attribution and parity tests.

### Corrected

- Missing comma concatenating K04480 and K14082.
- Methylotrophic classifier names not matching h2_dependent/h2_independent.
- Methanotrophy exclusion now applies to both methanogenesis evidence branches.
- Renamed misleading ko_fdh helper; retained its actual monooxygenase evidence.
- Signature graph nodes merging sibling OR groups under an AND: expressions now
  enforce independent strict signature groups.
- Explicit Boolean parsing replaces identity-based tests.
- Snakemake all is explicitly the default target to request final outputs.
- Removed tests/ from .gitignore so regression tests are tracked normally.
- Specialisation workflow parameters now trigger regeneration when genome
  grouping changes; supplemental files are read through their resolved paths.
- Numeric-looking genome IDs such as `001` survive Snakemake config parsing
  and child-job execution without losing leading zeros.

### Compatibility and limits

- Preserved 0.7/0.6 thresholds, partial reaction presence, conditional KO coverage,
  repeated KO occurrences, independent route maxima, original priorities,
  empty cells and None detected.
- All 1,485 captured default pathway comparisons agree. Example products and
  corrected-classifier tables agree byte-for-byte. Comammox signature diagnostics
  intentionally differ where the old graph merged nodes.
- Missing CAZy/refinement assessments are reported separately; KO-only outputs
  cannot reproduce calls depending on unavailable supplemental evidence.
- GraftM/dbCAN execution, raw taxonomy interpretation and raw CAZy distillation
  are outside this change. Supply interpreted marker/verified CAZy results.
- Initial fixtures executed unmodified legacy code with minimal adapters;
  independent verification with actual upstream libraries is now complete.

### Final validation — 2026-10-05

- Imported the verified bundle on `dev/metabolic-specialisation` and merged
  the existing post-base `main` commit's comments and regression tests into
  development. `main` remains at `db5a5ce`.
- Full unittest discovery in the normal editable Kolach environment: 89 tests
  passed (61 tracked tests plus 28 pre-existing local tests in the untracked
  `tests/test_fixes.py`). Includes 28 specialisation tests and a real optional
  Snakemake regression test for numeric IDs and genome-grouping reruns.
- Pinned AnnoGuild checkout `291b2329`: all 1,485 captured pathway calls and
  coverage/signature diagnostics independently verified with actual pandas,
  numpy, networkx, click, graphviz and pytest imports, without source adapters.
- Standalone full example legacy tables match expected outputs byte-for-byte;
  small-rule default and relaxed-threshold examples and existing-product mode
  also pass. Editable imports work without reinstalling.
- Wheel build and extracted-wheel execution pass; the installed package includes
  reaction/pathway definitions, manifest, GPL license and specialisation workflow.
- Snakemake 9.26.1 dry run and controlled execution pass using captured KOfam
  annotations, actual integration, marker/CAZy imports and all five reports.
  Unchanged inputs are up to date; changed grouping schedules regeneration.
  External annotation tools and live database downloads were not exercised.
- Python compilation and Git whitespace checks pass.
