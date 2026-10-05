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
- Initial fixtures execute unmodified legacy code with minimal adapters; a
  verifier using actual upstream libraries is included. The build environment
  lacks pandas/numpy and Snakemake, so the full old suite and live workflow need
  validation in an installed Kolach environment.
