# Integration regression checks

Run from the repository root with Python 3.11+ and its declared dependencies:

    PYTHONPATH=src python -m unittest discover -s tests -v

On PowerShell, set `$env:PYTHONPATH = "src"` before running Python.
The exhaustive suite takes a few minutes on the development machine.

`reference/integrate_db5a5ce.py` is an unmodified snapshot of
`src/kolach/integrate.py` at upstream commit
`db5a5ceb84c380d15dad18f77af1ae73d0d88643`.
Its Git blob SHA is `e851aa8a9861e8bd44822f6f821bd90aadba3edc`.
It is test-only and is excluded from the installed package by the existing
`src` package configuration. Do not update it merely to make comparisons pass.

## Coverage

- 2,129,396 comparisons enumerate every per-method category assignment over three
  representative KOs. Categories are absent, passing, below threshold, unknown,
  and (for KOfam) heuristic rescue. Each tool can also be inactive. Active tools
  with no calls are distinct from inactive tools. Both conflict policies and
  both eggNOG filtering settings run for each combination.
- 123,904 comparisons enumerate zero, one or two eggNOG source rows, each with
  any nonempty subset of those KOs and any source status. All external confident
  and candidate support subsets and both policies are covered.
- Larger fixtures exercise six KOs and repeated groups.
- The seven nonempty combinations of example tools run under all four policy
  combinations. Comparisons include full DataFrames, dtypes, ordering, the audit
  table in attrs, and both TSV files byte for byte.
- Boundary cases cover duplicate ranking and ties, missing domain metrics,
  threshold equality, malformed and nonfinite values, zero E-values, gzip,
  metadata, empty inputs, duplicate FASTA headers, extra genes, definition
  precedence, errors, and the pipeline entry-point signature.
- KOfam score inference and duplicate ranking cover 1,152 rows, including missing
  and nonfinite metrics, full/domain profiles, explicit assignments, and large integers.

This is exhaustive coverage of the declared finite category model, not a claim
to enumerate every arbitrary-size annotation table or every numeric input.
The frozen reference preserves the prior policy. The comparison harness adds only
unselected heuristic-rescued KOfam KOs to its expected alternatives before
comparing summaries and evidence; accepted calls and consensus decisions still
match the original policy. Explicit characterization assertions
and boundary fixtures also provide independent expected behavior.
The reference side builds its own category indexes rather than borrowing
indexes computed by the new implementation.

## Gene processing

The production pipeline constructs one `GeneAnnotation` at a time and runs:

1. Categorize its borrowed records and index independent supporting methods.
2. Disambiguate eggNOG sources while retaining their provenance.
3. Resolve the consensus and alternatives.
4. Emit the summary and evidence rows.

The object owns these steps through `group_evidence`, `filter_eggnog_calls`,
`choose_annotation`, `summary_row`, and `evidence_rows`. Annotation selection
separates candidates, a single confident tool, agreement, and conflicts. Early
returns make the rule order explicit. The old gene-level module helpers have
been removed rather than retained as compatibility wrappers.

## Benchmark

    PYTHONPATH=src python tests/benchmark_integration.py --genes 20000 --repeats 3

The script first checks full output parity for the deterministic larger fixture,
then alternates old/new runs in fresh processes. It reports median integration
time and process peak resident memory (working set on Windows), including
interpreter/import memory. Measurements are informational, not fragile timing
assertions in the test suite.

To compare against the main version used for this refactor, save its source
outside the checkout and supply `--baseline`:

    git show 9d0573a:src/kolach/integrate.py > ../baseline_main.py
    PYTHONPATH=src python tests/benchmark_integration.py --baseline ../baseline_main.py --genes 20000 --repeats 3

The supplied module is used for both the exact output comparison and the timing.

Compared against main at `9d0573a` on Windows with Python 3.12.14, pandas 3.0.1
and NumPy 2.3.5, using three fresh processes per version and dataset:

| Dataset | Main time | New time | Main peak memory | New peak memory |
| --- | ---: | ---: | ---: | ---: |
| Examples (3,095 genes) | 1.042 s | 1.010 s | 84.57 MiB | 84.56 MiB |
| Synthetic (20,000 genes) | 9.202 s | 9.075 s | 221.04 MiB | 221.00 MiB |

Median runtimes were approximately 3% and 1% lower, with essentially unchanged
peak memory. These small differences are specific to these fixtures and machine,
not performance guarantees. KOfam score classification uses pandas masks, with
object comparisons to preserve Python's treatment of integer thresholds above 2**53.

