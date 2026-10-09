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
  precedence, errors, and the existing public helper signatures.

This is exhaustive coverage of the declared finite category model, not a claim
to enumerate every arbitrary-size annotation table or every numeric input.
The reference preserves the prior policy; explicit characterization assertions
and boundary fixtures also provide independent expected behavior.
The reference side builds its own category indexes rather than borrowing
indexes computed by the new implementation.

## Gene processing

The production pipeline constructs one `GeneAnnotation` at a time and runs:

1. Categorize its borrowed records and index independent supporting methods.
2. Disambiguate eggNOG sources while retaining their provenance.
3. Resolve the consensus and alternatives.
4. Emit the summary and evidence rows.

The object is then released. It stores inputs and decisions, while plain
functions implement the rules. Existing public helpers remain compatibility
entry points; the main pipeline uses the gene object directly.

## Benchmark

    PYTHONPATH=src python tests/benchmark_integration.py --genes 20000 --repeats 3

The script first checks full output parity for the deterministic larger fixture,
then alternates old/new runs in fresh processes. It reports median integration
time and process peak resident memory (working set on Windows), including
interpreter/import memory. Measurements are informational, not fragile timing
assertions in the test suite.

Measured on Windows with Python 3.12.14, pandas 3.0.1 and NumPy 2.3.5,
using three fresh processes per version and dataset:

| Dataset | Old time | New time | Old peak memory | New peak memory |
| --- | ---: | ---: | ---: | ---: |
| Examples (3,095 genes) | 0.729 s | 0.752 s | 107.60 MiB | 108.01 MiB |
| Synthetic (20,000 genes) | 6.558 s | 6.880 s | 241.04 MiB | 243.65 MiB |

These measurements show approximately 3% and 5% longer runtimes and approximately
0.4% and 1.1% higher peak memory. They are specific to these fixtures and machine,
not performance guarantees. The initial larger slowdown was traced to an extra
pandas Series scan in KOfam parsing and removed by iterating its columns directly.
