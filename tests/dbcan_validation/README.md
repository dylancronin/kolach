# dbCAN / DRAM2 matched-database validation

The completed comparison is summarized in [RESULTS.md](RESULTS.md).

This checks annotation compatibility, not biological accuracy or a complete
DRAM2 annotation/distillation run. The reference is the unmodified DRAM2
`RUNDBCAN_CAZYMEANNOTATION` Nextflow process at commit
`4536cd93c97c07d783fab9e9ef99e56b271f3782`.

Both searches use dbCAN 5.2.9 and the same DIAMOND executable (2.2.6 here).
DRAM2's module environment specifies DIAMOND 2.2.5, but this comparison disables
environment/container creation so the same installed 2.2.6 binary is used by
both wrappers. This isolates workflow/conversion differences from tool versions.

The original 3,096 proteins are reused from `../anvio_validation/proteins.faa`.
Kolach prepares reversible search IDs to protect numeric IDs from upstream type
inference. DRAM2 receives that exact same sequence file, and the comparison maps
search IDs back to the original genes. No genes are predicted again.

## Reproduce

From the repository root, create an isolated Python 3.11 environment:

```sh
python -m venv tests/dbcan_validation/env
tests/dbcan_validation/env/bin/python -m pip install dbcan==5.2.9 polars scikit-bio
tests/dbcan_validation/env/bin/python tests/dbcan_validation/prepare_reference.py
curl -L --fail https://github.com/nextflow-io/nextflow/releases/download/v24.10.4/nextflow-24.10.4-dist \
  -o tests/dbcan_validation/nextflow
chmod +x tests/dbcan_validation/nextflow
```

Put this environment's executables first on PATH, followed by the existing
Kolach environment (providing Snakemake, HMMER and DIAMOND). Set `PYTHONPATH` to
this checkout's absolute `src` path. Download the fixed database bundle once:

```sh
kolach download --database-dir tests/dbcan_validation/db --databases dbcan
sbatch tests/dbcan_validation/run.sh
```

Submit from the repository root. The job requests four CPUs on PDS0325, with
memory following the CPU allocation, and runs the searches sequentially.
The environment location in `run.sh` is specific to this installation.
Native Kolach/Nextflow rerun handling reuses completed searches after failures.
Do not resubmit completed searches merely to inspect results; `compare.py` can
be rerun directly with the same isolated Python and `PYTHONPATH`.

## Comparisons

- Compare complete family HMM, subfamily HMM and DIAMOND rows as multisets.
- Require the entire upstream overview, including recommendations, to agree.
- Require every native domain/alignment to survive in `dbcan_evidence.tsv`.
- Invoke DRAM2's original `combine_annotations.py` and verify that every retained
  CAZy value exists in Kolach. List additional values preserved by Kolach in
  `additional_domains.tsv`; these reflect DRAM2's first-value aggregation.
- Reuse the existing matched KOfam/eggNOG/DeepKOALA tables, require all original
  gene columns and KO evidence to remain identical after adding dbCAN, and
  check retention of all 1,577 anvi'o reference pairs.

`validation.json` records comparison counts, versions, database hashes and source
provenance. `reference_sources.json` pins downloaded reference files. `trace.tsv`,
Slurm accounting and `/usr/bin/time -v` logs record resource usage. Large native
outputs, database files, environments and logs are ignored by Git.
