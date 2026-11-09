# dbCAN / CAZy annotation

## Design and reference implementation

Kolach uses the three searches in dbCAN 5.2.9: family HMMs, subfamily HMMs,
and DIAMOND against CAZy. This is annotation of supplied proteins; no genes,
clusters, secretion signals, or metabolic guilds are predicted here.

The reference is DRAM2 commit `4536cd93c97c07d783fab9e9ef99e56b271f3782`:

- [dbCAN module](https://github.com/BortonWrightonLabs/DRAM/blob/4536cd93c97c07d783fab9e9ef99e56b271f3782/modules/nf-core/rundbcan/cazymeannotation/main.nf)
  runs `run_dbcan CAZyme_annotation --mode protein`, with dbCAN 5.2.9 and
  DIAMOND 2.2.5 in its environment specification.
- [Table conversion](https://github.com/BortonWrightonLabs/DRAM/blob/4536cd93c97c07d783fab9e9ef99e56b271f3782/bin/combine_annotations.py)
  imports the family and subfamily HMM tables, then selects the first non-null
  value per gene. Multiple domains can consequently disappear from its summary.
  Kolach preserves every reported domain and summarizes distinct values.
- The old `assets/internal/prepare_databases.py` still names V12 resources;
  these are not the complete modern run_dbcan bundle and are not used here.

The modern [database definitions](https://github.com/bcb-unl/run_dbcan/blob/v5.2.9/dbcan/constants/databases_constants.py)
contain four CAZyme assets: `CAZy.dmnd`, `dbCAN.hmm`, `dbCAN-sub.hmm`, and
`fam-substrate-mapping.tsv`. The remote subfamily filename is `dbCAN_sub.hmm`;
the downloader saves it as `dbCAN-sub.hmm`. The fixed S3 release is
`db_v5-2-9_5-5-2026`; the HTTP `db_current` source is a moving snapshot.

We reuse the upstream searches, overlap filtering, and recommendations instead
of reimplementing them. All passing hits contribute annotations, independently
of the separate multi-method recommendation. Family/EC pairs are extracted only
within a subfamily domain; unrelated domains are never cross-joined. Evidence
is distinct from KO consensus and does not add a fourth KO voter.

These outputs are intended for future DRAM2 specialisation rules. CAZy families
alone do not establish substrate usage, and this feature does not reconstruct
AnnoGuild's verified CAZy product or import its specialisation classifier.

## Usage

Create the environment from the updated `environment.yaml` (dbCAN 5.2.9 is pinned).
Existing installations without dbCAN continue to support all KO-only commands.
For an existing installation, install dbCAN into a separate environment and put
its `run_dbcan` executable on PATH alongside Kolach's tools.

```sh
kolach download --databases dbcan --database-dir kolach-db
kolach annotate --databases kofam dbcan --protein-fasta proteins.faa \
  --database-dir kolach-db --output-dir results --threads 4
```

`--databases dbcan` alone also works: all input genes remain in the integrated
table, with empty KO calls. Existing defaults remain KOfam-only. Searches use
the upstream default E-values (DIAMOND `1e-102`, HMMs `1e-15`) and HMM coverage
(`0.35`). Only positive thread allocations are accepted. The downloader writes
its success marker after validating the HMM, DIAMOND, and mapping files. The
annotation DAG requires all four assets; it does not silently download them.

Download metadata (`dbcan/manifest.json`) records source URLs, the fixed release,
checksums and tool version. Run metadata records input/database hashes, versions,
command, thresholds and the evidence policy. No existing shared database is
modified during validation.

## Outputs and standalone integration

| Output | Contents |
|---|---|
| `dbcan_annotations.tsv` | One row per input protein, in FASTA order |
| `dbcan_evidence.tsv` | Every reported domain/alignment, including its full native row as JSON |
| `dbcan_metadata.json` | Versions, hashes, parameters and evidence policy |
| `dbcan/` | Native result files, search input, ID map and log |
| `kolach_annotations.tsv` | Existing KO fields followed by optional dbCAN summary columns |

The summary columns are `dbcan_id`, `dbcan_sub_id`, `dbcan_sub_composition`,
`dbcan_sub_ec`, `dbcan_sub_substrate`, `dbcan_family_ec`, `dbcan_cazy_id`,
`dbcan_methods`, `dbcan_tools_count`, and `dbcan_recommended`.
Distinct values are sorted and semicolon-separated, except correlated family/EC
pairs (pipe-separated, e.g. `GH13;EC:3.2.1.1|GH5;EC:3.2.1.4`). Recommendations
are copied from the upstream overview without reinterpreting its policy.
Missing values are `-`; method counts are zero for genes without hits.
Method names and counts describe the protein as a whole; they do not assert
that every reported family is supported by every method.

Subfamily EC frequency suffixes (`:count`) are preserved in `dbcan_sub_ec`.
They are removed only from explicit family/EC features. A domain with multiple
possible families does not generate correlated pairs. `dbcan_cazy_id` retains
DIAMOND family IDs separately from the family HMM evidence in `dbcan_id`.

dbCAN's pandas readers can coerce numeric-looking identifiers. Kolach uses
reversible search IDs, preserving every amino-acid sequence, and restores the
original IDs in both normalized tables. Native files retain the search IDs;
`dbcan/gene_ids.tsv` maps them to the original IDs, including leading zeros.
Failed searches retain their temporary `dbcan-run-*` directory and log.

Keep the database directory separate from the annotation output directory, and
keep input proteins outside the native `dbcan/` output folder. Kolach rejects
overlapping database/output paths and file aliases that could overwrite inputs.
Native output folders carry a Kolach ownership marker; an unrelated nonempty
folder is rejected. Reruns replace generated files individually and preserve
other files added to a previously owned folder. Missing evidence, metadata and
declared native outputs are regenerated by the annotation workflow.

Import a normalized summary alongside existing KO search tables:

```sh
kolach integrate --protein-fasta proteins.faa \
  --kofam-table results/kofam_annotations.tsv \
  --dbcan-table results/dbcan_annotations.tsv --output-file combined.tsv
```

Python callers can pass the trailing optional `dbcan_tsv` argument to
`integrate_annotations()`. dbCAN never changes KO votes, accepted calls,
alternatives or the KO evidence table. Omitting dbCAN preserves the old schema.
With FASTA supplied, unknown dbCAN genes are rejected; otherwise the supplied
tables' gene inventories are combined. Duplicate summary gene IDs are rejected.
KEGG pathway enrichment preserves the supplementary columns.

The specialisation engines are not included on this branch. The DRAM2-compatible
columns and correlated features are the interface for a later integration;
raw CAZy family presence must not be treated as a verified AnnoGuild CAZy product.

## Validation

Run `PYTHONPATH=src python -m unittest discover -s tests -v`.
The dbCAN tests need no database downloads. The optional Snakemake test runs
the real DAG with deterministic executable stand-ins and checks reruns and
pathway enrichment. Live reference validation is documented under
`tests/dbcan_validation/` and invokes DRAM2's unmodified dbCAN Nextflow process.
The [completed comparison](../tests/dbcan_validation/RESULTS.md) matched all
360 native records and upstream recommendations on the 3,096-protein test input,
while preserving 18 summary values dropped by DRAM2's combiner.
