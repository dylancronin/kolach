# dbCAN integration validation — 2026-11-09

Final verification passed **59 tests**, superseding the original implementation's
57-test count below; the rebuilt wheel passed 15 tests. The real-data comparison
was reproduced unchanged after the output-safety and recovery fixes.

The annotation foundation on `dev/dbcan-integration` reproduces DRAM2's dbCAN
search calls and recommendations while retaining additional domain annotations
in the integrated gene table. KO consensus and evidence are unchanged.

## Matched search comparison

Input: the same 3,096 anvi'o-predicted proteins used for the earlier KO check.
Protein FASTA SHA-256:
`fdd7e441e6f0d870f8da7df00c1c26a03a5480f25a042d3ad8ca5d65ca47529d`.
The reference uses DRAM2 commit `4536cd93c97c07d783fab9e9ef99e56b271f3782`,
its unmodified dbCAN Nextflow process, and its unmodified annotation combiner.
Source-file hashes are in `reference_sources.json`.

Both workflows used dbCAN 5.2.9, DIAMOND 2.2.6, four threads, and the same four
database files from S3 release `db_v5-2-9_5-5-2026`. DRAM2's environment file
names DIAMOND 2.2.5; using the same host 2.2.6 executable in both runs deliberately
controls for that version difference. Both receive identical amino-acid sequences
and reversible search IDs, which Kolach maps back to the original numeric IDs.

| Method | Kolach records | DRAM2 records | Missing | Additional |
|---|---:|---:|---:|---:|
| Family HMM | 106 | 106 | 0 | 0 |
| Subfamily HMM | 108 | 108 | 0 | 0 |
| CAZy DIAMOND | 146 | 146 | 0 | 0 |

All native rows agree as multisets, and the entire upstream overview agrees.
All 360 records survive in `dbcan_evidence.tsv`. All 3,096 original genes retain
their order in the normalized table; 169 have at least one hit and 96 have an
upstream recommended call. These are annotation results, not validated biological
phenotypes.

## Conversion and KO preservation

DRAM2's combiner retains only the first non-null value in each annotation field.
Kolach retains **18 additional values across three proteins**: three family IDs,
five subfamily IDs, and ten composition values. Every value retained by DRAM2
also appears in Kolach. See `additional_domains.tsv` for the exact differences.
No additional native search hits were introduced.

The pre-existing matched KOfam, eggNOG and DeepKOALA search tables were reused.
After adding dbCAN, every original KO summary column and every KO evidence row
is identical to integration without dbCAN. All **1,577** anvi'o reference gene/KO
pairs remain present across accepted and alternative calls. No KO search or
gene-prediction job was repeated.

## Tests, packaging and resources

- Full suite: **57 tests passed** in 334 seconds, including the exhaustive
  integration parity checks. The frozen reference was not edited; the public
  signature check permits only the new optional trailing `dbcan_tsv=None`.
- Latest dbCAN-specific suite: **14 tests passed**, including a real Snakemake
  DAG with executable stand-ins, mixed KOfam/dbCAN integration, pathway enrichment,
  unchanged-input reuse, database-change reruns, and missing-asset detection.
- The built wheel contains the wrapper and both dbCAN rules; **13 tests passed**
  against the extracted wheel. It excludes metabolic-specialisation code.
- The real download workflow recognizes the validated bundle as up to date.
- Slurm job **57193888** completed successfully on PDS0325: four CPUs, 18,224 MiB
  allocated, 15m58s elapsed, peak batch RSS approximately 3.0 GiB. Kolach took
  8m58s including database validation and provenance hashing; the reference
  workflow took 6m52s. `accounting.tsv` records scheduler measurements.

`validation.json` contains the machine-readable comparison and database hashes.
The scripts and reproduction instructions are in `README.md`. Large environments,
databases, native search outputs and logs remain available locally but are ignored
by Git. The specialisation engines and a biological accuracy benchmark are outside
this annotation-foundation change.
