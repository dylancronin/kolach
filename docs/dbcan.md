# dbCAN annotation and CAZy metabolism

This branch adds released **dbCAN 5.2.9** (the current successor to dbCAN3) as
an optional annotation method. It runs upstream's DIAMOND, dbCAN HMM and
dbCAN-sub HMM methods, then imports the upstream domain recommendations.
Search thresholds and overlap resolution stay in dbCAN, rather than being
reimplemented in Kolach.

Install in a Linux annotation environment with DIAMOND available:

```sh
pip install '.[cazy]'
kolach download --database-dir kolach-db --databases dbcan
kolach annotate --protein-fasta MAG001.faa --database-dir kolach-db \
  --databases kofam dbcan --output-dir results --threads 8 \
  --add-specialisations --genome-id MAG001
```

The downloader retrieves the four CAZyme references from upstream's pinned S3
release `db_v5-2-9_5-5-2026`: `CAZy.dmnd`, `dbCAN.hmm`, `dbCAN-sub.hmm`, and
`fam-substrate-mapping.tsv`. Files are streamed to temporary paths, checked for
nonempty contents, then published. A manifest records URLs, sizes and SHA-256s.
Reusing a file requires its recorded release and hash to match. These hashes
record downloaded bytes; upstream does not supply trusted checksums at these URLs.
Externally managed reference directories can also be used; the runner records
their actual hashes without claiming they are the pinned downloaded release.

`--databases dbcan` can run alone. Every protein still appears in the integrated
table; KOs remain unannotated until a KO method supplies evidence. dbCAN families
never vote in the existing KO consensus. The three raw method outputs and
`overview.tsv` remain under `results/dbcan/`, alongside command/version/reference
hash metadata. Each execution uses a fresh temporary directory so failed or empty
reruns cannot consume stale search hits.

## Import existing annotations

For already computed dbCAN 5 overview output:

```sh
kolach cazy --annotation-table kolach_annotations.tsv --overview overview.tsv \
  --genome-id MAG001 --output-dir cazy-results
kolach specialise --annotation-table cazy-results/kolach_cazy_annotations.tsv \
  --genome-id MAG001 --cazy-product cazy-results/kolach_cazy_product.tsv \
  --output-dir metabolism
```

`kolach cazy --protein-fasta ... --database-dir kolach-db/dbcan` runs the searches
instead of importing an overview. `kolach integrate --dbcan-overview ...` attaches
CAZy columns during an existing table integration. `kolach specialise
--dbcan-overview ...` performs the import and CAZy calculation before the legacy
guild classification. Combined genome tables can use `--genome-column` or an
explicit `--gene-genome-map`. dbCAN joins require globally unique gene IDs;
rename duplicate IDs before annotating a mixed protein FASTA.

The importer reads columns by name and rejects unknown/duplicate genes, invalid
family names and inconsistent tool counts. At least two tools plus an upstream
recommendation are required for confirmed family evidence. Single-method hits
remain in candidate columns. Proteins missing from a successfully computed
overview receive zero calls while remaining in the annotation inventory.
An imported overview is an assertion that the upstream search completed; Kolach
cannot establish that from the table alone. Its import metadata records hashes
of both the original annotation table and overview.

Normalized evidence includes `dbcan_id` (base families), `dbcan_sub_id` (selected
subfamilies), EC/substrate fields, tool counts and raw candidate columns.
`dbcan_family_ec` retains explicit family/EC pairs from selected subfamily domains.
Pairs follow the matching ordered domain/EC lists in dbCAN's overview, rather
than a cross-product of families and ECs. For multi-name domains, compound pairs
are omitted when their association is ambiguous. Substrate summaries are omitted
if any subfamily domain was discarded by upstream recommendations. The raw
overview remains available for auditing that conservative choice.

After merging both metabolism PRs, the normalized annotation table can feed
`kolach specialise --ruleset dram2`, including its newer CAZy product rules.
The two branches are independently based on `dev/metabolic-specialisation`.

## Fill the original CAZy metabolism gap

`kolach_cazy_product.tsv` feeds the existing `--cazy-product` input directly.
The bundled definitions are copied unchanged from AnnoGuild commit
`291b232907e6568597178e21384ff0f964993a82`. A substrate call requires every
DRAM1 cleavage row plus an observed family unique to that substrate, excluding
AnnoGuild's low-specificity list. This preserves the original signature gate;
raw family presence alone does not imply a macromolecule-degrader guild.

`kolach_cazy_pathways.tsv` explains cleavage counts, signature evidence and
supporting features. `kolach_cazy_metadata.json` records the input and definition
hashes and grouping parameters. Existing specialisation files retain their
original schemas, with CAZy now marked assessed when this evidence is supplied.

## Reproduce validation

```sh
PYTHONPATH=src python -m unittest discover -s tests
git clone --branch v5.2.9 https://github.com/bcb-unl/run_dbcan.git upstream-dbcan
curl -L https://raw.githubusercontent.com/dylancronin/AnnoGuild/291b232907e6568597178e21384ff0f964993a82/scripts/append_cazy.py -o append_cazy.py
PYTHONPATH=src python scripts/verify_dbcan_parity.py upstream-dbcan append_cazy.py
PYTHONPATH=src python scripts/smoke_dbcan_search.py
```

The parity script checks exact upstream source commits/blobs and runs dbCAN's
real overview generator over nine synthetic search profiles, including all
method combinations and overlapping/multiple domains. It also passed **3,382**
comparisons against AnnoGuild's original CAZy gate. Checked-in snapshots support
offline regression tests. Unit tests cover failure-safe reruns/downloads,
gene inventory, weak hits, signatures, ambiguous KO exclusion and domain/EC
association.

The Linux CI smoke test builds tiny synthetic HMM and DIAMOND databases, runs
all three actual searches through Snakemake, and checks integrated annotations
and a CAZy metabolic call. Synthetic sequences are a plumbing test, not a
biological benchmark. Production reference searches need a suitable Linux
environment and the full downloaded references.
