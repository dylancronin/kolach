# Results — 2026-11-17

Annotation searches ran at Kolach `41affce` (`run_syncomm.sh`). The GraftM and
specialisation jobs were then re-run alone with the final code
(`rerun_specialisation.sh`, which regenerates its workflow config and records
the commit in `syncomm/kolach_commit.txt`); a Snakemake dry run listed only
`annotate_graftm` and `specialise_genomes`. The final re-run used `b1b3963` (CAZy rules extended to DRAM2's ecosystem
sheets); earlier re-runs at `f0ed054` and `5a6d08b` used the common sheet only.
Compute: 26.7 allocated CPU-hours over five jobs (`accounting.tsv`).

## 1. GraftM with HMMER 3.4 (`graftm_hmmer_check.json`)

GraftM 0.15.1 pip-installed beside HMMER 3.4 reproduced every placement that
was recorded with HMMER 3.2.1:

| Check | Placements | Identical |
|---|---:|---:|
| Temporal-paper test sets (`emerge_read_tax.tsv`, 8 packages) | 546 | 546 |
| SynComm, all 29,327 proteins, 9 packages (AnnoGuild's run) | 103 | 103 |

The same check against the shared `kolach` environment after GraftM was added
to it (2026-11-19) first gave 546/546 and 100/103 with pplacer `1.1.alpha23`,
which conda chose by default: three borderline narH/nxrB placements moved
between `Root` and `Root; NxrB` (`graftm_hmmer_mismatches_shared_env_alpha23.tsv`).
None were proteins Kolach gates. With pplacer pinned to `1.1.alpha19` (as in the
validation clone and the original GraftM/AnnoGuild environment) the shared
environment reproduced 546/546 and 103/103 (`graftm_hmmer_check_shared_env.json`).
`environment.yaml` now pins `pplacer=1.1.alpha19`.

### New-user installation (2026-11-19)

A fresh clone of the branch (`7a123c6`) and `conda env create -f environment.yaml`
into an empty prefix succeeded in about 4 minutes (packages cached) with
`pip check` clean. It provides kolach (editable), GraftM 0.15.1, HMMER 3.4,
pplacer/guppy 1.1.alpha19, mafft 7.526, FastTree 2.2.0, dbCAN 5.2.9,
Snakemake 9.27.2, eggNOG-mapper and DeepKOALA. Two environment fixes were
needed: `extern` is installed by pip (the conda-forge build is tagged cp38 and
fails `pip check`), and `snakemake-interface-common<1.23.1` (1.23.1 requires
`packaging>=26.1`, Snakemake 9 requires `<26`). In that environment,
`kolach download --databases graftm` installed the published bundle (SHA-256
matches), and GraftM on the SynComm proteins reproduced the validated
`kolach_graftm_markers.tsv` byte for byte. A clean Python 3.11 virtualenv with
`pip install -e ".[graftm]"` also passed `pip check` and ran `graftM --version`.

### Shared environment end to end (2026-11-19)

GraftM packages were installed into the shared database directory
(, via ; existing
databases untouched).  then ran
 in the shared 
environment on 29 Afacilis proteins (job 57663813, 4 CPUs, 00:09:45): real KOfam, dbCAN,
GraftM and specialisation completed without errors. GraftM refined the same three
proteins as the full SynComm run (narG, narH, qNor → norZ), and the genome was
named from the FASTA ().

## 2. Logic parity with AnnoGuild on real data

Kolach's evaluator and classifier on AnnoGuild's own SynComm inputs
(post-GraftM KO table; gene rows duplicated by AnnoGuild's DRAM/MicrobeAnnotator
merge were combined per gene):

| Check | Compared | Mismatches |
|---|---:|---:|
| L1: 55 pathways × 6 genomes vs `EMERGE_20260226_product_refined.tsv` | 330 | 0 |
| L2: guild labels vs `EMERGE_20260226_genome_specialisations.tsv` | 18 | 0 |

The 1,485 captured pathway snapshots were also re-verified against the real
AnnoGuild code at `291b232` (`scripts/verify_annoguild_parity.py`).

## 3. Kolach end to end versus AnnoGuild (SynComm)

Kolach used KOfam + eggNOG (fast) + DeepKOALA + dbCAN + GraftM on the same
29,327 proteins (the DRAM `genes.faa` files AnnoGuild used).

### GraftM (`syncomm_E1_graftm.tsv`)

68 proteins were refined by either pipeline: 3 identical, 4 different,
61 refined only by AnnoGuild. **None of the 61 carries the package KO in
Kolach**; AnnoGuild placed every protein with every package. Kolach annotates
them as, for example, formate dehydrogenase (K00123), periplasmic nitrate
reductase NapA (K02567), assimilatory sulfite reductases (K00392, K00381),
NADH nitrite reductase (K00362), DMSO/thiosulfate reductases and urease
(K01429). These are the off-target hits the gating removes.

The 4 differences: three proteins that AnnoGuild called both narG and narH
(Kolach: narG only; their KOs gate only the narG_nxrA package, so the narH
package never placed them), and one qNor protein (AnnoGuild: norB; Kolach:
norZ, by design).

### Pathways (`syncomm_E3_pathways.tsv`)

34 of 330 calls differ. Nitrite oxidation was called by AnnoGuild in 5 of 6
genomes, every time from the off-target nxrA/nxrB hits; Kolach calls none.
Gtundricola's nitrate reduction likewise rested on off-target narG/narH.
The remaining differences follow differences in KO calls (e.g. assimilatory
sulfate reduction, sugar degradation, fermentation routes), in both directions.

The qNor→norZ mapping changed no call here: Afacilis's qNor (norZ) is present
with nosZ but without nirK/nirS, so denitrification is not called either way.

### Labels (`syncomm_E3_labels.tsv`)

| Genome | AnnoGuild | Kolach | Kolach signature-only CAZy |
|---|---|---|---|
| Afacilis DSM649 | fermenter / nitrate_reducer / sulfur_oxidiser | fermenter / nitrate_reducer / sulfur_oxidiser | fermenter |
| Blathyri DSM26656 | generalist / denitrifier | fermenter / denitrifier | fermenter |
| Dlutea DSM19970 | generalist / None detected / sulfur_oxidiser | generalist / None detected / sulfur_oxidiser | generalist |
| Gtundricola DSM23138 | generalist / nitrate_reducer | generalist / None detected | generalist |
| Troseus DSM18391 | generalist / nitrite_oxidiser | generalist / None detected | generalist |
| run1 bin 001 | generalist / nitrate_reducer | generalist / nitrate_reducer | fermenter |

Adding DRAM2's ecosystem-sheet substrates (host glycans, β-glucans, marine
algal polysaccharides, galactomannan) changed one CAZy call and no label:
Gtundricola gained CAZy-Hemicellulose through one α-galactosidase pair
(`GH27;EC:3.2.1.22`, DRAM2's galactomannan rule). No SynComm genome matched any
host-glycan or marine-polysaccharide rule, so those groups are covered by unit
tests only.

Nitrogen labels changed only where AnnoGuild relied on off-target GraftM hits.
Blathyri drops from generalist because DRAM2's family/EC rules give it one
polymer substrate (starch/glycogen) where AnnoGuild's legacy CAZy product gave
four. CAZy categories differ between the two rule sets, so substrate columns
are compared only through the labels.

## Not covered

- No biological ground truth: these compare pipelines, not phenotypes.
- Six genomes only; no methanogens, methanotrophs, comammox or reverse-dsr
  organisms, so those GraftM refinements and labels were not exercised on
  real data beyond the temporal-paper placement checks.
- Much of the new GraftM logic was not exercised on real data: no package KO
  appeared among SynComm `alternative_kos`, and all 7 gated proteins were
  accepted-gated and refined. Alternative gating, `superseded`, `excluded`,
  `unresolved`/`not_placed` and `unresolved_refinement` are covered by unit
  and workflow tests only.
- GraftM check A grafted only sequences already placed under HMMER 3.2.1, so it
  cannot reveal new placements; check B (all SynComm proteins) can, and found none.
- The bundle was published on 2026-11-18 as GitHub release `graftm-packages-v1`
  (dylancronin/kolach). `download_graftm()` with the pinned URL and SHA-256
  installed it and validated all nine packages against the placement table.
