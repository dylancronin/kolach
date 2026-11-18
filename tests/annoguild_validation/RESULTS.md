# Results — 2026-11-17

Kolach commit for the final specialisation run: `f0ed054` (annotation searches
from the same run at `41affce`; only the GraftM and specialisation jobs were
re-run, as a dry run confirmed). Compute: 26.6 allocated CPU-hours over three jobs
(`accounting.tsv`).

## 1. GraftM with HMMER 3.4 (`graftm_hmmer_check.json`)

GraftM 0.15.1 pip-installed beside HMMER 3.4 reproduced every placement that
was recorded with HMMER 3.2.1:

| Check | Placements | Identical |
|---|---:|---:|
| Temporal-paper test sets (`emerge_read_tax.tsv`, 8 packages) | 546 | 546 |
| SynComm, all 29,327 proteins, 9 packages (AnnoGuild's run) | 105 | 105 |

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
(Kolach: narG only, its accepted package), and one qNor protein (AnnoGuild:
norB; Kolach: norZ, by design).

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
- `kolach download --databases graftm` from the GitHub release URL is untested
  (the release does not exist yet); the bundle was installed from the local
  file through the same code path.
