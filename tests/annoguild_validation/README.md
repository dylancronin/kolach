# AnnoGuild specialisation: real-data validation

Checks for `kolach annotate --add-specialisations` against AnnoGuild and the
temporal-paper GraftM workflow. Results are in [RESULTS.md](RESULTS.md).

| File | Purpose |
|---|---|
| `check_graftm_hmmer.py`, `.sh` | GraftM 0.15.1 with HMMER 3.4 versus recorded HMMER 3.2.1 placements (temporal-paper test sets; SynComm, all proteins) |
| `compare_syncomm.py` | Logic parity on AnnoGuild's own SynComm features (L1, L2) and Kolach end-to-end comparison (E1–E3) |
| `run_syncomm.sh` | Full Kolach run on the six SynComm genomes, then `compare_syncomm.py` |
| `graftm_hmmer_check.json`, `syncomm_comparison.json`, `syncomm_*.tsv` | Saved results |

Local, untracked folders (see `.gitignore`):

- `env/`: clone of the shared `kolach` environment plus the GraftM tools from
  `environment.yaml` and `pip install graftm==0.15.1` (HMMER stays 3.4).
- `code/`: worktree pinned to the commit used for the SynComm run.
- `db/`: symlinks to the shared `/fs/scratch/PDS0325/Cronin/kolach-db`
  (kofam, eggnog, deepkoala, dbcan, ko00001.json; not modified) plus `graftm/`,
  installed from the locally built bundle through `download_graftm()`.
- `syncomm/`: `proteins.faa` (the six DRAM `genes.faa` files AnnoGuild used;
  29,327 proteins, same gene set as its `all_genes.faa`), `gene_genome_map.tsv`,
  and `kolach_output/`.

Reference inputs (read only):

- AnnoGuild SynComm run: `/fs/scratch/PDS0325/Cronin/syncomm/annoguild`
  (2026-02-26, before AnnoGuild `291b232`).
- Temporal-paper GraftM analysis:
  `/fs/scratch/PDS0325/Cronin/temporal-paper/temporal_git/temporal_paper_202505/GraftM-analysis`.
