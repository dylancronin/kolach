# Slurm executor validation (`--executor slurm`)

Run on OSC Pitzer on 2026-11-19, with the shared `kolach` env and `snakemake-executor-plugin-slurm` 2.6.1 (pip). `PYTHONPATH` pointed to this checkout's `src`. Per-job accounting is in `accounting.tsv`. The annotation outputs are not committed (see `.gitignore`); their SHA-256 prefixes are listed below.

Input: the anvi'o validation proteins (`../anvio_validation/proteins.faa`, 3,096 genes) with the shared database dir `/fs/scratch/PDS0325/Cronin/kolach-db`. The tool settings match `../anvio_validation/run.sh`: `--eggnog-sensmode fast --eggnog-dmnd-block-size 4 --eggnog-dmnd-index-chunks 2 --deepkoala-release 202608 --device cpu`.

## 1. Parallel annotation (`run.sh`, controller 57663838)

The controller had 2 CPUs and ran `kolach annotate --databases kofam eggnog deepkoala --executor slurm --slurm-account PDS0325` with the bundled profile. Each tool ran as its own Slurm job, and all three started within 5 s of each other:

| Rule | Job | CPUs | Mem | Limit | Elapsed |
|---|---|---|---|---|---|
| annotate_kofam | 57663839 | 8 | 32000M | 6 h | 13:23 |
| annotate_deepkoala | 57663840 | 8 | 32000M | 6 h | 0:33 |
| annotate_eggnog | 57663841 | 16 | 62.5G | 12 h | 14:17 |

- `integrate_annotations` and `all` ran as localrules inside the controller. The whole run took 16:31.
- Every tool table is byte-identical to the original sequential run in `../anvio_validation/kolach_output/`:
  - kofam: `3f2766599313a64a`
  - eggnog: `f652bb2a0651fff5`
  - deepkoala: `679fe2eb51c477a4`
- `kolach_annotations.tsv` (`0552dc478e76c021`) and `kolach_evidence.tsv` (`1c2f1130849737cd`) are byte-identical to a local `kolach integrate --database-dir ...` of the original tables (`local_integrate/`).

This run used the code from before the review fixes. Runs 2–4 below use the final code.

## 2. Custom profile with `default-resources` (`run_final.sh`)

- **57663844 (cancelled after 14:39):** `final_profile/config.yaml` then set `slurm_extra: '--comment=...'`.
  - The plugin forbids executor-managed options in `slurm_extra`, such as `--comment`, `--qos` and `--mem`.
  - It raises only at submission, and the controller then hung at `Execute 1 jobs...` with no `sbatch` call. This was reproduced with a stub `sbatch`.
  - kolach now rejects these options before running Snakemake (`check_slurm_extra`).
- **57663847 (completed, 1:26):** the same run with `slurm_extra: '--mail-type=NONE'`.
  - DeepKOALA ran as job 57663850. Its table matches the original byte for byte (`679fe2eb51c477a4`).
  - Integration ran locally.
  - The per-run config file was removed afterwards, and no lock was left.

## 3. Cancelling the controller (`cancel_check/run_cancel.sh`)

In each case the controller was `scancel`led 30 s after its kofam job started running.

- **57663848, batch script runs `kolach` without `exec`:** Slurm killed kolach together with the batch shell. The kofam job 57663849 **kept running** until it was cancelled by hand.
- **57663851, batch script runs `exec kolach`:** kolach forwarded SIGTERM to Snakemake as SIGINT, and Snakemake cancelled kofam job 57663852 within the same second.
  - Snakemake's shutdown takes about 30 s, which loses to KillWait (30 s). The `.snakemake/locks` and the per-run config file therefore remained.
  - The README documents `exec` and the stale-lock recovery (`rm -r .snakemake/locks`).

`run_cancel.sh` is the `exec` version that was used for 57663851.

## 4. Stub-`sbatch` checks (no allocation)

A fake `sbatch` and `scancel` on `PATH` captured the submissions. The checks confirmed:
- CPUs, `--mem`, `-t`, `-A`, `-p` and `--gpus=1` for DeepKOALA with `--device cuda`.
- The `--configfiles` path that submitted jobs read.
- That light rules are never submitted.
- That `scancel` fires about 10 ms after SIGTERM reaches kolach.

An independent review found four medium issues: a lost profile `default-resources`, resource values evaluated as Python, a shared config file overwritten by a second run, and missing cancellation cleanup. All four are fixed and covered by `tests/test_cli.py` and `tests/test_workflows.py` (174 tests pass).
