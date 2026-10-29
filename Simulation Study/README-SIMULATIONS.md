# Simulation study (paper Section 3, Supplement Section S4)

Replicated study of the estimation error of the GP method and its competitors, for one covariate (`1d_skn`) and two covariates (`2d_hole`), n = 50, 250, 500, 1000 replicates, 50 datasets per scenario and n (400 datasets). The methods themselves are in `../R/` (see `../R/README-R.md`); this folder contains the scripts that run the study, its saved datasets, its written specification and its results.

Run every script from this folder (`Data Analysis/Simulation Study`), either with `Rscript` in the Terminal (quote the path because of the spaces) or from the R console (section "How to run from the R console"). Every runner can be interrupted: rerunning the same command resumes from the checkpoints, and completed fits are skipped.

## Scripts, in the order they are used

| Step | Script | What it does | Writes to `Results/` |
|---|---|---|---|
| 1 | `run_replicate_study.R` | The main study: for each dataset, the competitors (kernel with oracle and CV bandwidth, log-polynomials of degree 1-4) and one GP chain (60000 sweeps, 20000 burn-in, thinning 10); errors on a common reporting grid. Run in two stages: datasets 1-5, then 6-50. | `Replicate_Study/` |
| 2 | `run_long_fits.R` | Long two-chain fits of datasets 1-5 per scenario and n, extended until R-hat < 1.01 and bulk/tail ESS >= 200: reference fits and hyper-parameter summaries. | `Long_Fits/` |
| 3 | `validate_replicate_study_cap.R` | Paired check of the one-chain fits against the long fits (posterior means within 1%, error change < 0.002), to confirm that 60000 sweeps suffice. Requires all 40 comparisons (datasets 1-5 of every scenario and n) and stops with an error otherwise. No sampling. | `Replicate_Study/Validation/` |
| 4 | `run_sensitivity.R` | Sensitivity analysis (Supplement S4.5): one factor at a time around the replicate-study baseline, datasets 1-5, n = 50 and 1000 (300 fits, including the 20 baseline fits of the quadrature block S3). | `Sensitivity/` |
| 5 | `report_hyperparameters.R` | Hyper-parameter table of the paper from the long fits. No sampling. | `Paper_Tables/` |
| 6 | `plot_paper_figures.R` | Simulation figures of Section 3.3 from the replicate study and the long fits. No sampling. | `Paper_Figures/` (previews) |

`plot_paper_figures.R` also writes `sim1d_fits.pdf` and `sim2d_maps.pdf` directly into the manuscript folder (`Submitted Manuscripts/Round 2/BA 2026_10_xx/Figures_Main/`), replacing the figures currently in the paper.

Typical sequence (the commands are also in `SIMULATION_PROTOCOL.md`, Section 7):

```sh
Rscript run_replicate_study.R --smoke                                   # workflow check
Rscript run_replicate_study.R --cells=all --reps=1:5 --workers=4        # stage 1
Rscript run_long_fits.R --scenario=1d_skn --n=50 --reps=1:5 --workers=4 --compact   # every scenario and n
Rscript validate_replicate_study_cap.R
Rscript run_replicate_study.R --cells=all --reps=6:50 --workers=4       # stage 2
Rscript run_sensitivity.R --blocks=all --workers=4
Rscript report_hyperparameters.R
Rscript plot_paper_figures.R
```

## Protocols

The written specification of the study, in two files next to this README:

- `SIMULATION_PROTOCOL.md` (version 2.2): the specification of the study. What changed with respect to the first-round manuscript (new prior, extended kernel bandwidth grid, log-polynomial competitors, 20x20 GP grid in 2D), the design (scenarios, n, datasets, methods), the data generation and seeds, the prior and GP representation, the sampling settings of the three stages (replicate study, long fits, validation and its tolerances), the competitors and error metric, and the commands.
- `SENSITIVITY_SPEC.md` (version 3): the sensitivity analysis. Which referee points it answers, the one-factor-at-a-time design, the factors with their levels (S1 GP discretisation, S2 likelihood binning, S3 spatial quadrature, S4a truncation of rho*, S4b rate of the rho* prior, S5 tilt lambda, S6 rate b_ell of the length-scale prior), what is reported, and the commands.

## How to run from the R console

In a fresh R or RStudio session, set the working directory and define a helper `rs()` that runs a script with `Rscript` (`caffeinate -i` keeps the Mac awake while a run is going, as long as the lid is open; elsewhere use `rs <- function(...) system2(file.path(R.home("bin"), "Rscript"), c(...))`). Each `rs()` call waits until its script has finished, so several lines can be pasted at once and run in sequence. `--workers` is the number of fits (or chains) run in parallel, one core each.

```r
setwd("/Users/matteogiordano/Desktop/Research/Preprints/(2026, arXiv) Dolmeta, Giordano/Data Analysis/Simulation Study")
rs <- function(...) system2("caffeinate", c("-i", file.path(R.home("bin"), "Rscript"), c(...)))
```

**1. Replicate study** (one chain per dataset; stage 1 = datasets 1-5, stage 2 = datasets 6-50, into the same folder):

```r
rs("run_replicate_study.R", "--smoke")                                     # workflow check, a few minutes
rs("run_replicate_study.R", "--cells=all", "--reps=1:5", "--workers=3")    # stage 1, about 1.5 hours
rs("run_replicate_study.R", "--cells=all", "--reps=6:50", "--workers=3")   # stage 2
rs("run_replicate_study.R", "--report")                                    # rebuild the tables, no sampling
```

**2. Long fits** (two chains per dataset, extended until convergence; one call per scenario and n):

```r
rs("run_long_fits.R", "--smoke")
for (sc in c("1d_skn", "2d_hole")) for (n in c(50, 250, 500, 1000))
  rs("run_long_fits.R", paste0("--scenario=", sc), paste0("--n=", n), "--reps=1:5", "--workers=2", "--compact")
rs("run_long_fits.R", "--report")                                          # summary table over all finished fits
```

Each call fits the datasets of one cell with two chains each, so `--workers=2` keeps both chains of a dataset busy (more workers help only with several datasets per call). The 2D fits are the slow ones (200,000-300,000 sweeps per chain, 1-2 hours each with 2 workers); to use one more core, run the 2D n = 1000 call in a second session or in RStudio's Terminal tab (`caffeinate -i Rscript run_long_fits.R --scenario=2d_hole --n=1000 --reps=1:5 --workers=1 --compact`) and skip it in the loop.

**3. Validation and hyper-parameter table** (no sampling, a few seconds):

```r
rs("validate_replicate_study_cap.R")    # stops with an error unless all 40 comparisons are present and pass
rs("report_hyperparameters.R")
```

**4. Sensitivity analysis** (one chain per fit; datasets 1-5 at n = 50 and 1000; needs the replicate-study fits of those datasets as baselines):

```r
rs("run_sensitivity.R", "--smoke")
rs("run_sensitivity.R", "--blocks=all", "--reps=1:5", "--workers=3")      # all blocks, 300 fits (incl. the S3 baseline)
rs("run_sensitivity.R", "--blocks=S1,S2", "--reps=1:5", "--workers=3")    # or selected blocks
rs("run_sensitivity.R", "--report")                                        # rebuild the tables, no sampling
```

**5. Figures** (no sampling; overwrites the simulation figures in the manuscript folder):

```r
rs("plot_paper_figures.R")
```

A subset of datasets can be run first and the rest added later with the same command and a different `--reps`, into the same output folder: completed fits are skipped, and the report and validation are rebuilt over everything present. Changing any of the scripts or the files in `../R/` that a runner records in its `config.rds` requires a new `--out` folder, so keep the code unchanged until a study is complete.

## `Data/`

The saved simulation datasets, the input of every script. One file `<scenario>_n<n>_r<r>.rds` per dataset (400 for the study), plus `sens_<scenario>_m200_n<n>_r<r>.rds` (data simulated at 200x200 pixels for the quadrature block S3 of the sensitivity analysis) and a few tiny `n6` files from the smoke tests. Each dataset is simulated once from its seed, saved with a checksum (`data_id`) and reloaded by every method and every script, so all methods analyse exactly the same observations. The datasets do not depend on the prior. If deleted, they are regenerated identically.

## `Results/`

Empty for now: the folders below are created by the scripts. The results with the old prior (first-round manuscript) are in `../archive.zip` (`archive/simulation_results_old_prior/`).

```
Results/
├── Replicate_Study/        run_replicate_study.R (+ Validation/ from validate_replicate_study_cap.R)
├── Long_Fits/              run_long_fits.R
├── Sensitivity/            run_sensitivity.R
├── Paper_Tables/           report_hyperparameters.R
├── Paper_Figures/          plot_paper_figures.R
└── *_Smoke/                smoke tests of the three runners (can be deleted once they pass)
```

### `Replicate_Study/`

```
Replicate_Study/
├── config.rds, job_log.csv, STATUS.md, sessionInfo_*.txt
├── Source_Snapshot/        copy of the code used
├── Tables/
├── Validation/             cap_validation.csv, SUMMARY.md (validate_replicate_study_cap.R)
└── 1d_skn/ 2d_hole/ → n0050/ n0250/ n0500/ n1000/ → r01/ ... r50/
```

- `config.rds`: settings, prior and checksums of the code (a rerun after a code change stops and asks for a new `--out`); `job_log.csv`: one row per dataset job; `STATUS.md`: progress and error table at a glance.
- `Tables/`: `error_table.csv` (mean and SD of the relative L2 error of every method, per scenario and n: the paper's error table), `method_errors.csv` (every method on every dataset), `paired_differences_vs_GP.csv`, `dataset_summary.csv`, `computation_table.csv`, `hyperparameter_descriptive.csv`, and `unresolved_technical_failures.csv` if any job failed.
- Per dataset (`<scenario>/n<NNNN>/r<RR>/`): `summary.csv` (errors, MC noise, flags), `points.csv` (truth and all methods on the reporting grid), `errors.csv`, `diagnostics.csv`, `parameters.csv`, `chain.csv`, `timings.csv`, `concentration.csv`, `kernel_scores.csv` (CV and oracle scores over the bandwidth grid), `compact_fit.rds` (all of the above in one object), `fit.pdf/png` and `trace.pdf/png`, `verification.json` and `final_status.json`. The `Work/` subfolder holds the worker logs; the full chain and the checkpoint are deleted once the compact record has been written and verified.

### `Long_Fits/`

```
Long_Fits/
├── Tables/                 long_fit_results.csv (one row per dataset), long_fit_table.csv
└── 1d_skn/ 2d_hole/ → n<NNNN>/ → r01/ ... r05/
```

Per dataset: `config.rds`, `chains.csv`, `diagnostics.csv` (R-hat, ESS of rho*, ell, the mean intensity and rho at probe points), `final_status.csv`, `monitored_draws.rds` (draws of the monitored quantities, used by `report_hyperparameters.R`), `pooled_fit.rds` (posterior mean and pointwise 95% band, used by the validation and the 1D figure), and the chain logs. With `--compact` the full chains and checkpoints are deleted once a dataset passes.

### `Sensitivity/`

```
Sensitivity/
├── config.rds, job_log.csv, Source_Snapshot/
├── Tables/                 sensitivity_summary.csv, sensitivity_by_dataset.csv
└── S1/ ... S6/ → <setting>/ → <scenario>/ → n<NNNN>/ → r<RR>/
```

`Tables/sensitivity_summary.csv` gives, for every setting, scenario and n, the median [min, max] over the five datasets of the relative change in error against the baseline, the change in posterior mean against the Monte Carlo envelope, and the posterior means of rho* and ell (the Supplement S4.5 tables). Per fit: `summary.csv`, `points.csv`, `parameters.csv`, `result.rds` and `final_status.json`.

### `Paper_Tables/` and `Paper_Figures/`

`Paper_Tables/`: `hyperparameter_results.csv` (per dataset), `hyperparameter_table.csv` (medians over datasets) and `hyperparameter_table_rows.tex` (rows for the LaTeX table). `Paper_Figures/`: PNG previews `sim1d_fits.png` and `sim2d_maps.png` of the figures written to the manuscript folder.
