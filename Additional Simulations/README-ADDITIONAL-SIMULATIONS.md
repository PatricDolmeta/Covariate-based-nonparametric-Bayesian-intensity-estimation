# Additional simulation experiments (Supplement S4.4–S4.6)

This folder re-runs the additional experiments of Supplement Sections S4.4–S4.6 with the revised method. Only the data-generating settings are new. Everything else comes unchanged from `../R/`, exactly as in the main replicate study (`../Simulation Study`):

- the prior (a_ℓ, b_ℓ, λ) = (1, 2.5, 0.1), with ρ\* ~ Gamma(1, 0.1) truncated at 100 + log n;
- one GP chain per dataset (60,000 sweeps, 20,000 burn-in, thinning 10), with V = 200 nodes in 1D and a 20×20 grid in 2D;
- the competitors: the kernel estimator with oracle and CV bandwidths, and full log-polynomials of degree 1–4;
- the relative L² errors on the common reporting grid, flags and the verified compact export.

Nothing in `../R/` or `../Simulation Study/` is modified.

## Files

| File | Role |
|---|---|
| `additional_scenarios.R` | The new ground truths, the data loaders (including the deterministic covariate and the over-parametrised refit of the main-study data), and the hook `repstudy_load_data()` that plugs them into the engine of `../R/replicate_study.R` |
| `run_additional_simulations.R` | The runner, with the same controller, workers, checkpoints and resume logic as `run_replicate_study.R`. It writes the main study's tables plus `ADDITIONAL_SUMMARY.md` (supplement-style tables) |
| `plot_additional_figures.R` | Preview figures analogous to those of S4.4–S4.6, from the saved results (no sampling) |

## Scenarios

The window is W = [0,1]² (the paper's [−1/2,1/2]², shifted). Covariates are the GP fields of the main study (length-scales 0.005 and 0.05) on 50×50 pixels. Default sample sizes are n = 50, 250, 500, 1000, as in the main study. The ‖ρ₀‖ column is the L² norm on [0,1]^d, computed from the formulas; `Tables/truth_norms.csv` recomputes it.

| Scenario | Supplement | Truth | ‖ρ₀‖ |
|---|---|---|---|
| `1d_negexp` | S4.4, (1DTruth2) | 2 exp(3(1−z)−1) | 6.03 |
| `1d_minmax` | S4.4, (1DTruth3) | 2 + 2P(z; 3/4, 3/8) − 2P(z; 1/4, 3/8), smootherstep 6t⁵−15t⁴+10t³ | 2.27 |
| `2d_slope` | S4.5, (2DTruth3) | max{0, 30 − 90 f_SN(z; (.3,.3), .5I, (−1,−1))} | 17.86 |
| `2d_twobump` | S4.5, (2DTruth2) | 6 f_SN(z; (.3,.8), .03I, (−1,−1)) + 14 f_SN(z; (.7,.2), .05I, (3,−2)) | 22.60 |
| `2d_detcov` | S4.6, (2DTruthCoord) | 15 z₂ f_SN(z₁; .8, .3, −5), with Z₁ random and Z₂(x) = x₁ deterministic | 11.12 |
| `2d_overparam` | S4.6, over-parametrisation | 5 f_SN(z₁; .8, .3, −5): the main study's `1d_skn` datasets (same observations), refitted with the extra covariate Z₂(x) = x₁ | 6.42 |

Notes on the truths:

- **Source of the formulas.** The truths implement the formulas displayed in the Supplement, not the first-submission code. The old code used extra scale factors (×10 for the slope truth; different weights, centres and τ = 1 for the two-bump truth).
- **Smootherstep typo.** The "0.6t⁵" in the Supplement's definition of p(t) is a typo for 6t⁵.
- **Norms in the Supplement captions.** The norms stated in the captions of Tables S1DTab2, S2DResults3 and S2DResultsCoord do not all match these formulas; use `Tables/truth_norms.csv`.
- **Deterministic covariate.** Z₂(x) = x₁ is evaluated at the pixel centres and 16-bit coded like the random fields.
- **Over-parametrisation.** `2d_overparam` reads the datasets from `../Simulation Study/Data`, so the over-parametrised fit can be compared with the main study's 1D fit of the same observations. The comparison, paired by `data_id`, is in `Tables/overparam_vs_main_1d.csv`. A dataset not yet present there is simulated from the same seed into `Data/` here.

Seeds: data seed = scenario base + 1000 n + rep. The bases are 11e6, 13e6, 21e6, 23e6 and 25e6, separate from the main study's seeds. `2d_overparam` uses the main study's `1d_skn` seeds. MCMC seed = data seed + 70000007.

## How to run (R console or RStudio)

```r
setwd("<path to the repository>/Additional Simulations")
rs <- function(...) system2("caffeinate", c("-i", file.path(R.home("bin"), "Rscript"), c(...)))

rs("run_additional_simulations.R", "--smoke")                                    # workflow check, a few minutes
rs("run_additional_simulations.R", "--cells=all", "--reps=1:5", "--workers=3")   # stage 1
rs("run_additional_simulations.R", "--cells=all", "--reps=6:50", "--workers=3")  # stage 2
rs("run_additional_simulations.R", "--report")                                   # rebuild tables, no sampling
rs("plot_additional_figures.R")                                                  # figures, no sampling
```

- `--cells` selects scenarios or single cells, e.g. `--cells=1d_negexp,1d_minmax` (all four n) or `--cells=2d_slope:250,2d_detcov:1000`. Any n in 2..1000 can be given explicitly, e.g. `--cells=2d_slope:10`.
- Interrupting is safe: rerun the same command and completed datasets are skipped.
- After any change to the scripts here or to `../R/`, use a new `--out` folder (the run checks the code checksums).

**Time** (M1, 3 workers, from the main-study timings): about 1.8 min per 1D fit and 10–14 min per 2D fit including the competitors. That makes about 70 minutes per dataset index over all 24 cells: about 6 hours for stage 1 and 2–2.5 days for all 50 datasets. The 1D scenarios alone (`--cells=1d_negexp,1d_minmax`) take about 4 hours for 50 datasets.

## Outputs (`Results/Additional_Simulations/`)

The layout is the same as `../Simulation Study/Results/Replicate_Study/`: `config.rds`, `job_log.csv`, `STATUS.md`, `Source_Snapshot/`, and per-dataset folders `<scenario>/n<NNNN>/r<RR>/`, which hold `points.csv` (truth, all methods, posterior band), `summary.csv`, `errors.csv`, `fit.pdf`, `trace.pdf`, etc.

In addition:

- `ADDITIONAL_SUMMARY.md`: per scenario, mean (SD) of the relative L² error of every method for each n, the truth norm, and the over-parametrisation comparison.
- `Tables/truth_norms.csv` and `Tables/overparam_vs_main_1d.csv`, besides the main study's tables (`error_table.csv`, `paired_differences_vs_GP.csv`, `computation_table.csv`, …).
- `Figures/` (from `plot_additional_figures.R`): `1d_negexp.pdf`, `1d_minmax.pdf`, `2d_*_maps.pdf`, `2d_slope_diagonals.pdf`, `2d_detcov_slices.pdf`, `2d_overparam_slices.pdf`.
- `Data/` (next to this README): the saved datasets, one `.rds` per scenario, n and dataset.
