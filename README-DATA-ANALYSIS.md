# Data Analysis: code and results

Code for the simulation study (paper Section 3, Supplement Section S4) and the wildfire application (paper Section 4, Supplement Section S5) of *A nonparametric Bayesian analysis of i.i.d. observations of covariate-driven Poisson processes*. Consolidated 12 October 2026.

```
Data Analysis/
├── R/                  shared library (one copy of each component)
│   ├── gp_sampler.R        prior, Metropolis-within-Gibbs sampler (1, 2 or 3 covariates)
│   ├── data_generation.R   covariate fields, ground truths, saved replicated datasets
│   ├── competitors.R       kernel estimator (oracle / CV bandwidth), log-polynomial models
│   ├── chain_diagnostics.R multi-chain summaries (long fits)
│   ├── replicate_study.R   one-chain fits and tables of the replicate study (simulation)
│   ├── canada_data.R       wildfire data loading and pre-processing
│   ├── canada_posterior.R  posterior summaries and convergence assessment (application)
│   └── README-R.md         what each file does, which script uses it, relation to the old code
├── Simulation Study/   simulation study: runners, protocol, results
├── Application/        wildfire application: runners, data, results
└── archive.zip         earlier code, pilot scripts, dated run records, verbatim copies of the
                        files before consolidation (paths 'archive/...' in the docs refer to its contents)
```

## Model and prior (both parts)

`rho(z) = rho* sigmoid(w(z))`, Gaussian-process prior on w with covariance `exp(-sum_j ell_j^2 (z_j - z'_j)^2)`, length-scale hyper-prior `prod_j ell_j^(a-1) e^(-b ell_j) * exp(-lambda prod_j (1 + ell_j))` with **(a, b, lambda) = (1, 2.5, 0.1)**, and `rho* ~ Gamma(1, 0.1)` truncated to [0, 100 + log n]. The likelihood integral uses pixel quadrature; covariate values are binned; w is discretised on a grid and interpolated.

## Simulation study (`Simulation Study/`)

Specification: `Simulation Study/SIMULATION_PROTOCOL.md` (version 2.2); sensitivity analysis: `Simulation Study/SENSITIVITY_SPEC.md` (version 3). Run from `Data Analysis/Simulation Study`:

```sh
Rscript run_replicate_study.R --smoke                  # 1. workflow check
Rscript run_replicate_study.R --cells=all --reps=1:5 --workers=4
Rscript run_long_fits.R --scenario=1d_skn --n=50 --reps=1:5 --workers=4 --compact     # all scenarios and n
Rscript validate_replicate_study_cap.R
Rscript run_replicate_study.R --cells=all --reps=6:50 --workers=4
Rscript run_sensitivity.R --blocks=all --workers=4
Rscript report_hyperparameters.R
Rscript plot_paper_figures.R
```

From the RStudio console (working directory `Simulation Study/`), set `REPSTUDY_ARGS`, `LONG_ARGS` or `SENS_ARGS` to the same options and `source()` the script. Saved datasets are in `Simulation Study/Data/` (generated once, checksummed; regenerated identically from their seeds if deleted).

Results (to be produced): `Results/Replicate_Study/`, `Results/Long_Fits/`, `Results/Sensitivity/`, and the paper tables and figure previews in `Results/Paper_Tables/` and `Results/Paper_Figures/`. The results with the old prior (first-round manuscript) are in `archive.zip` (`archive/simulation_results_old_prior/`).

## Application (`Application/`)

See `Application/README-APPLICATION.md`. Run from `Data Analysis/Application`: `Rscript run_canada.R --workers=4`, then `Rscript plot_canada.R`.

## Software

R 4.6.1 (tested on macOS with Apple Accelerate BLAS). Install the required packages once:

```r
install.packages(c("Matrix", "coda", "sn", "jsonlite", "posterior", "processx", "ps", "ggplot2", "patchwork"))
```

| Package | Used for |
|---|---|
| `Matrix` | sparse interpolation and binning matrices (sampler, kernel estimator) |
| `coda` | single-chain effective sample sizes |
| `posterior` | multi-chain R-hat and ESS (long fits, sensitivity analysis, application) |
| `sn` | skew-normal ground truth of the 1D simulation scenario |
| `jsonlite` | run records and verification files of the simulation study |
| `processx`, `ps` | starting and monitoring the worker processes |
| `ggplot2`, `patchwork` | application figures (`Application/plot_canada.R`) |

Each runnable script checks at startup for the packages it needs and, if any is missing, stops with the `install.packages(...)` command for the missing ones. Optional, not used in the paper: `RTriangle` and `geometry` (triangular-mesh option of the sampler), `spatstat.geom` (`keep_im = TRUE` in the data generation).

Each runner starts one R process per chain or dataset with single-threaded BLAS; `--workers` sets how many run at once.

The code also runs on Windows and Linux. R on Windows uses the reference BLAS by default, so runs are slower than on macOS. Draws may differ from macOS runs in the last digits (the saved datasets are identical), and a run cannot be resumed on a different platform from the one that started it.
