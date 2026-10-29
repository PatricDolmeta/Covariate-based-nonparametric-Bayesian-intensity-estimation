# Code structure and relation to the first-submission code

`R/` is the master code of the simulation study and of the data analysis: it is the only place where the methods are defined (sampler, data generation, competitors, diagnostics, Canada pre-processing). Every script in `Simulation Study/` and `Application/` loads its functions from here with `source("../R/...")`, so a change in `R/` reaches every run. The same `gp_sampler.R` runs both the simulations and the Canada analysis, and `competitors.R` computes the kernel curves in both, so the CV kernel estimator is identical in the two parts.

`R/` is a rewrite for the revision, not a cleaned-up version of the first-submission scripts. The old code no longer runs anything; it is kept for reference in `archive.zip` (`old_code_first_submission/` and `canada_original/`).

The revision changed the prior (ARD length-scale hyper-prior, truncated rho*) and the algorithm (pCN on a discretised GP), so the sampler is new code. The data-generation and pre-processing logic was carried over from the old scripts and re-implemented.

## Files in `R/`

| File | Role | Relation to the first-submission code |
|---|---|---|
| `gp_sampler.R` | The MCMC sampler (pCN for w, ell updates, rho*, binning, Kronecker trick), shared by the simulations and the application | **New.** Replaces `MCMC_LatentVar.R` and `Adaptive_Anisotropic_GP.R`. Written for the simulations in the revision, later extended to 3D for Canada. |
| `data_generation.R` | Covariate fields, true intensities, replicated point patterns | **Re-implementation** of the data generation in `CODE_DataGeneration_Analysis_R.R`, with the same covariate GPs (length-scales 0.005 and 0.05). Adds reproducibility fixes: Cholesky instead of eigendecomposition, exact 16-bit coding, each dataset saved once and reloaded. |
| `competitors.R` | Averaged kernel estimator (oracle and CV bandwidth) and log-polynomial Poisson models of degree 1-4 | The kernel part keeps the **idea** of the old code (average of the per-replicate ratio estimators, oracle and CV bandwidths), re-implemented with binned Gaussian smoothing on the extended bandwidth grid. The log-polynomial models are **new**. |
| `chain_diagnostics.R` | R-hat and ESS through the `posterior` package | **New.** |
| `replicate_study.R` | Engine of the replicated simulation study: per dataset, competitors -> GP chain -> summaries -> compact export | **New.** |
| `canada_data.R` | Loading and pre-processing of the wildfire data | Pre-processing **taken from** `CODE_Canada_analysis.R`: ECDF rescaling to [0,1], pixel quadrature, axis-label coordinates. Adds one rule: events on NA border pixels go to the nearest valid pixel instead of being dropped. |
| `canada_posterior.R` | Diagnostics and posterior summaries (curves, slices, maps) from the stored draws | **New.** |

## Required packages

`Matrix`, `coda`, `sn`, `jsonlite`, `posterior`, `processx`, `ps`, plus `ggplot2` and `patchwork` for the application figures. Install line and what each package is for: `../README-DATA-ANALYSIS.md`, section Software. Each runnable script checks for the packages it needs at startup.

## Which script uses which file

| Script | Loads from `R/` |
|---|---|
| `Simulation Study/run_replicate_study.R` | data_generation, gp_sampler, competitors, chain_diagnostics, replicate_study |
| `Simulation Study/run_long_fits.R` | data_generation, gp_sampler, chain_diagnostics, replicate_study |
| `Simulation Study/run_sensitivity.R` | data_generation, gp_sampler, replicate_study |
| `Simulation Study/validate_replicate_study_cap.R` | data_generation |
| `Application/run_canada.R` | gp_sampler, canada_data, canada_posterior |
| `Application/plot_canada.R` | gp_sampler, canada_data, canada_posterior, competitors |
| `Application/profile_canada.R` | gp_sampler, canada_data |

`Simulation Study/plot_paper_figures.R` and `Simulation Study/report_hyperparameters.R` load nothing from `R/`; they only read finished results.

## Other folders

The other folders hold no algorithms of their own:

* `Simulation Study/` and `Application/` contain the scripts you actually run. They handle options, worker processes, outputs and figures.
* `SIMULATION_PROTOCOL.md` and `SENSITIVITY_SPEC.md` (inside `Simulation Study/`) document the design of the simulation study.
* `Results/` holds the outputs. Each run also stores a `Source_Snapshot/` of the code it used.

Intermediate versions written earlier in the revision (`simulation_R/`, `revised_code_early_copy/`, pilot scripts) are in `archive.zip`; everything still useful from them has been merged into `R/`.
