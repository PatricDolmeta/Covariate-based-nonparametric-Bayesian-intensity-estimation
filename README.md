# Covariate-based nonparametric Bayesian intensity estimation

Code and results for

> P. Dolmeta and M. Giordano (2026). *A Bayesian analysis of covariate-driven Poisson processes.* [journal / arXiv link]

The code of the first submission is available under the tag
[`first-submission`](../../tree/first-submission).

## Contents

| Folder | Contents | Paper |
|---|---|---|
| `R/` | Shared code: data generation, GP sampler (pCN + non-centred length-scale updates), competitors (kernel, log-polynomial), MCMC diagnostics, wildfire data and posterior summaries. See `R/README-R.md`. | all |
| `Simulation Study/` | Main simulation study, long two-chain validation fits, sensitivity analysis, paper tables and figures. See `README-SIMULATIONS.md`, `SIMULATION_PROTOCOL.md`, `SENSITIVITY_SPEC.md`. | Sec. 3; Suppl. S4.1–S4.3 |
| `Additional Simulations/` | Further ground truths, deterministic covariates, over-parametrisation. See `README-ADDITIONAL-SIMULATIONS.md`. | Suppl. S4.4–S4.6 |
| `Application/` | Ontario wildfire application, data and prior tuning. See `README-APPLICATION.md`. | Sec. 4 |
| `Additional Applications/` | Saskatchewan and British Columbia. See `README-ADDITIONAL-APPLICATIONS.md`. | Suppl. S5 |

## Paper items and the scripts producing them

| Item | Script (run from its folder) | Output |
|---|---|---|
| Table 1 (Tab:Errors) | `run_replicate_study.R`, then `make_error_table_tex.R` | `Results/Paper_Tables/` |
| Table 2 (Tab:Hyper) | `run_long_fits.R`, then `report_hyperparameters.R` | `Results/Paper_Tables/` |
| Figures `sim1d_fits`, `sim2d_maps` | `plot_paper_figures.R` | `Results/Paper_Figures/` |
| Tab:Validation (Suppl.) | `validate_replicate_study_cap.R`, then `make_validation_table_tex.R` | `Results/Paper_Tables/` |
| Trace plots (Suppl.) | `run_long_fits.R` | `Results/Long_Fits/<setting>/n<NNNN>/r<RR>/traces.pdf` |
| Tab:Sensitivity (Suppl.) | `run_sensitivity.R`, then `make_sensitivity_table_tex.R` | `Results/Paper_Tables/` |
| Suppl. S4.4–S4.6 tables and figures | `run_additional_simulations.R`, then `make_supplement_tables_tex.R`, `plot_supplement_figures.R` | `Results/Additional_Simulations/Paper_Tables/`, `…/Paper_Figures/` |
| Wildfire application (Sec. 4) | `run_canada.R`, then `plot_canada.R` | `Application/Results/` |
| Suppl. S5 | `run_additional_regions.R`, then `plot_2d_regions.R`, `plot_spatial_maps.R` | `Additional Applications/Results/` |

Exact commands, options (`--cells`, `--reps`, `--workers`, `--report`, `--smoke`) and run times are in each folder's README.

## What is included

* All code and documentation.
* Final tables (`Tables/`, `Paper_Tables/`), status reports, validation results and session information (`sessionInfo_*.txt`).
* Small per-dataset summaries (`summary.csv`, `errors.csv`, `parameters.csv`, `chain.csv`, `diagnostics.csv`, `kernel_scores.csv`, `timings.csv`), from which all tables and the numbers quoted in the text can be recomputed without sampling.
* The wildfire data (`Application/Data/`, see `README-DATA.md` for sources).

Not included, because of their size: the simulated datasets (about 6 GB), the MCMC chains and checkpoints. The datasets are regenerated exactly from the seeds in the code (`R/data_generation.R`), and every runner checks them against stored checksums. [Optional: they are also available at Zenodo, DOI: …]

## Software

R 4.6.1 with `sn` 2.1.3, `posterior`, `processx` and the packages listed in `README-DATA-ANALYSIS.md`. The dataset checksums require R ≥ 4.6. The results in the paper were computed on macOS (Apple silicon) and on Linux (Ubuntu 26.04 under WSL2, OpenBLAS 0.3.32). Results on another platform or BLAS agree within Monte Carlo error, not to the last digit.

## Licence

[keep the repository's existing licence, or e.g. MIT]