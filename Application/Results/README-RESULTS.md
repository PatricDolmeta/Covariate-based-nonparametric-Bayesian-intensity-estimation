# Results of the wildfire application

| Folder | Content | Written by |
|---|---|---|
| `Canada/` | The 13 fits of the paper with the final prior (a_ell, b_ell, lambda) = (1, 2.5, 0.1), their tables and the figures of Section 4 and Supplement S5 (to be produced) | `run_canada.R`, `plot_canada.R` |
| `Profile/` | Binning accuracy and timing profile with the final configuration (12 October 2026), which confirms the binning used (2000, 200² and 50³ bins) | `profile_canada.R` |
| `Smoke/` | Short end-to-end check (`--smoke`); can be deleted once it has passed | `run_canada.R --smoke` |

Run all scripts from `Data Analysis/Application` (see `../README-APPLICATION.md`). The prior comparison on Ontario is in `../Prior Tuning/`.

## `Canada/`

```
Canada/
├── config.rds          all settings and the prior (and MD5 checksums of the code and data files)
├── job_log.csv         one row per chain job: fit, chain, target sweeps, status, start/end times, minutes
├── Tables/
│   └── computation_summary.csv   one row per fit: timings, acceptance rates, R-hat, ESS, MC noise,
│                                 posterior median and 95% interval of rho* and ell, passed/final
├── Fits/
│   ├── Ontario_temp/  Ontario_prec/  Ontario_wind/  Ontario_temp_prec/  Ontario_full/
│   ├── Saskatchewan_temp/  Saskatchewan_prec/  Saskatchewan_wind/  Saskatchewan_full/
│   └── British_Columbia_temp/  British_Columbia_prec/  British_Columbia_wind/  British_Columbia_full/
└── Figures/            PDF figures of the paper (plot_canada.R)
```

Each folder in `Fits/` contains:

- `chain<k>.rds`: the result of chain k, including the stored draws of the GP node values and rho*. All posterior summaries and figures are computed from these, so the figures can be redrawn at any time without rerunning the chains;
- `chain<k>_checkpoint.rds`: the sampler state, used to resume after an interruption or to extend the fit;
- `chain<k>_target.txt`: number of sweeps completed by chain k;
- `status.rds`: state of the fit used by `run_canada.R` (current target, finished, criteria passed, summary row);
- `chains.csv`: per chain: seed, sweeps, run time, acceptance rates, final step sizes;
- `diagnostics.csv`: per monitored quantity (rho*, each ell_j, log-likelihood, rho at the probe points): posterior mean, sd, R-hat, bulk and tail ESS, MC standard error, 95% interval, criteria passed;
- `traces.pdf`: trace plots of rho*, ell and the log-likelihood for both chains;
- `Logs/`: console output and error stream of each chain job (`chain<k>_<date>T<time>.log`, `..._stderr.log`).

`Figures/` (file names as in the manuscript's `figures/` folder, so the new figures can replace the old ones directly):

| File | Content | Paper |
|---|---|---|
| `ontario_1d_smoth.pdf` | 1D posterior curves (temperature, precipitation, wind) with the CV kernel estimates | Fig. 5 |
| `ontario_temp.pdf` | spatial maps from the temperature fit, 2013, 2015, 2021 | Fig. 6 |
| `2d_slicing.pdf` | (temperature, precipitation) slices of the 3D fit at the 5%, 50%, 95% wind quantiles | Fig. 7 |
| `ontario_3d.pdf` | spatial maps from the 3D fit, 2013, 2015, 2021 | Fig. 8 |
| `ontario_2d_estim.pdf`, `ontario_2d_spat.pdf` | 2D (temperature, precipitation): observed covariates and events, GP surface on the whole square, kernel estimate inside the hull of the observed covariates (white line), and the spatial maps | Supplement |
| `ontario_3d_add.pdf`, `ontario_3d_1.pdf`, `ontario_3d_2.pdf`, `ontario_3d_3.pdf` | 3D maps for the other years | Supplement |
| `saska_1d.pdf`, `saska_3d.pdf` | Saskatchewan: 1D curves, 3D maps for nine years | Supplement |
| `bt_columbia_1d.pdf`, `BT_3d.pdf` | British Columbia: 1D curves, 3D maps for nine years | Supplement |

Size: the 1D and 2D fits take about 10-30 MB each; the 3D fits store 8000 node values per saved draw, a few hundred MB per chain (about 1-2 GB for the three 3D fits). If needed once the paper is final, the checkpoints can be deleted: the figures can still be redrawn from `chain<k>.rds`, but the fits can no longer be extended.

## `Profile/`

Output of `profile_canada.R` (Ontario), before the final runs:

- `data_summary.csv`: years, valid pixels, pixel area, events (2 to 130 per year), border events moved to the nearest valid pixel;
- `binning_accuracy.csv`: error of the binned log-likelihood against the exact one, per covariate dimension and number of bins, with the cost of one evaluation;
- `timing.csv`: a 300-sweep chain per covariate dimension, extrapolated to 60000 sweeps;
- `session.txt`: R version and BLAS.

Run with the final prior (1, 2.5, 0.1) and GP grids (200 nodes, 20×20, 20³). Over 20 prior draws, the binning error of the log-likelihood has SD 0.001 nats in 1D (2000 bins), 0.02 in 2D (200²) and 0.07 in 3D (50³; largest difference between draws 0.32), far below the 1-nat scale that matters for the Metropolis-Hastings ratios. In 3D, 50 bins per axis are about as accurate as 100 (SD 0.07 vs 0.05) at less than half the cost (about 28 against 66 minutes per 60000-sweep chain).
