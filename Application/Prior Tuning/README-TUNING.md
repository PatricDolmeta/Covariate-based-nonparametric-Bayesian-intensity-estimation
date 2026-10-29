# Prior tuning on Ontario (11 October 2026)

Comparison of length-scale hyper-priors (a_ell, b_ell, lambda) on the three Ontario 1D fits (temperature, precipitation, wind), which led to the final prior **(1, 2.5, 0.1)** used in `../Results/Canada`. Motivation: under the earlier prior (2, 0.5, 0.05) the data pulled ell to 11-28, the curves followed year-specific structure, and the 2D and 3D fits showed competing modes.

All runs: 2 chains, 60000 sweeps, 20000 burn-in, thinning 10, rho* ~ Gamma(1, 0.1) truncated to [0, 100 + log n]. Most runs were capped at 60000 sweeps (`--max_sweeps=60000`), so a failed convergence check (✗: R-hat >= 1.01, ESS < 200 or MC noise >= 1%) partly reflects the cap; the runs with more sweeps were allowed to extend. The contents of each folder are described at the end of this file.

## 1D fits: posterior median of ell and max R-hat

| Folder | (a_ell, b_ell, lambda) | Sweeps | Temperature | Precipitation | Wind |
|---|---|---|---|---|---|
| `Ontario_a2_b2` | (2, 2, 0.05) | 60000–260000 | 1.68; R-hat 1.007 ✗ | 2.088; R-hat 1.005 | 2.939; R-hat 1.003 |
| `Ontario_a3_b3` | (3, 3, 0.05) | 60000 | 1.2; R-hat 1.004 ✗ | 1.57; R-hat 1.010 ✗ | 2.737; R-hat 1.003 |
| `Ontario_a4_b4_v2` | (4, 4, 0.05) | 60000 | 1.072; R-hat 1.017 ✗ | 1.337; R-hat 1.023 ✗ | 2.559; R-hat 1.011 ✗ |
| `Ontario_a2.5_b2.5` | (2.5, 2.5, 0.05) | 60000 | 1.153; R-hat 1.008 ✗ | 1.984; R-hat 1.075 ✗ | 2.82; R-hat 1.005 |
| `Ontario_a1_b1` | (1, 1, 0.05) | 60000 | 9.485; R-hat 1.002 | 4.268; R-hat 1.048 ✗ | 3.293; R-hat 1.001 |
| `Ontario_a1_b1.5` | (1, 1.5, 0.05) | 60000 | 2.829; R-hat 1.448 ✗ | 1.917; R-hat 1.059 ✗ | 2.994; R-hat 1.002 |
| `Ontario_a1_b2` | (1, 2, 0.05) | 60000 | 1.044; R-hat 1.025 ✗ | 1.517; R-hat 1.013 ✗ | 2.759; R-hat 1.007 |
| `Ontario_a1_b2_l0.1` | (1, 2, 0.1) | 60000 | 2.738; R-hat 1.203 ✗ | 1.386; R-hat 1.022 ✗ | 2.793; R-hat 1.005 |
| `Ontario_a1_b2.5_l0.1` | (1, 2.5, 0.1) | 60000–210000 | 0.907; R-hat 1.008 | 1.231; R-hat 1.002 | 2.626; R-hat 1.009 |

## Outcome

- The smoothing is mainly controlled by the right tail of the ell prior, i.e. by the rate b_ell: a larger b_ell pulls ell down and gives smoother curves.
- With heavier right tails (a_ell = 1, b_ell <= 2) the temperature fit reaches large ell: a wiggly curve following year-specific structure at (1, 1, 0.05) (ell about 9.5), and bimodal chains at (1, 1.5, 0.05) and (1, 2, 0.1) (R-hat 1.45 and 1.20).
- (1, 2.5, 0.1) gives unimodal, well-mixing fits for all three covariates (all checks passed once the chains were extended) without flattening the wind peak. It was adopted as the final prior.
- `Ontario_a1_b2.5_l0.1_2d3d`: capped (60000 sweeps) check of the final prior on the 2D (temperature, precipitation) and 3D fits, posterior medians of ell (0.7835, 1.012) (2D) and (0.8263, 0.9696, 1.683) (3D); mixing much improved over the earlier prior, the full runs extend until convergence.

To rerun one configuration, from `Data Analysis/Application`:

```sh
Rscript run_canada.R --fits=ontario:temp,ontario:prec,ontario:wind --a_ell=1 --b_ell=2.5 --lambda=0.1 --max_sweeps=60000 --out="Prior Tuning/Ontario_a1_b2.5_l0.1_rerun"
```

## Contents of each `Ontario_*` folder

Each folder is the `--out` directory of one `run_canada.R` call, followed by `plot_canada.R`. The chain files (`chain<k>.rds` with the stored draws, and the checkpoints) were deleted to save space, so the fits cannot be resumed or re-plotted; everything else is as written by the scripts.

```
Ontario_<prior>/
├── config.rds              all settings of the run
├── job_log.csv             one row per chain job
├── Tables/
│   └── computation_summary.csv
├── Figures/                PDF figures
└── Fits/
    └── Ontario_<model>/    one folder per fit (temp, prec, wind; temp_prec and full in _2d3d)
        ├── status.rds
        ├── chain1_target.txt, chain2_target.txt
        ├── chains.csv
        ├── diagnostics.csv
        ├── traces.pdf
        └── Logs/
```

- **`config.rds`**: R list with every setting of the run (chains, sweeps, burn-in, thinning, extension rule, GP nodes and bins, convergence thresholds, seed base), the prior (a_ell, b_ell, lambda and the rho* prior), MD5 checksums of the source files and data files (a resumed run stops if any of them changed), creation time, R version and BLAS.
- **`job_log.csv`**: one row per chain job launched by `run_canada.R`: fit id, chain, target number of sweeps, status (`complete` or `FAILED`), start and end times, wall-clock minutes. A fit that was extended appears once per extension (target 60000, then 110000, ...).
- **`Tables/computation_summary.csv`**: one row per fit (written by `run_canada.R --report`): GP nodes and covariate bins, number of years and events, sweeps reached, minutes per chain and seconds per 1000 sweeps, acceptance rates of the w and ell updates, max R-hat, min bulk and tail ESS, relative MC noise of the posterior mean of rho, posterior median and 95% interval of rho* and of each ell_j, and whether the convergence criteria were passed (`passed`) and the fit is finished (`final`). This is the table summarised above.
- **`Figures/`**: the PDF figures of `plot_canada.R` for the fits present, named as the corresponding manuscript figures. 1D folders: `ontario_1d_smoth.pdf` (posterior mean and 95% bands of rho for temperature, precipitation and wind, with the CV kernel estimate dashed) and `ontario_temp.pdf` (spatial maps of the fitted intensity rho(Z(x)) under the temperature model for 2013, 2015, 2021). `_2d3d` folder: `ontario_2d_estim.pdf` (posterior mean surface of the (temperature, precipitation) model, with the kernel estimate), `ontario_2d_spat.pdf` (its spatial maps), `2d_slicing.pdf` ((temperature, precipitation) slices of the 3D fit at the 5%, 50% and 95% quantiles of wind), `ontario_3d.pdf` (spatial maps of the 3D fit for 2013, 2015, 2021), `ontario_3d_add.pdf` (the same maps for six further years) and `ontario_3d_1.pdf`, `ontario_3d_2.pdf`, `ontario_3d_3.pdf` (all the other years of 2004-2022).
- **`Fits/Ontario_<model>/`**, one folder per fit:
  - `status.rds`: state of the fit used by `run_canada.R` to resume and extend (current target, `final`, `passed`, and the summary row of `computation_summary.csv`);
  - `chain<k>_target.txt`: number of sweeps completed by chain k;
  - `chains.csv`: one row per chain: seed, sweeps, run time, acceptance rates, final pCN step size `beta` and ell proposal scale `s_ell`;
  - `diagnostics.csv`: one row per monitored quantity (rho*, each ell_j, the log-likelihood, and rho at a fixed set of probe points, used for the MC-noise criterion): posterior mean, sd, R-hat, bulk and tail ESS, MC standard error, 2.5% and 97.5% quantiles, and whether it passes the criteria;
  - `traces.pdf`: trace plots of rho*, ell and the log-likelihood for both chains;
  - `Logs/`: console output of each chain process (`chain<k>_<date>T<time>.log`: settings, progress every 1000 sweeps, checkpoints) and its error stream (`..._stderr.log`, normally empty), one pair per job, so a fit extended three times has four pairs per chain.

`run_canada.R` and `plot_canada.R` write these capitalised folder names since 12 October 2026 (these runs originally used `fits/ontario__<model>`, `tables/`, `figures/`, `logs/`, renamed afterwards). The `fit_id` column of the CSV files keeps the internal name `ontario__<model>`, which is also the name used in the `--fits=` option.
