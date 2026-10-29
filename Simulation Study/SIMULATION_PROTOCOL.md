# Simulation-study protocol

**Version 2.2, 12 October 2026.** Current specification of the replicated simulation study (paper Section 3, Supplement Section S4). Version 2.1 (8 October 2026) is preserved verbatim in `archive/pre_consolidation_copies/SIMULATION_PROTOCOL.md`; its results (old prior) are in `archive.zip` (`archive/simulation_results_old_prior/Replicate_Study_v2_1/`); the version-2.2 results go to `Results/`.

**What changed in version 2.2, and why.**
- **Length-scale hyper-prior** $(a_\ell, b_\ell, \lambda)$: (2, 0.5, 0.05) → **(1, 2.5, 0.1)**, the prior chosen for the wildfire application, so that one prior is used throughout the paper. On datasets 1–5 of every cell (sensitivity check of 11 October 2026, `archive.zip` (`archive/simulation_results_old_prior/Sensitivity_appprior/`)), the change altered the GP error by −1.4% to +3.4% (median over datasets) from n = 250 on, and by +16% (2D) to +23% (1D) at n = 50. From n = 250 the GP remained the most accurate method.
- **Kernel bandwidth grid**: 1D 1000 bins and 25 bandwidths in [0.003, 1] (was 500 and 20 in [0.01, 1]); 2D 160×160 bins and 12×12 bandwidths in [0.008, 2]² (was 80×80 and 8×8 in [0.02, 1]²). In 2D the old lower edge was often selected at large n. On a check of 8 datasets (11 October 2026) no selection reached the new edge, and the kernel errors changed by a few percent in both directions.
- **Parametric competitors**: full log-polynomial Poisson models of total degree 1–4 (all interactions in 2D) are part of the replicate study. In version 2.1 they were add-ons; the additive log-cubic is dropped.
- **2D GP nodes**: 20×20 tensor grid (was 25×25), as in the application (20 nodes per axis in 2D and 3D). Not checked in advance: block S1 of the sensitivity analysis measures it (15×15 and 35×35 around the new baseline). In version 2 (old prior), 15×15 raised the 2D error at n = 1000 by 2.8% relative to 25×25, and 35×35 changed it by −0.3%.
- **Code**: one shared library in `Data Analysis/R/` (the sampler is the 3-covariate version used for the application, which reproduces the version-2.1 sampler exactly in 1D and 2D).

## 1. Study design

| Item | Specification |
|---|---|
| Scenarios | `1d_skn`, `2d_hole` |
| Sample sizes | n = 50, 250, 500, 1000 (both scenarios) |
| Repetitions | 50 independent datasets per scenario and n: 400 datasets |
| Methods | GP posterior mean; averaged ratio-form kernel with oracle and with leave-one-replicate-out CV bandwidth; log-polynomial Poisson models of degree 1–4 |
| Primary output | Relative L2 error of each method on a common reporting grid |

n is the number of independent covariate-field/point-pattern pairs in one dataset. Dataset r is an independently simulated set of n pairs.

## 2. Observations: one saved dataset for all methods

Spatial window [0,1]², m = 50 pixel centres per axis, all 2500 pixels used for quadrature. For each covariate independently:
1. `K[a,b] = exp(-(x[a]-x[b])^2/(2*ls))`, `B Bᵀ = K + j I`, j = 10⁻⁶;
2. `X = B E Bᵀ` with independent standard-normal entries of E;
3. `Z = Φ(X/(1+j))` (population standardisation);
4. Z rounded to multiples of 1/65535;
5. in each pixel, Poisson(rho0(Z)/m²) events, placed uniformly within the pixel.

| Scenario | Covariate ls | Truth rho0 |
|---|---|---|
| 1D | 0.005 | `5 * sn::dsn(z, xi=0.8, omega=0.3, alpha=-5)` |
| 2D | (0.005, 0.05) | `pmax(0, 10 - 10*g(0.8,0.3) + 10*g(0.3,0.8))`, `g(a,b) = exp(-0.5*((z1-a)^2/0.08^2 + (z2-b)^2/0.5^2))` |

Key `<scenario>_n<n>_r<r>`, data seed `1000000*d + 1000*n + r`, generated once by `get_replicates()` and saved in `Data/`. The saved record and its checksum (`data_id`) are authoritative; every method analyses the same observations. The datasets of version 2.1 are reused unchanged.

## 3. Prior and GP representation

`rho(z) = rho_star * sigmoid(w(z))`, GP covariance `exp(-sum_j ell_j^2 (z_j - z'_j)^2)` (ell = inverse length-scales), length-scale prior

`pi(ell) ∝ prod_j ell_j^(a_ell-1) exp(-b_ell*ell_j) * exp(-lambda*prod_j(1+ell_j))`, with **a_ell = 1, b_ell = 2.5, lambda = 0.1**,

and `rho_star ~ Gamma(shape 1, rate 0.1)` truncated to [0, 100 + log n].

| Numerical component | 1D | 2D |
|---|---|---|
| GP nodes | 200 equally spaced on [0,1] | 20×20 tensor grid |
| Interpolation | piecewise linear | bilinear |
| Covariance factorisation | dense | Kronecker |
| Likelihood covariate bins | 2000 | 200×200 |
| Spatial quadrature | 50×50 pixels | 50×50 pixels |
| Reporting points | 500 midpoints | 60×60 midpoints |

## 4. Sampling

Metropolis-within-Gibbs: ten pCN moves for w per sweep; truncated-gamma Gibbs step for rho_star; sequential non-centred log-random-walk updates of ell_j; Robbins–Monro adaptation during burn-in only (targets 0.30 for pCN, 0.44 for ell), frozen at the average over the second half of burn-in; cold start drawn from the prior.

### 4.1 Replicate study (paper error table): `run_replicate_study.R`

| Item | Rule |
|---|---|
| Budget | One chain per dataset: 60000 sweeps, 20000 burn-in, thinning 10. Fixed, no extensions. |
| Seed | `mcmc_seed = data_seed + 70000007` |
| Order | Stage 1: datasets 1–5 in every cell; stage 2: datasets 6–50. |
| Saved per dataset | Posterior mean and its MCSE at all reporting points; all competitor curves; errors; ESS/MCSE and traces of rho_star, ell, log-likelihood, mean intensity and rho(z) probes; selected bandwidths (with grid-edge indicator); acceptance, tuning, seeds, start values; timings; flags; provenance. |
| Inspection flags | relative RMS posterior-mean MC noise > 1%, or within-chain split R-hat > 1.05 on the rho(z) probes or the mean intensity. Flags mark datasets for inspection; they never exclude, redraw or extend a dataset. |
| Output | `Results/Replicate_Study/` (tables in `Tables/`, status in `STATUS.md`) |

### 4.2 Long two-chain fits (hyper-parameter table, diagnostics, validation): `run_long_fits.R`

Datasets 1–5 of every cell. Two chains from independent prior draws, burn-in 50000, thinning 10, initial 100000 sweeps (1D) or 200000 (2D; 250000 at n = 1000), extended by 50000 sweeps up to 400000 until rank-normalised R-hat < 1.01 and bulk and tail ESS ≥ 200 for all monitored quantities. `mcmc_seed = data_seed + 7 + 10000 * chain`. Output: `Results/Long_Fits/`.

### 4.3 Validation of the replicate-study budget: `validate_replicate_study_cap.R`

Each replicate-study fit with a long reference is compared with it: ‖replicate-study mean − reference mean‖/‖truth‖ < 1% and absolute error change < 0.002. All 40 comparisons (2 scenarios × 4 sample sizes × datasets 1–5) are required: `SUMMARY.md` opens with a verdict and lists any missing or failed comparison, and the script stops with an error unless all 40 are present and pass. Results are reported in full; a failure triggers review of the budget, never the exclusion of a dataset. (Version 2.1: 23 of 24 comparisons passed; the one exception, 1D n = 50, was within Monte Carlo noise.)

## 5. Competitors and error metric (`R/competitors.R`)

- **Kernel**: plain average over replicates of the per-replicate ratio-form estimates, Gaussian kernel, covariate values binned (1D 1000 bins, 2D 160×160); bandwidth grids 25 log-spaced values in [0.003, 1] (1D) and 12×12 in [0.008, 2]² (2D). Oracle: minimal error to the truth on the bin centres (a truth-informed benchmark). CV: maximal leave-one-replicate-out Poisson likelihood. Evaluated at the reporting points by linear (1D) or bilinear (2D) interpolation.
- **Log-polynomials**: intercept plus all monomials of total degree 1–4, exact Poisson maximum likelihood on the pixel quadrature (Newton–Raphson).
- **Error**: `sqrt(mean((estimate-truth)^2))/sqrt(mean(truth^2))` on the reporting points. Per cell: mean, between-dataset SD and simulation MCSE (SD/√50); paired differences to the GP.

## 6. Sensitivity analysis

Specification: [SENSITIVITY_SPEC.md](SENSITIVITY_SPEC.md) (version 3). Baseline: the version-2.2 replicate-study fits.

## 7. Execution

From `Data Analysis/Simulation Study` (Terminal), or from the RStudio console with `REPSTUDY_ARGS`, `LONG_ARGS`, `SENS_ARGS`:

```sh
Rscript run_replicate_study.R --smoke                          # workflow check (Results/Replicate_Study_Smoke)
Rscript run_replicate_study.R --cells=all --reps=1:5 --workers=4
Rscript run_long_fits.R --scenario=1d_skn  --n=50 --reps=1:5 --workers=4 --compact   # and n = 250, 500, 1000
Rscript run_long_fits.R --scenario=2d_hole --n=50 --reps=1:5 --workers=4 --compact   # and n = 250, 500, 1000
Rscript validate_replicate_study_cap.R                         # paired check against the long fits
Rscript run_replicate_study.R --cells=all --reps=6:50 --workers=4
Rscript run_sensitivity.R --blocks=all --workers=4
Rscript report_hyperparameters.R                          # paper Table Tab:Hyper
Rscript plot_paper_figures.R                              # paper figures sim1d_fits, sim2d_maps
```

All runners resume interrupted work and refuse to mix results from different code or settings (a new `--out` is needed after any change).
