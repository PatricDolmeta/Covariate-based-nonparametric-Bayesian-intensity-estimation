# Sensitivity study: specification

**Version 3, 12 October 2026.** Same design as version 2 (10 October 2026), re-centred on the version-2.2 baseline: prior $(a_\ell, b_\ell, \lambda) = (1, 2.5, 0.1)$, with blocks S5 and S6 moved around it. Version 2 and its results (old prior) are preserved in `archive/pre_consolidation_copies/SENSITIVITY_SPEC_DRAFT.md` and `archive.zip` (`archive/simulation_results_old_prior/Sensitivity_v2_1/`).

## Scope

| Question | Requested by / claimed in | Blocks |
|---|---|---|
| Effect of GP discretisation, likelihood binning and spatial quadrature; rule of thumb for the quadrature | AE 1, R1.3, R2 (numerical approximations); manuscript Sec. 3.2 | S1–S3 |
| Is the n-dependent truncation of ρ\* needed; effect of a "standard" prior on ρ\* | AE 2, R1.2, R2 (ρ\* prior) | S4a–b |
| Hyper-parameters affect only constants; robustness to their values | Manuscript Sec. 2.1 | S4b, S5, S6 |

## Design

- One factor at a time around the replicate-study baseline (`Results/Replicate_Study`): same truths, sampler, budget (one chain, 60000 sweeps, 20000 burn-in, thinning 10) and seeds (`mcmc_seed = data_seed + 70000007`, common random numbers with the replicate study).
- Both scenarios, n = 50 and 1000, datasets 1–5. GP only.
- Baseline: the replicate-study fit of the same dataset (read, never rerun); S3 fits its own baseline (q = 50) on new data simulated at 200×200 pixels.

## Factors (baseline in bold)

| Block | Factor | 1D | 2D |
|---|---|---|---|
| S1 | GP discretisation | V = 100, **200**, 400 | 15×15, **20×20**, 35×35 |
| S2 | Likelihood binning | 500, **2000**, 8000 bins | 100², **200²**, 400² |
| S3 | Spatial quadrature | q = 25, **50**, 100 (data at 200×200) | same |
| S4a | Truncation of ρ\* | c = **100**, ∞ | same |
| S4b | Rate of the ρ\* prior | b_ρ = 0.01, **0.1**, 2 | same |
| S5 | Tilt | λ = 0.05, **0.1**, 0.25 | same |
| S6 | Rate of the exponential ℓ hyper-prior (a_ℓ = 1) | b_ℓ = 1.5, **2.5**, 4 | same |

S6 moves the right tail of the length-scale prior, which is what controls smoothing (prior comparison on the Ontario data, 11 October 2026): b_ℓ = 1.5 gives a heavier tail (rougher fits allowed), b_ℓ = 4 a lighter one.

S3 data: keys `sens_<scenario>_m200_n<n>_r<r>`, data seed `5000000 + 1000000*d + 1000*n + r`; events from the full 200×200 intensity; quadrature on q×q subsampled pixel centres with weight 1/q². Rule of thumb (R1.3): quadrature cell size 1/q against the correlation length of the roughest covariate (√0.005 ≈ 0.07).

## Reported

Per fit and baseline: relative change in error (err − err_base)/err_base; change in posterior mean ‖mean − mean_base‖/‖ρ0‖ against the Monte Carlo envelope √(mcse² + mcse_base²)/‖ρ0‖, with mcse the RMS Monte Carlo standard error of each fit's posterior mean over the reporting grid (same scale as the change; the two chains are treated as independent, so the envelope is a rough yardstick, not a test); posterior means of ρ\* and ℓ (S4: also the 95% interval of ρ\* and the fraction of draws near the truncation bound). Supplement S4.5: one table per scenario, median [min, max] over the five datasets.

Version-2 findings (old prior, 10–11 October 2026), for reference: discretisation on the plateau (only the 15×15 grid in 2D at n = 1000 raised the error, by 2.8%); quadrature q = 25–100 within about 1%; truncation never active (c = ∞ identical); b_ρ = 2 degrades the estimates (prior mean far below the intensity), b_ρ ≤ 0.1 harmless; λ and (a_ℓ, b_ℓ) changes of a few percent, fading with n.

## Running

From `Data Analysis/Simulation Study`, after the replicate-study run (stage 1 at least for datasets 1–5):

```sh
Rscript run_sensitivity.R --list
Rscript run_sensitivity.R --smoke
Rscript run_sensitivity.R --blocks=all --workers=4          # 280 fits, output Results/Sensitivity
Rscript run_sensitivity.R --blocks=S6                       # one block
Rscript run_sensitivity.R --report
```
