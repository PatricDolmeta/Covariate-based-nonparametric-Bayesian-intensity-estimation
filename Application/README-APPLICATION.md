# Wildfire application (paper Section 4, Supplement Section S5)

Bayesian analysis of June wildfire occurrences in Ontario (main text), Saskatchewan and British Columbia (Supplement), 2004–2022 (n = 19 years), with average June temperature, precipitation and wind speed as covariates. Uses the shared library in `../R/` (sampler `gp_sampler.R`, data `canada_data.R`, posterior summaries `canada_posterior.R`, kernel estimator `competitors.R`).

Run every script from this folder (`Data Analysis/Application`), either with `Rscript` in the Terminal as below, or from the R console (section "How to run from the R console").

```sh
Rscript run_canada.R --smoke                       # short end-to-end check (Results/Smoke)
Rscript plot_canada.R --out=Results/Smoke --kernel=FALSE
Rscript run_canada.R --workers=4                   # the 13 fits of the paper (Results/Canada)
Rscript run_canada.R --report                      # computation/diagnostics table
Rscript plot_canada.R                              # figures (Results/Canada/Figures)
Rscript profile_canada.R                           # optional: binning accuracy and timings
```

Selected fits: `--fits=ontario:temp,ontario:prec,ontario:wind,ontario:temp_prec,ontario:full` (regions `ontario`, `saskatchewan`, `british_columbia`; models `temp`, `prec`, `wind`, `temp_prec`, `full`). Rerunning a command resumes from the checkpoints; changing any setting or source file requires a new `--out`.

## How to run from the R console

In a fresh R or RStudio session, set the working directory and define a helper `rs()` that runs a script with `Rscript` (`caffeinate -i` keeps the Mac awake during the run, as long as the lid is open):

```r
setwd("/Users/matteogiordano/Desktop/Research/Preprints/(2026, arXiv) Dolmeta, Giordano/Data Analysis/Application")
rs <- function(...) system2("caffeinate", c("-i", file.path(R.home("bin"), "Rscript"), c(...)))
```

Then run the scripts with the same options as in the Terminal, e.g. the five Ontario fits followed by the figures:

```r
rs("run_canada.R", "--fits=ontario:temp,ontario:prec,ontario:wind,ontario:temp_prec,ontario:full", "--workers=3")
rs("plot_canada.R")
```

Each `rs()` call waits until the script has finished, so several lines can be pasted at once and run in sequence; the console is busy meanwhile (progress is printed there). `--workers` is the number of chains run in parallel, one core each. Rerunning an interrupted `run_canada.R` command resumes from the checkpoints. `caffeinate` exists only on macOS: elsewhere, define `rs <- function(...) system2(file.path(R.home("bin"), "Rscript"), c(...))`.

## Data and pre-processing (`../R/canada_data.R`, data in `Data/`)

- Years 2004–2022 selected by name. Covariates mapped to [0, 1] by the empirical c.d.f. of their pixel values pooled over the 19 years (as in the original analysis; the manuscript's "c.d.f. of the N(0,10) distribution" needs correcting).
- Likelihood integral: every valid pixel of every year is a quadrature node with weight the pixel area (original pixel units, so intensities are of order 1e-3 to 1e-2 per pixel).
- Events take the covariate values of the pixel containing them; border events on NA pixels go to the nearest valid pixel (3 in Ontario, 31 in Saskatchewan, 11 in British Columbia, all within 2.3 pixel widths; the original code dropped them).
- Ontario: 21615 pixels, 648 events. Saskatchewan: 14300 pixels, 2610 events. British Columbia: 20795 pixels, 994 events.

## Model, prior and sampler

Same model and sampler as the simulation study (`../R/gp_sampler.R`), extended to three covariates (20×20×20 tensor grid, trilinear interpolation, Kronecker factorisation in each coordinate).

- Prior: `rho(z) = rho* sigmoid(w(z))`, kernel `exp(-sum_j ell_j^2 (z_j - z'_j)^2)`, length-scale hyper-prior with **(a_ell, b_ell, lambda) = (1, 2.5, 0.1)** (exponential gamma part, tilt 0.1), `rho* ~ Gamma(1, 0.1)` truncated to [0, 100 + log n].
- Why this prior: under the earlier (2, 0.5, 0.05) the data pulled ell to 11–28 and the curves followed year-specific structure (the CV kernel shows the same bumps), with competing modes in 2D and 3D. A comparison of ten priors on the Ontario 1D fits (11 October 2026; tables and figures in `Prior Tuning/`, see `Prior Tuning/README-TUNING.md`) showed that the smoothing is controlled by the right tail of the ell prior (rate b_ell); with a_ell = 1, b_ell = 2.5 is the smallest rate tried that gives unimodal, well-mixing fits for all three covariates without flattening the wind peak.
- Numerics: GP nodes 200 (1D), 20×20 (2D), 20³ (3D); covariate bins 2000, 200², 50³ (3D: 50 bins per axis, about 2.5× faster than 100 with the same binning accuracy at the fitted length-scales).
- MCMC: 2 chains per fit, 60000 sweeps, 20000 burn-in, thinning 10, extended by 50000 sweeps up to 260000 until R-hat < 1.01, bulk and tail ESS ≥ 200 and relative MC noise of the posterior mean < 1%.

## Kernel curves (`plot_canada.R`)

The dashed curves (1D) and the third 2D panel are the averaged ratio-form kernel estimator of the simulation study with leave-one-year-out likelihood-CV bandwidth (`kernel_competitor()` in `../R/competitors.R`, extended bandwidth grid). The selected bandwidth is printed for each curve, with a warning if it lies on the edge of the grid.

## Results

- `Results/Canada/` (to be produced): the 13 fits of the paper with the final prior.
- `Prior Tuning/`: prior comparison on Ontario (compact records; chain files deleted), summarised in `Prior Tuning/README-TUNING.md`.
- `Results/Profile/`: binning accuracy and timing profile of 11 October 2026 (`profile_canada.R`), on which the 3D binning (50 bins per axis) is based.
- The verification of the 3D sampler extension (11 October 2026, all checks passed) is in `archive.zip`: `archive/application/verify_nd_sampler.R` and its report `verify_20261011T063453.txt`.
