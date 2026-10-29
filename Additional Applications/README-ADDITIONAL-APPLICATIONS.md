# Additional applications: Saskatchewan and British Columbia (Supplement S5)

This folder runs the wildfire analysis for the two provinces of Supplement S5 with the revised prior and sampler. The Ontario analysis of the paper (Section 4) is in `../Application/` and is left untouched.

`run_additional_regions.R` is a launcher. It executes the unchanged application code (`../Application/run_canada.R`, `../Application/plot_canada.R`, and through them `../R/`) with working directory `../Application`, which is where the data files are. It points `--out` to `Results/` in this folder.

- **Nothing is written to `../Application`.** The launcher stops if `--out` would lie inside it.
- **Same configuration as Ontario.** Data pre-processing, prior (a_ℓ, b_ℓ, λ) = (1, 2.5, 0.1) with ρ\* ~ Gamma(1, 0.1) truncated at 100 + log n, GP grids (200 nodes, 20³), binning (2000 bins in 1D, 50³ in 3D), seeds and convergence rules.
- **Convergence rules:** 2 chains, 60,000 sweeps, extended by 50,000 up to 260,000 until R-hat < 1.01, bulk and tail ESS ≥ 200 and MC noise < 1%.
- **Integrity checks:** each run records MD5 checksums of the code and data files and refuses to continue if they change.

## Fits

These are the eight fits of Supplement S5:

| Province | Fits |
|---|---|
| Saskatchewan (14,300 pixels, 2,610 events) | temperature, precipitation, wind speed (1D); full model (3D) |
| British Columbia (20,795 pixels, 994 events) | temperature, precipitation, wind speed (1D); full model (3D) |

Other fits can be selected with `--fits`. For example, `--fits=saskatchewan:temp_prec,british_columbia:temp_prec` adds the 2D (temperature, precipitation) model, which the Supplement shows for Ontario only.

## How to run (R console or RStudio)

```r
setwd("/Users/matteogiordano/Desktop/Research/Preprints/(2026, arXiv) Dolmeta, Giordano/Data Analysis/Additional Applications")
rs <- function(...) system2("caffeinate", c("-i", file.path(R.home("bin"), "Rscript"), c(...)))

rs("run_additional_regions.R", "--smoke")             # short end-to-end check, a few minutes
rs("run_additional_regions.R", "--smoke", "--plot")   # figures of the smoke run
rs("run_additional_regions.R", "--workers=3")         # the 8 fits
rs("run_additional_regions.R", "--report")            # computation and diagnostics table
rs("run_additional_regions.R", "--plot")              # figures
```

Interrupting is safe: rerunning the same command resumes from the checkpoints, and finished fits are skipped.

**Time:** about 1.5–2 hours with 3 workers, based on the Ontario timings. Each 1D chain takes about 1.5 minutes per 60,000 sweeps, with a few extensions. Each 3D chain takes about 30 minutes per 60,000 sweeps; the Ontario 3D fit needed one extension, about 55 minutes per chain in total. The four 3D chains dominate, so `--workers=4` brings it closer to 1 hour.

## Outputs (`Results/Canada_Regions/`)

The layout is the same as `../Application/Results/Canada/` (see `../Application/Results/README-RESULTS.md`): `config.rds`, `job_log.csv`, `Tables/computation_summary.csv`, and `Fits/<Province>_<model>/` with the chains, diagnostics and `traces.pdf`.

`Figures/` uses the file names of the manuscript's `figures/` folder, so the files can replace the old ones directly once checked:

| File | Content | Supplement |
|---|---|---|
| `saska_1d.pdf` | Saskatchewan: posterior curves for temperature, precipitation and wind speed, with the CV kernel estimates | Fig. S1DSaska |
| `saska_3d.pdf` | Saskatchewan: spatial plug-in posterior means from the 3D fit, 2006, 2008, 2010, 2013, 2015, 2016, 2018, 2021, 2022 | Fig. S3DsasSpatial |
| `bt_columbia_1d.pdf` | British Columbia: 1D posterior curves with the CV kernel estimates | Fig. S1DBTColumbia |
| `BT_3d.pdf` | British Columbia: 3D spatial maps for the same nine years | Fig. S3DBTSpatial3 |

The covariate illustrations (`cov_saska_13/21.pdf`, `cov_bt_13/21.pdf`) show only data, which are unchanged, so they are not redrawn here.
