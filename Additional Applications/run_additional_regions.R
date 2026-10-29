###############################################################################
# run_additional_regions.R -- wildfire analysis for Saskatchewan and British Columbia
# (Supplement Section S5), with the revised prior and sampler.
#
# This is a launcher: it runs the UNCHANGED application code (../Application/run_canada.R and
# ../Application/plot_canada.R, which use ../R/) with working directory ../Application (where the
# data are) and an absolute --out that points into THIS folder. Nothing in ../Application is
# written or modified; the settings, prior, seeds and convergence rules are exactly those of the
# Ontario analysis of the paper (2 chains, 60000 sweeps, extended by 50000 up to 260000 until
# R-hat < 1.01, bulk/tail ESS >= 200 and MC noise < 1%).
#
# Default fits: the eight fits of Supplement S5 (temperature, precipitation, wind speed and the
# full three-covariate model, for each of the two provinces).
#
# Run from 'Data Analysis/Additional Applications':
#   Rscript run_additional_regions.R --smoke                  # short end-to-end check (Results/Smoke)
#   Rscript run_additional_regions.R --smoke --plot           # figures of the smoke run
#   Rscript run_additional_regions.R --workers=3              # the 8 fits (Results/Canada_Regions)
#   Rscript run_additional_regions.R --report                 # computation/diagnostics table
#   Rscript run_additional_regions.R --plot                   # figures (Results/Canada_Regions/Figures)
#   Rscript run_additional_regions.R --fits=saskatchewan:temp_prec --workers=2   # any other fit
# All other options of run_canada.R (--chains, --sweeps, ...) are passed through unchanged; changing
# a setting requires a new --out folder, exactly as for run_canada.R.
#
# From the RStudio console (working directory 'Data Analysis/Additional Applications'):
#   ADDAPP_ARGS <- c("--workers=3"); source("run_additional_regions.R")
###############################################################################
if (!requireNamespace("processx", quietly = TRUE))
  stop("Missing package processx. Install with install.packages(\"processx\")", call. = FALSE)

APP_DIR <- normalizePath(file.path("..", "Application"), mustWork = FALSE)
if (!file.exists(file.path(APP_DIR, "run_canada.R")) || !file.exists(file.path(APP_DIR, "plot_canada.R")))
  stop("Run this script from 'Data Analysis/Additional Applications' (../Application/run_canada.R not found)")

REGION_FITS <- paste0(rep(c("saskatchewan", "british_columbia"), each = 4), ":",
                      c("temp", "prec", "wind", "full"))
SMOKE_FITS <- c("saskatchewan:temp", "british_columbia:full")

args <- if (exists("ADDAPP_ARGS", envir = globalenv(), inherits = FALSE))
  get("ADDAPP_ARGS", envir = globalenv()) else commandArgs(TRUE)
smoke <- "--smoke" %in% args
plot <- "--plot" %in% args
report <- "--report" %in% args
out_arg <- args[startsWith(args, "--out=")]
if (length(out_arg) > 1) stop("Duplicate option --out")
out <- if (length(out_arg)) sub("^--out=", "", out_arg) else
  file.path("Results", if (smoke) "Smoke" else "Canada_Regions")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)                      # absolute: the scripts run in ../Application
if (startsWith(out, paste0(APP_DIR, .Platform$file.sep)) || identical(out, APP_DIR))
  stop("--out must not lie inside ../Application")

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run_in_app <- function(script, script_args) {
  cat(sprintf("Running %s in %s\n  output: %s\n", script, APP_DIR, out)); flush.console()
  res <- processx::run(rscript, c(script, script_args), wd = APP_DIR, echo = TRUE,
                       stderr_to_stdout = TRUE, error_on_status = FALSE, cleanup_tree = TRUE)
  if (!identical(as.integer(res$status), 0L)) stop(script, " exited with status ", res$status, call. = FALSE)
  invisible(res$status)
}

if (plot) {
  # plot_canada.R draws every figure whose fits exist in --out: here the Saskatchewan and British
  # Columbia figures (the Ontario ones are skipped with a message, as their fits are not in --out).
  run_in_app("plot_canada.R", c(paste0("--out=", out), if (smoke) "--kernel=FALSE"))
} else {
  pass <- setdiff(args[!startsWith(args, "--out=")], "--plot")
  if (!report && !any(startsWith(pass, "--fits=")))
    pass <- c(pass, paste0("--fits=", paste(if (smoke) SMOKE_FITS else REGION_FITS, collapse = ",")))
  run_in_app("run_canada.R", c(pass, paste0("--out=", out)))
}
