###############################################################################
# validate_replicate_study_cap.R -- paired check of the one-chain replicate-study fits (60000 sweeps,
# 20000 burn-in) against the long two-chain reference fits of the same saved datasets
# (run_long_fits.R). Every dataset with both a completed replicate-study record and a passing long
# fit is compared (both scenarios, all n).
#
# Tolerances (protocol, Section 4.3):
#   ||replicate-study mean - reference mean|| / ||truth|| < 0.01  and  |error change| < 0.002.
# A failure is reported, never used to exclude or redraw a dataset. No sampling.
# Completeness: 40 comparisons are expected (2 scenarios x 4 sample sizes x datasets 1-5).
# SUMMARY.md starts with a verdict (all passed / failed / incomplete) and lists any missing or
# failed comparison; the script stops with an error unless all 40 are present and pass, so a
# successful run means the 60000-sweep budget is validated.
#
# Usage (from 'Data Analysis/Simulation Study'):
#   Rscript validate_replicate_study_cap.R [repstudy_dir] [long_fit_dir]
#   defaults: Results/Replicate_Study  Results/Long_Fits
# Writes <repstudy_dir>/Validation/cap_validation.csv and SUMMARY.md.
###############################################################################
# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("sn")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
source("../R/data_generation.R")
args <- commandArgs(TRUE)
out <- if (length(args) >= 1) args[1] else "Results/Replicate_Study"
long <- if (length(args) >= 2) args[2] else "Results/Long_Fits"
stopifnot(dir.exists(out), dir.exists(long))

EXPECTED <- expand.grid(rep = 1:5, n = c(50L, 250L, 500L, 1000L), scenario = c("1d_skn", "2d_hole"),
                        stringsAsFactors = FALSE)[, c("scenario", "n", "rep")]
rows <- list(); missing <- character(0)
for (k in seq_len(nrow(EXPECTED))) {
  e <- EXPECTED[k, ]; id <- sprintf("%s n=%d r=%d", e$scenario, e$n, e$rep)
  ldir <- file.path(long, e$scenario, sprintf("n%04d", e$n), sprintf("r%02d", e$rep))
  pdir <- file.path(out, e$scenario, sprintf("n%04d", e$n), sprintf("r%02d", e$rep))
  f <- file.path(ldir, "final_status.csv")
  if (!file.exists(f)) { missing <- c(missing, paste0(id, ": long reference fit not finished")); next }
  st <- read.csv(f, colClasses = c(data_id = "character"))
  if (!isTRUE(st$pass)) { missing <- c(missing, paste0(id, ": long reference fit did not pass its convergence rule")); next }
  if (!file.exists(file.path(pdir, "final_status.json"))) { missing <- c(missing, paste0(id, ": replicate-study fit not completed")); next }
  repfit <- readRDS(file.path(pdir, "compact_fit.rds"))
  ref <- readRDS(file.path(dirname(f), "pooled_fit.rds"))
  pz <- as.matrix(repfit$points[, if (st$scenario == "1d_skn") "z" else c("z1", "z2"), drop = FALSE])
  stopifnot(identical(ref$data_id, repfit$summary$data_id),
            max(abs(pz - as.matrix(ref$z))) < 1e-12,
            max(abs(repfit$points$truth - ref$truth)) < 1e-12)
  truth <- repfit$points$truth; tn <- sqrt(mean(truth^2))
  diff <- sqrt(mean((repfit$points$GP - ref$post_mean)^2)) / tn
  e_rep <- rel_L2(repfit$points$GP, truth); e_ref <- rel_L2(ref$post_mean, truth)
  rows[[length(rows) + 1L]] <- data.frame(scenario = st$scenario, n = st$n, rep = st$rep,
    data_id = ref$data_id, reference_effort = sprintf("%d chains x %d sweeps", ref$n_chains, st$n_iter),
    reference_error = e_ref, replicate_study_error = e_rep, absolute_error_change = abs(e_rep - e_ref),
    mean_difference_relative_to_truth = diff,
    replicate_study_MC_noise = repfit$summary$relative_RMS_MC_noise,
    replicate_study_flagged = repfit$summary$flag_for_inspection,
    passes = diff < 0.01 && abs(e_rep - e_ref) < 0.002)
}
res <- if (length(rows)) do.call(rbind, rows) else NULL
if (!is.null(res)) res <- res[order(res$scenario, res$n, res$rep), ]
vdir <- file.path(out, "Validation"); dir.create(vdir, showWarnings = FALSE)
if (!is.null(res)) write.csv(res, file.path(vdir, "cap_validation.csv"), row.names = FALSE)

n_exp <- nrow(EXPECTED); n_cmp <- if (is.null(res)) 0L else nrow(res); n_pass <- if (is.null(res)) 0L else sum(res$passes)
failed <- if (is.null(res)) character(0) else
  with(res[!res$passes, , drop = FALSE], sprintf("%s n=%d r=%d", scenario, n, rep))
verdict <- if (n_pass == n_exp) sprintf("**VALIDATED: all %d of %d comparisons present and passing.**", n_pass, n_exp) else
  sprintf("**NOT VALIDATED: %d of %d comparisons passing (%d missing, %d failed).**",
          n_pass, n_exp, length(missing), length(failed))
lines <- c("# One-chain cap validation against long two-chain fits", "", verdict, "",
  sprintf("Generated %s from %s and %s. No sampling.", format(Sys.time(), "%Y-%m-%d %H:%M"),
          normalizePath(out), normalizePath(long)), "",
  "Tolerances: mean difference < 1% of the truth norm and absolute relative-L2 error change < 0.002.",
  sprintf("Expected: %d comparisons (2 scenarios x 4 sample sizes x datasets 1-5).", n_exp), "")
if (length(missing)) lines <- c(lines, "Missing comparisons:", paste0("- ", missing), "")
if (length(failed)) lines <- c(lines, "Failed comparisons:", paste0("- ", failed), "")
if (n_cmp) {
  lines <- c(lines,
    "| Scenario | n | Dataset | Reference | Ref. error | Replicate-study error | Mean diff. (%) | Error change | Pass |",
    "|---|---:|---:|---|---:|---:|---:|---:|---|")
  for (i in seq_len(n_cmp)) {
    x <- res[i, ]
    lines <- c(lines, sprintf("| %s | %d | %d | %s | %.5f | %.5f | %.3f | %.5f | %s |", x$scenario,
      x$n, x$rep, x$reference_effort, x$reference_error, x$replicate_study_error,
      100 * x$mean_difference_relative_to_truth, x$absolute_error_change, x$passes))
  }
  lines <- c(lines, "", sprintf("Datasets compared: %d; passing: %d; largest mean difference %.3f%%; largest error change %.5f.",
    n_cmp, n_pass, 100 * max(res$mean_difference_relative_to_truth), max(res$absolute_error_change)))
}
writeLines(lines, file.path(vdir, "SUMMARY.md"))
cat(lines, sep = "\n"); cat("\n")
if (n_pass < n_exp)
  stop(sprintf("Cap validation not complete: %d of %d comparisons passing (%d missing, %d failed). See %s.",
               n_pass, n_exp, length(missing), length(failed), file.path(vdir, "SUMMARY.md")), call. = FALSE)
