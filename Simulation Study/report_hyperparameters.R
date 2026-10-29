# report_hyperparameters.R
# Paper Table Tab:Hyper: posterior summaries of the hyper-parameters from the long two-chain
# fits (run_long_fits.R), univariate and bivariate settings, n = 50, 250, 500, 1000, datasets 1-5.
# Reads the saved monitored traces only; no sampling.
# Run from 'Data Analysis/Simulation Study':  Rscript report_hyperparameters.R [long_fit_dir]
#   (default Results/Long_Fits)
# Writes Results/Paper_Tables/hyperparameter_results.csv (one row per dataset),
#        Results/Paper_Tables/hyperparameter_table.csv (medians over datasets, by setting and n),
#        Results/Paper_Tables/hyperparameter_table_rows.tex (rows for the LaTeX table).
# Last column of the paper table: minimal bulk ESS and maximal R-hat over the hyper-parameters
# (rho_star and the length-scales). Convergence rule (agreed 9 October 2026): rank-normalised
# R-hat < 1.01 and bulk and tail ESS >= 200 for all monitored quantities.

# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("posterior")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
args <- commandArgs(TRUE)
long <- if (length(args)) args[1] else file.path("Results", "Long_Fits")
out_dir <- file.path("Results", "Paper_Tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
settings <- c("1d_skn" = "univariate", "2d_hole" = "bivariate")
ns <- c(50L, 250L, 500L, 1000L); reps <- 1:5

rows <- list()
for (sc in names(settings)) for (n in ns) for (r in reps) {
  dir <- file.path(long, sc, sprintf("n%04d", n), sprintf("r%02d", r))
  if (!file.exists(file.path(dir, "final_status.csv"))) { message("Long fit not available: ", dir); next }
  st <- read.csv(file.path(dir, "final_status.csv"), colClasses = c(data_id = "character"))
  a <- readRDS(file.path(dir, "monitored_draws.rds"))
  stopifnot(length(dim(a)) == 3, dim(a)[2] == 2)
  repfit <- sprintf("Results/Replicate_Study/%s/n%04d/r%02d/summary.csv", sc, n, r)
  if (file.exists(repfit)) stopifnot(identical(read.csv(repfit, colClasses = c(data_id = "character"))$data_id, st$data_id))
  dg <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a),
    "mean", "rhat", "ess_bulk", "ess_tail"))
  ell <- grepl("^ell[0-9]+$", dg$variable); hyp <- ell | dg$variable == "rho_star"
  g <- function(v) if (v %in% dg$variable) dg$mean[dg$variable == v] else NA_real_
  rows[[length(rows) + 1]] <- data.frame(setting = settings[[sc]], n = n, rep = r, data_id = st$data_id,
    iterations = st$n_iter, ell1 = g("ell1"), ell2 = g("ell2"), rho_star = g("rho_star"),
    ell_min_ess_bulk = min(dg$ess_bulk[ell]), ell_max_rhat = max(dg$rhat[ell]),
    hyper_min_ess_bulk = min(dg$ess_bulk[hyp]), hyper_max_rhat = max(dg$rhat[hyp]),
    all_max_rhat = max(dg$rhat), all_min_ess_bulk = min(dg$ess_bulk), all_min_ess_tail = min(dg$ess_tail),
    passes_rule = all(is.finite(dg$rhat) & dg$rhat < 1.01 & dg$ess_bulk >= 200 & dg$ess_tail >= 200))
}
if (!length(rows)) stop("No long fits found in ", long)
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "hyperparameter_results.csv"), row.names = FALSE)
tab <- do.call(rbind, lapply(split(res, list(res$setting, res$n), drop = TRUE), function(a) data.frame(
  setting = a$setting[1], n = a$n[1], datasets = nrow(a), reps = paste(a$rep, collapse = ","),
  ell1 = median(a$ell1), ell2 = median(a$ell2), rho_star = median(a$rho_star),
  ell_min_ess_bulk = min(a$ell_min_ess_bulk), ell_max_rhat = max(a$ell_max_rhat),
  hyper_min_ess_bulk = min(a$hyper_min_ess_bulk), hyper_max_rhat = max(a$hyper_max_rhat),
  iterations_min = min(a$iterations), iterations_max = max(a$iterations),
  all_pass_rule = all(a$passes_rule))))
tab <- tab[order(tab$setting != "univariate", tab$n), ]
write.csv(tab, file.path(out_dir, "hyperparameter_table.csv"), row.names = FALSE)
tex <- with(tab, sprintf("\t%s & %d & %.1f & %s & %.1f & %.0f / %.5f \\\\",
  ifelse(duplicated(setting), "", ifelse(setting == "univariate", "Univariate", "Bivariate")), n, ell1,
  ifelse(is.na(ell2), "--", sprintf("%.1f", ell2)), rho_star, hyper_min_ess_bulk, hyper_max_rhat))
writeLines(tex, file.path(out_dir, "hyperparameter_table_rows.tex"))
print(tab, row.names = FALSE, digits = 4)
cat("\nLaTeX rows (", file.path(out_dir, "hyperparameter_table_rows.tex"), "):\n", sep = ""); writeLines(tex)
if (any(!res$passes_rule)) {
  cat("\nDatasets failing the rule:\n")
  print(res[!res$passes_rule, c("setting", "n", "rep", "all_max_rhat", "all_min_ess_bulk", "all_min_ess_tail")], row.names = FALSE)
} else cat("\nAll reported datasets satisfy R-hat < 1.01 and bulk/tail ESS >= 200 with two chains.\n")
