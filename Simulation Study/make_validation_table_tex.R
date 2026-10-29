# make_validation_table_tex.R
# LaTeX body of the Supplement's validation table (Tab:Validation) from the per-dataset comparisons
# written by validate_replicate_study_cap.R (Results/Replicate_Study/Validation/cap_validation.csv).
# No sampling. Writes Results/Paper_Tables/validation_table_rows.tex.
# Per scenario and n: lengths of the long two-chain runs (range over the compared datasets), the maximal
# difference of the posterior means (as computed by the validation script: RMS difference on the
# reporting grid divided by the RMS of the true intensity), and the maximal absolute change in the
# relative L2 error. Comparisons that exceed the protocol tolerances are listed in a LaTeX comment.
# Added 19 October 2026; not part of the checksummed runner sources.
#
# Run from 'Data Analysis/Simulation Study':  Rscript make_validation_table_tex.R

v <- read.csv(file.path("Results", "Replicate_Study", "Validation", "cap_validation.csv"), stringsAsFactors = FALSE)
v$sweeps <- as.integer(sub("^.* x ([0-9]+) sweeps$", "\\1", v$reference_effort))
big <- function(x) gsub(",", "{,}", formatC(x, format = "d", big.mark = ","))
out <- c("\\begin{tabular}{ll|ccc}",
         "\t& & Long runs & Max.~difference of & Max.~change \\\\",
         "\t& $n$ & (iterations) & posterior means & in estimation error \\\\",
         "\t\\hline")
for (s in list(c("1d_skn", "Univariate"), c("2d_hole", "Bivariate"))) {
  ns <- sort(unique(v$n[v$scenario == s[1]]))
  for (i in seq_along(ns)) {
    a <- v[v$scenario == s[1] & v$n == ns[i], ]
    r <- range(a$sweeps)
    len <- if (r[1] == r[2]) sprintf("$2\\times%s$", big(r[1])) else sprintf("$2\\times(%s$--$%s)$", big(r[1]), big(r[2]))
    lab <- c(sprintf("\\textbf{%s}", s[2]), "\\textbf{setting}", "", "")[i]
    out <- c(out, sprintf("\t%s & %d & %s & $%.2f\\%%$ & $%.4f$ \\\\  %% %d datasets: %s", lab, ns[i], len,
                          100 * max(a$mean_difference_relative_to_truth), max(a$absolute_error_change),
                          nrow(a), paste(sort(a$rep), collapse = ",")))
  }
  out <- c(out, "\t\\hline")
}
out <- c(out, "\\end{tabular}")
fail <- v[!v$passes, ]
out <- c(out, sprintf("%% %d comparisons, %d passing the protocol tolerances (mean difference < 1%%, error change < 0.002).",
                      nrow(v), sum(v$passes)))
if (nrow(fail)) out <- c(out, sprintf("%% Exceeding a tolerance: %s n=%d r=%d (difference %.2f%%, error change %.4f)",
                                       fail$scenario, fail$n, fail$rep, 100 * fail$mean_difference_relative_to_truth,
                                       fail$absolute_error_change))
dir.create(file.path("Results", "Paper_Tables"), recursive = TRUE, showWarnings = FALSE)
f <- file.path("Results", "Paper_Tables", "validation_table_rows.tex")
writeLines(out, f)
cat(out, sep = "\n"); cat("\nWritten:", f, "\n")
