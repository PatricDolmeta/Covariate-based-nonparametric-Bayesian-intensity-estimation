# make_error_table_tex.R
# LaTeX body of the paper's error table (Tab:Errors) from the replicate-study results.
# No sampling; reads Results/Replicate_Study/Tables/error_table.csv (written by run_replicate_study.R)
# and writes Results/Paper_Tables/error_table_rows.tex.
# Cells: mean relative L2 error over the datasets (between-dataset SD); bold = smallest mean in each
# column within each setting. Added 16 October 2026; not part of the checksummed runner sources.
#
# Run from 'Data Analysis/Simulation Study':  Rscript make_error_table_tex.R

e <- read.csv(file.path("Results", "Replicate_Study", "Tables", "error_table.csv"), stringsAsFactors = FALSE)
ns <- c(50L, 250L, 500L, 1000L)
k <- unique(e$datasets)
if (length(k) != 1) warning("Cells have different numbers of datasets: ", paste(sort(unique(e$datasets)), collapse = ", "))
rows <- list(c("GP", "Posterior mean"), NULL, c("Kernel_CV", "Kernel, CV"), c("Kernel_oracle", "Kernel, oracle"), NULL,
             c("Log_poly1", "Log-polynomial, $k=1$"), c("Log_poly2", "Log-polynomial, $k=2$"),
             c("Log_poly3", "Log-polynomial, $k=3$"), c("Log_poly4", "Log-polynomial, $k=4$"))
sep <- "\t\\noalign{{\\color{gray!60}\\hrule height 0.3pt}}"
out <- c("\\begin{tabular}{l|cccc}", "\t& $n=50$ & $n=250$ & $n=500$ & $n=1000$ \\\\", "\t\\hline")
for (s in list(c("1d_skn", "Univariate setting"), c("2d_hole", "Bivariate setting"))) {
  x <- e[e$scenario == s[1], ]
  best <- vapply(ns, function(n) min(x$mean_relative_L2[x$n == n]), 0)
  out <- c(out, "\t\\noalign{\\vskip 3pt}", sprintf("\t\\multicolumn{5}{c}{\\textbf{%s}}\\\\", s[2]), "\t\\noalign{\\vskip 3pt}")
  for (r in rows) {
    if (is.null(r)) { out <- c(out, sep); next }
    cells <- vapply(seq_along(ns), function(j) {
      y <- x[x$n == ns[j] & x$method == r[1], ]
      stopifnot(nrow(y) == 1)
      m <- sprintf("%.3f", y$mean_relative_L2)
      if (abs(y$mean_relative_L2 - best[j]) < 1e-12) m <- sprintf("\\textbf{%s}", m)
      sprintf("%s (%.3f)", m, y$sd_relative_L2)
    }, "")
    out <- c(out, sprintf("\t%s & %s \\\\", r[2], paste(cells, collapse = " & ")))
  }
  out <- c(out, "\t\\hline")
}
out <- c(out, "\\end{tabular}",
  sprintf("%% Caption: Mean relative $L^2$-errors (and between-dataset standard deviations) over %s independent datasets for each setting and sample size.",
          paste(sort(unique(e$datasets)), collapse = "/")))
dir.create(file.path("Results", "Paper_Tables"), recursive = TRUE, showWarnings = FALSE)
f <- file.path("Results", "Paper_Tables", "error_table_rows.tex")
writeLines(out, f)
cat(out, sep = "\n"); cat("\nWritten:", f, "\n")
