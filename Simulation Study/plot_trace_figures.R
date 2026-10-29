###############################################################################
# plot_trace_figures.R
# Trace plots for Supplement Section S4.2 (validation of the sampling setup and MCMC
# diagnostics), from the saved long two-chain fits (run_long_fits.R). No sampling.
# Run from 'Data Analysis/Simulation Study':  Rscript plot_trace_figures.R
# Writes <BA folder>/Figures_Suppl/trace_1d.pdf and trace_2d.pdf (dataset 1, n = 1000),
# with PNG previews in Results/Paper_Figures/.
# Shown: rho*, the length-scale(s) and the Poisson log-likelihood, both chains overlaid,
# after burn-in (50000) and thinning (10), as set in run_long_fits.R.
###############################################################################

LONG <- file.path("Results", "Long_Fits")
fig_dir <- file.path("..", "..", "Submitted Manuscripts", "Round 2", "BA 2026_10_xx", "Figures_Suppl")
prev_dir <- file.path("Results", "Paper_Figures")
dir.create(fig_dir, showWarnings = FALSE); dir.create(prev_dir, showWarnings = FALSE)
BURN <- 50000L; THIN <- 10L                      # as in run_long_fits.R
COLS <- c("#B3343F", "#2C6496")                  # chain 1 (red), chain 2 (blue)

# PDF device with Unicode support, for the symbol ell: Quartz on macOS (no XQuartz needed), else cairo;
# without either, plain pdf() with an italic l.
PDF_DEV <- if (isTRUE(capabilities("aqua"))) "quartz" else if (isTRUE(capabilities("cairo"))) "cairo" else "pdf"
ELL <- if (PDF_DEV == "pdf") "l" else "\u2113"
open_pdf <- function(f, width, height) switch(PDF_DEV,
  quartz = quartz(type = "pdf", file = f, width = width, height = height),
  cairo  = cairo_pdf(f, width = width, height = height),
  pdf    = pdf(f, width = width, height = height))

var_label <- function(v, d) switch(v,
  rho_star = expression(rho^"*"),
  poisson_loglik = "Log-likelihood",
  ell1 = if (d == 1) bquote(italic(.(ELL))) else bquote(italic(.(ELL))[1]),
  ell2 = bquote(italic(.(ELL))[2]),
  v)

plot_traces <- function(scenario, file, ncol, width, height, n = 1000L, rep = 1L) {
  dir <- file.path(LONG, scenario, sprintf("n%04d", n), sprintf("r%02d", rep))
  if (!file.exists(file.path(dir, "monitored_draws.rds"))) stop("Long fit not available: ", dir)
  a <- readRDS(file.path(dir, "monitored_draws.rds"))   # draws x chains x variables
  stopifnot(length(dim(a)) == 3, dim(a)[2] == 2)
  v <- dimnames(a)[[3]]
  show <- c("rho_star", grep("^ell[0-9]+$", v, value = TRUE), "poisson_loglik")
  d <- sum(grepl("^ell[0-9]+$", v))
  it <- BURN + THIN * seq_len(dim(a)[1])
  draw <- function() {
    op <- par(mfrow = c(ceiling(length(show) / ncol), ncol), mar = c(3.4, 5.2, 0.8, 0.6),
              mgp = c(3.7, 0.6, 0), oma = c(1.6, 0, 0, 0), cex = 0.85, cex.lab = 1.25, cex.axis = 1.05,
              las = 1)
    on.exit(par(op))
    for (nm in show) {
      y <- a[, , nm]
      plot(range(it) / 1000, range(y), type = "n", xlab = "", ylab = var_label(nm, d))
      title(xlab = expression("Iteration" ~ (""%*%1000)), line = 2.3)
      for (ch in 2:1) lines(it / 1000, y[, ch], col = adjustcolor(COLS[ch], alpha.f = 0.45), lwd = 0.5)
    }
    # common legend in the bottom outer margin, below the panels
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("bottom", legend = c("Chain 1", "Chain 2"), col = COLS, lwd = 2, horiz = TRUE, bty = "n",
           cex = 1, inset = 0.005)
  }
  open_pdf(file.path(fig_dir, file), width = width, height = height); draw(); invisible(dev.off())
  png(file.path(prev_dir, sub("\\.pdf$", ".png", file)), width = width, height = height,
      units = "in", res = 150); draw(); invisible(dev.off())
  cat(sprintf("%s: %s, n = %d, dataset %d, %d draws per chain (iterations %d-%d)\n",
              file, scenario, n, rep, dim(a)[1], min(it), max(it)))
}

plot_traces("1d_skn",  "trace_1d.pdf", ncol = 3, width = 10, height = 2.9)
plot_traces("2d_hole", "trace_2d.pdf", ncol = 2, width = 10, height = 5.1)
