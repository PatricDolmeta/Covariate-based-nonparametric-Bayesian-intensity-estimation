###############################################################################
# plot_additional_figures.R -- preview figures for Supplement S4.4-S4.6 (no sampling)
#
# Run from 'Data Analysis/Additional Simulations' after run_additional_simulations.R:
#   Rscript plot_additional_figures.R                      # reads Results/Additional_Simulations
#   Rscript plot_additional_figures.R --out=DIR            # another results folder
# Writes PDFs to <out>/Figures/. Nothing is written to the manuscript folder.
#
# Figures (as in the Supplement):
#   1d_<truth>.pdf        S4.4: posterior mean, pointwise 95% band, kernel CV and truth, n = 250, 500, 1000
#   2d_<truth>_maps.pdf   S4.5-S4.6: posterior means for n = 50, 250, 1000 and the truth
#   2d_slope_diagonals.pdf, 2d_detcov_slices.pdf, 2d_overparam_slices.pdf: one-dimensional
#                         projections (diagonals; slices z2 = 0.8 and z1 = 0.4 or 0.3)
# For each scenario and n the dataset shown is the one whose errors are closest to the average
# (sum over GP, kernel CV, kernel oracle, log-polynomial of degree 4 of |error / mean error - 1|),
# as in ../Simulation Study/plot_paper_figures.R. Off-grid slices (z = 0.8, 0.4, 0.3 lie between
# lines of the 60 x 60 reporting grid) are bilinear interpolations of the stored posterior mean and
# pointwise quantiles; the truth is evaluated exactly.
###############################################################################
args <- commandArgs(TRUE)
OUT <- if (any(startsWith(args, "--out="))) sub("^--out=", "", args[startsWith(args, "--out=")][1]) else
  file.path("Results", "Additional_Simulations")
source("../R/data_generation.R"); source("additional_scenarios.R")
ef <- file.path(OUT, "Tables", "method_errors.csv")
if (!file.exists(ef)) stop("No results yet: ", ef, " (run run_additional_simulations.R first)")
errs <- read.csv(ef)
fig_dir <- file.path(OUT, "Figures"); dir.create(fig_dir, showWarnings = FALSE)
shown <- c("GP", "Kernel_CV", "Kernel_oracle", "Log_poly4")
col_gp <- "#3274A1"; col_band <- adjustcolor("#3274A1", 0.25); col_k <- "#C44E52"

representative <- function(s, n) {
  e <- errs[errs$scenario == s & errs$n == n & errs$method %in% shown, ]
  if (!nrow(e)) return(NA_integer_)
  m <- tapply(e$relative_L2, e$method, mean)
  score <- tapply(abs(e$relative_L2 / m[e$method] - 1), e$rep, sum)
  as.integer(names(score)[which.min(score)])
}
read_points <- function(s, n) {
  r <- representative(s, n)
  if (is.na(r)) return(NULL)
  p <- read.csv(file.path(OUT, s, sprintf("n%04d", n), sprintf("r%02d", r), "points.csv"))
  attr(p, "rep") <- r
  p
}
available <- function(s, ns) Filter(function(n) any(errs$scenario == s & errs$n == n), ns)

# Bilinear interpolation of a column of a 2D points table (regular reporting grid) at (q1, q2).
grid_interp <- function(p, col, q1, q2) {
  xx <- sort(unique(p$z1)); k <- length(xx)
  M <- matrix(NA_real_, k, k); M[cbind(match(p$z1, xx), match(p$z2, xx))] <- p[[col]]
  loc <- function(q) {
    t <- (pmin(pmax(q, xx[1]), xx[k]) - xx[1]) / (xx[2] - xx[1])
    i <- pmin(floor(t + 1e-9), k - 2); list(i = i + 1, f = pmax(t - i, 0))
  }
  a <- loc(q1); b <- loc(q2)
  (1 - a$f) * (1 - b$f) * M[cbind(a$i, b$i)] + a$f * (1 - b$f) * M[cbind(a$i + 1, b$i)] +
    (1 - a$f) * b$f * M[cbind(a$i, b$i + 1)] + a$f * b$f * M[cbind(a$i + 1, b$i + 1)]
}

# ---------------------------------------------------------------- S4.4: univariate truths
for (s in c("1d_negexp", "1d_minmax")) {
  ns <- available(s, c(250L, 500L, 1000L)); if (!length(ns)) next
  P <- lapply(ns, function(n) read_points(s, n))
  ylim <- range(0, unlist(lapply(P, function(p) c(p$truth, p$lower, p$upper, p$Kernel_CV))))
  pdf(file.path(fig_dir, paste0(s, ".pdf")), width = 4 * length(ns), height = 3.4)
  par(mfrow = c(1, length(ns)), mar = c(4, 4, 2.5, 1))
  for (k in seq_along(ns)) {
    p <- P[[k]]
    plot(p$z, p$truth, type = "n", ylim = ylim, xlab = "z", ylab = "Intensity",
         main = sprintf("n = %d (dataset %d)", ns[k], attr(p, "rep")))
    polygon(c(p$z, rev(p$z)), c(p$lower, rev(p$upper)), col = col_band, border = NA)
    lines(p$z, p$GP, col = col_gp, lwd = 2); lines(p$z, p$Kernel_CV, col = col_k, lwd = 1.5)
    lines(p$z, p$truth, lwd = 1.5)
    if (k == 1) legend("topright", c("Truth", "Posterior mean", "95% band", "Kernel CV"),
                       col = c("black", col_gp, col_band, col_k), lwd = c(1.5, 2, 8, 1.5), bty = "n", cex = 0.8)
  }
  dev.off()
}

# ---------------------------------------------------------------- S4.5-S4.6: maps
ns2 <- c(50L, 250L, 1000L)
for (s in c("2d_slope", "2d_twobump", "2d_detcov", "2d_overparam")) {
  ns <- available(s, ns2); if (!length(ns)) next
  P <- lapply(ns, function(n) read_points(s, n))
  xx <- sort(unique(P[[1]]$z1))
  to_mat <- function(p, col) { M <- matrix(NA_real_, length(xx), length(xx))
    M[cbind(match(p$z1, xx), match(p$z2, xx))] <- p[[col]]; M }
  lim <- range(0, P[[1]]$truth, unlist(lapply(P, `[[`, "GP")))
  pal <- hcl.colors(100, "YlOrRd", rev = TRUE)
  pdf(file.path(fig_dir, paste0(s, "_maps.pdf")), width = 3.6 * (length(ns) + 1) + 0.8, height = 3.6)
  layout(matrix(seq_len(length(ns) + 2), 1), widths = c(rep(1, length(ns) + 1), 0.22))
  par(mar = c(4, 4, 2.5, 0.8))
  for (k in seq_along(ns))
    image(xx, xx, to_mat(P[[k]], "GP"), zlim = lim, col = pal, asp = 1, useRaster = TRUE,
          xlab = expression(z[1]), ylab = expression(z[2]),
          main = sprintf("n = %d (dataset %d)", ns[k], attr(P[[k]], "rep")))
  image(xx, xx, to_mat(P[[1]], "truth"), zlim = lim, col = pal, asp = 1, useRaster = TRUE,
        xlab = expression(z[1]), ylab = expression(z[2]), main = "Truth")
  par(mar = c(4, 0.5, 2.5, 3))
  image(c(0, 1), seq(lim[1], lim[2], length.out = 101), matrix(seq(lim[1], lim[2], length.out = 100), 1),
        col = pal, axes = FALSE, xlab = "", ylab = ""); axis(4, las = 1)
  dev.off()
}

# ---------------------------------------------------------------- S4.5-S4.6: projections
profiles <- list(
  "2d_slope" = list(file = "2d_slope_diagonals.pdf",
    lines = list(list(label = "Diagonal z2 = z1", q = function(t) cbind(t, t)),
                 list(label = "Anti-diagonal z2 = 1 - z1", q = function(t) cbind(t, 1 - t)))),
  "2d_detcov" = list(file = "2d_detcov_slices.pdf",
    lines = list(list(label = "Slice z2 = 0.8", q = function(t) cbind(t, 0.8), axis = "z1"),
                 list(label = "Slice z1 = 0.4", q = function(t) cbind(0.4, t), axis = "z2"))),
  "2d_overparam" = list(file = "2d_overparam_slices.pdf",
    lines = list(list(label = "Slice z2 = 0.8", q = function(t) cbind(t, 0.8), axis = "z1"),
                 list(label = "Slice z1 = 0.3", q = function(t) cbind(0.3, t), axis = "z2"))))
for (s in names(profiles)) {
  ns <- available(s, ns2); if (!length(ns)) next
  P <- lapply(ns, function(n) read_points(s, n)); pr <- profiles[[s]]
  t <- sort(unique(P[[1]]$z1))
  pdf(file.path(fig_dir, pr$file), width = 4 * length(ns), height = 3.2 * length(pr$lines))
  par(mfrow = c(length(pr$lines), length(ns)), mar = c(4, 4, 2.5, 1))
  for (L in pr$lines) {
    q <- L$q(t); truth <- ADDSIM_SCENARIOS[[s]]$rho0(q)
    V <- lapply(P, function(p) sapply(c("GP", "lower", "upper"), function(cn) grid_interp(p, cn, q[, 1], q[, 2])))
    ylim <- range(0, truth, unlist(V))
    for (k in seq_along(ns)) {
      v <- V[[k]]
      plot(t, truth, type = "n", ylim = ylim, xlab = if (is.null(L$axis)) expression(z[1]) else L$axis,
           ylab = "Intensity", main = sprintf("%s, n = %d", L$label, ns[k]))
      polygon(c(t, rev(t)), c(v[, "lower"], rev(v[, "upper"])), col = col_band, border = NA)
      lines(t, v[, "GP"], col = col_gp, lwd = 2); lines(t, truth, lwd = 1.5)
    }
  }
  dev.off()
}
cat("Figures written to", normalizePath(fig_dir), "\n")
