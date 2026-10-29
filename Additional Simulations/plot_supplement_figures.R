###############################################################################
# plot_supplement_figures.R -- Supplement S4.4-S4.6 figures in the style of the main-text figures
# (../Simulation Study/plot_paper_figures.R: sim1d_fits.pdf and sim2d_maps.pdf). No sampling.
#
# Run from 'Data Analysis/Additional Simulations':  Rscript plot_supplement_figures.R [--out=DIR]
#   (default results folder Results/Additional_Simulations)
# Writes to ../BA 2026_10_xx/Figures_Suppl/ (previews: <out>/Paper_Figures/*.png):
#   sim1d_additional_fits.pdf  Fig:1DResults2   rows: truth (1DTruth2) and (1DTruth3); n = 50, 250, 1000
#   sim2d_slope_maps.pdf       Fig:2DResults3   top: posterior mean n = 50, 250, 1000 and truth;
#   sim2d_twobump_maps.pdf     Fig:2DResults2        bottom: kernel (CV), kernel (oracle),
#   sim2d_detcov_maps.pdf      Fig:2DResultsCoord    log-polynomial (k = 4), n = 1000, same dataset
#   sim2d_overparam_maps.pdf   Fig:2DResultsMiss
#   sim2d_slope_diagonals.pdf  Fig:2Dmarginals        projections on the diagonal and anti-diagonal
#   sim2d_detcov_slices.pdf    Fig:2DResultsCoordMarg slices z2 = 0.8 and z1 = 0.4
#   sim2d_overparam_slices.pdf Fig:2DResultsMissMarg  slices z2 = 0.8 and z1 = 0.3
# Shown datasets: for each scenario and n, the dataset whose errors are closest to the average for
# all shown methods (sum over GP, kernel CV, kernel oracle, log-polynomial of degree 4 of
# |error / mean error - 1|), among all completed datasets, as in plot_paper_figures.R. The 95% bands
# are the pointwise bands of the single replicate-study chain (there are no long fits here).
# Projections at off-grid values (0.8, 0.4, 0.3 lie between the lines of the 60 x 60 reporting
# grid) are bilinear interpolations of the stored estimates; the truth is evaluated exactly.
# Added 19 October 2026; not part of the checksummed runner sources.
###############################################################################
args <- commandArgs(TRUE)
OUT <- if (any(startsWith(args, "--out="))) sub("^--out=", "", args[startsWith(args, "--out=")][1]) else
  file.path("Results", "Additional_Simulations")
source("../R/data_generation.R"); source("additional_scenarios.R")
fig_dir <- file.path(OUT, "Paper_Figures", "pdf"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
prev_dir <- file.path(OUT, "Paper_Figures")
stopifnot(dir.exists(fig_dir)); dir.create(prev_dir, showWarnings = FALSE)
errs <- read.csv(file.path(OUT, "Tables", "method_errors.csv"))
shown <- c("GP", "Kernel_CV", "Kernel_oracle", "Log_poly4")
ns <- c(50L, 250L, 1000L)

representative <- function(s, n) {
  e <- errs[errs$scenario == s & errs$n == n & errs$method %in% shown, ]
  if (!nrow(e)) stop("No results for ", s, " n = ", n)
  m <- tapply(e$relative_L2, e$method, mean)
  score <- tapply(abs(e$relative_L2 / m[e$method] - 1), e$rep, sum)
  as.integer(names(score)[which.min(score)])
}
read_points <- function(s, n) {
  r <- representative(s, n)
  pd <- file.path(OUT, s, sprintf("n%04d", n), sprintf("r%02d", r))
  p <- read.csv(file.path(pd, "points.csv")); sm <- read.csv(file.path(pd, "summary.csv"))
  stopifnot(sm$n == n, sm$rep == r)
  attr(p, "rep") <- r
  p
}
used <- character()
note <- function(s, P) used <<- c(used, sprintf("%s: %s", s, paste(ns, vapply(P, attr, 0L, "rep"), sep = ":", collapse = ", ")))

# Style of plot_paper_figures.R
col_gp <- "#1F6FB4"; col_band <- adjustcolor(col_gp, alpha.f = 0.25)
col_ker <- "#C0392B"; col_pol <- "#2E8B57"
par_panels <- function() par(mar = c(4, 4.2, 2.2, 0.8), mgp = c(2.4, 0.7, 0), las = 1,
                             cex.axis = 1.3, cex.lab = 1.45, cex.main = 1.9)
legend_row <- function() {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", cex = 1.65, x.intersp = 0.6, text.width = NA,
    legend = c("Ground truth", "Posterior mean", "95% credible band", "Kernel (CV)",
               "Kernel (oracle)", "Log-polynomial (k = 4)"),
    col = c("black", col_gp, col_band, col_ker, col_ker, col_pol),
    lwd = c(2, 2, 10, 2, 2, 2), lty = c(1, 1, 1, 1, 2, 1), seg.len = 1.8)
}
curve_panel <- function(t, v, ylim, xlab, ylab, main) {
  plot(NA, xlim = c(0, 1), ylim = ylim, xaxs = "i", xlab = xlab, ylab = ylab, main = main)
  abline(h = pretty(ylim), v = seq(0, 1, 0.25), col = "grey90")
  polygon(c(t, rev(t)), c(v$lower, rev(v$upper)), col = col_band, border = NA)
  lines(t, v$truth, lwd = 2, col = "black")
  lines(t, v$Kernel_CV, lwd = 1.5, col = col_ker)
  lines(t, v$Kernel_oracle, lwd = 1.5, lty = 2, col = col_ker)
  lines(t, v$Log_poly4, lwd = 1.5, col = col_pol)
  lines(t, v$GP, lwd = 2, col = col_gp)
  box()
}
row_label <- function(lab) { par(mar = c(4, 0, 2.2, 0)); plot.new(); text(0.6, 0.5, lab, srt = 90, cex = 1.7) }
save_fig <- function(name, draw, w, h) {
  pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h); draw(); dev.off()
  png(file.path(prev_dir, paste0(name, ".png")), width = 200 * w, height = 200 * h, res = 200); draw(); dev.off()
  cat("Wrote", name, "\n")
}
cols <- c("truth", "GP", "lower", "upper", "Kernel_CV", "Kernel_oracle", "Log_poly4")

## ---------------------------------------------------------------- S4.4: univariate truths
rows1d <- list(list(s = "1d_negexp"), list(s = "1d_minmax"))
P1 <- lapply(rows1d, function(r) { P <- lapply(ns, function(n) read_points(r$s, n)); note(r$s, P); P })
draw1d <- function() {
  K <- length(ns); R <- length(rows1d)
  layout(rbind(matrix(seq_len(K * R), R, byrow = TRUE), rep(K * R + 1, K)), heights = c(rep(1, R), 0.16))
  par_panels()
  for (i in seq_len(R)) {
    P <- P1[[i]]
    ylim <- c(0, 1.03 * max(sapply(P, function(p) max(p[cols]))))
    for (k in seq_len(K))
      curve_panel(P[[k]]$z, P[[k]], ylim, "z", if (k == 1) expression(rho(z)) else "", bquote(n == .(ns[k])))
  }
  legend_row()
}
save_fig("sim1d_additional_fits", draw1d, 12, 8.2)

## ---------------------------------------------------------------- S4.5-S4.6: maps
to_grid <- function(p, v) {
  x <- sort(unique(p$z1)); y <- sort(unique(p$z2))
  M <- matrix(NA_real_, length(x), length(y)); M[cbind(match(p$z1, x), match(p$z2, y))] <- v
  stopifnot(!anyNA(M)); list(x = x, y = y, z = M)
}
pal2d <- hcl.colors(200, "Plasma")
draw_map <- function(g, main, zlim) {
  idx <- 1 + floor((length(pal2d) - 1) * (g$z - zlim[1]) / diff(zlim))
  cm <- matrix(pal2d[pmin(pmax(idx, 1), length(pal2d))], nrow(g$z))   # rows: z1, cols: z2
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = expression(z[1]), ylab = expression(z[2]), main = main)
  rasterImage(as.raster(t(cm)[ncol(cm):1, ]), 0, 0, 1, 1, interpolate = TRUE)
  box()
}
maps <- list()
for (s in c("2d_slope", "2d_twobump", "2d_detcov", "2d_overparam")) {
  P <- lapply(ns, function(n) read_points(s, n)); note(s, P); maps[[s]] <- P
  p1000 <- P[[3]]
  panels <- list(
    list(g = to_grid(P[[1]], P[[1]]$GP), main = "n = 50"),
    list(g = to_grid(P[[2]], P[[2]]$GP), main = "n = 250"),
    list(g = to_grid(p1000, p1000$GP), main = "n = 1000"),
    list(g = to_grid(p1000, p1000$truth), main = "Ground truth"),
    list(g = to_grid(p1000, p1000$Kernel_CV), main = "Kernel (CV)"),
    list(g = to_grid(p1000, p1000$Kernel_oracle), main = "Kernel (oracle)"),
    list(g = to_grid(p1000, p1000$Log_poly4), main = "Log-polynomial (k = 4)"))
  zlim <- c(0, max(sapply(panels, function(a) max(a$g$z))))
  draw2d <- local({ panels <- panels; zlim <- zlim; function() {
    wl <- 0.14
    layout(matrix(c(9, 1, 2, 3, 4, 10, 5, 6, 7, 8), 2, byrow = TRUE), widths = c(wl, 1, 1, 1, 1))
    par(mar = c(4, 4.2, 2.4, 1), mgp = c(2.4, 0.7, 0), las = 1, cex.axis = 1.3, cex.lab = 1.45, cex.main = 1.9)
    for (a in panels) draw_map(a$g, a$main, zlim)
    yy <- seq(zlim[1], zlim[2], length.out = 200)
    par(mar = c(4, 11, 2.4, 9))
    image(1, yy, matrix(yy, 1), col = pal2d, axes = FALSE, xlab = "", ylab = "",
          main = expression(rho(z)), cex.main = 1.9, useRaster = TRUE)
    axis(4, cex.axis = 1.6); box()
    par(mar = c(4, 0, 2.4, 0))
    for (lab in c("Posterior mean", "Other methods, n = 1000")) { plot.new(); text(0.6, 0.5, lab, srt = 90, cex = 1.8, font = 2) }
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
    plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    xs <- (wl + 3) / (wl + 4)
    segments(wl / (wl + 4), 0.5, 1, 0.5, lty = 3, lwd = 1.5, col = "grey35")
    segments(xs, 0.01, xs, 0.99, lty = 3, lwd = 1.5, col = "grey35")
  } })
  save_fig(paste0("sim", s, "_maps"), draw2d, 12.4, 6.6)
}

## ---------------------------------------------------------------- S4.5-S4.6: projections
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
profiles <- list(
  "2d_slope" = list(file = "sim2d_slope_diagonals",
    lines = list(list(label = expression(bold("Diagonal")), xlab = expression(z[1]), q = function(t) cbind(t, t)),
                 list(label = expression(bold("Anti-diagonal")), xlab = expression(z[1]), q = function(t) cbind(t, 1 - t)))),
  "2d_detcov" = list(file = "sim2d_detcov_slices",
    lines = list(list(label = expression(bold("Slice")~bolditalic(z)[2]~bold("= 0.8")), xlab = expression(z[1]), q = function(t) cbind(t, 0.8)),
                 list(label = expression(bold("Slice")~bolditalic(z)[1]~bold("= 0.4")), xlab = expression(z[2]), q = function(t) cbind(0.4, t)))),
  "2d_overparam" = list(file = "sim2d_overparam_slices",
    lines = list(list(label = expression(bold("Slice")~bolditalic(z)[2]~bold("= 0.8")), xlab = expression(z[1]), q = function(t) cbind(t, 0.8)),
                 list(label = expression(bold("Slice")~bolditalic(z)[1]~bold("= 0.3")), xlab = expression(z[2]), q = function(t) cbind(0.3, t)))))
for (s in names(profiles)) {
  P <- maps[[s]]; pr <- profiles[[s]]
  t <- seq(0, 1, length.out = 301)
  V <- lapply(pr$lines, function(L) {
    q <- L$q(t); tr <- ADDSIM_SCENARIOS[[s]]$rho0(q)
    lapply(P, function(p) { v <- lapply(setdiff(cols, "truth"), function(cn) grid_interp(p, cn, q[, 1], q[, 2]))
      names(v) <- setdiff(cols, "truth"); v$truth <- tr; v })
  })
  drawp <- local({ pr <- pr; V <- V; t <- t; function() {
    K <- length(ns); R <- length(pr$lines); wl <- 0.12
    m <- rbind(cbind(K * R + 1 + seq_len(R), matrix(seq_len(K * R), R, byrow = TRUE)), c(0, rep(K * R + 1, K)))
    layout(m, heights = c(rep(1, R), 0.16), widths = c(wl, rep(1, K)))
    par_panels()
    for (i in seq_len(R)) {
      ylim <- c(0, 1.03 * max(sapply(V[[i]], function(v) max(unlist(v)))))
      for (k in seq_len(K))
        curve_panel(t, V[[i]][[k]], ylim, pr$lines[[i]]$xlab, if (k == 1) expression(rho(z)) else "", bquote(n == .(ns[k])))
    }
    legend_row()
    for (i in seq_len(R)) row_label(pr$lines[[i]]$label)
  } })
  save_fig(pr$file, drawp, 12.6, 8.2)
}
writeLines(c("Datasets shown (n:dataset):", used), file.path(prev_dir, "datasets_shown.txt"))
cat("Datasets shown (n:dataset):\n"); writeLines(used)
cat("Figures written to", normalizePath(fig_dir), "\n")
