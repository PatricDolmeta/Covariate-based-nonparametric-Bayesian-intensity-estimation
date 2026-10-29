# plot_paper_figures.R
# Figures for Section 3.3 of the paper, regenerated from saved results (no sampling).
# Run from 'Data Analysis/Simulation Study':  Rscript plot_paper_figures.R
# Writes <BA folder>/Figures_Main/sim1d_fits.pdf and sim2d_maps.pdf, with previews in
# Results/Paper_Figures/.
#
# Sources (version 2.2):
#   replicate-study fits, Results/Replicate_Study (all methods on the common reporting grid);
#   long two-chain fits, Results/Long_Fits (1D posterior mean and pointwise 95% band).
# Representative datasets: for each n, the dataset whose errors are closest to the average for
# all shown methods (sum over GP, kernel CV, kernel oracle, log-polynomial of degree 4 of
# |error / mean error - 1|); in 1D among datasets 1-5 with a finished long fit, in 2D among all
# completed datasets.

REPSTUDY <- file.path("Results", "Replicate_Study"); LONG <- file.path("Results", "Long_Fits")
fig_dir <- file.path("Results", "Paper_Figures", "pdf"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
prev_dir <- file.path("Results", "Paper_Figures")
stopifnot(dir.exists(fig_dir)); dir.create(prev_dir, showWarnings = FALSE)
shown <- c("GP", "Kernel_CV", "Kernel_oracle", "Log_poly4")
errs <- read.csv(file.path(REPSTUDY, "Tables", "method_errors.csv"))
representative <- function(scenario, n, candidates) {
  e <- errs[errs$scenario == scenario & errs$n == n & errs$method %in% shown & errs$rep %in% candidates, ]
  m <- tapply(e$relative_L2, e$method, mean)
  score <- tapply(abs(e$relative_L2 / m[e$method] - 1), e$rep, sum)
  as.integer(names(score)[which.min(score)])
}

ns <- c(50L, 250L, 1000L)
# 1D: the credible band comes from the long two-chain fit, so the representative dataset is chosen
# among datasets 1-5 that have a finished long fit (all five once the long fits are complete).
has_long <- function(n) Filter(function(r) file.exists(file.path(LONG, "1d_skn", sprintf("n%04d", n),
                                 sprintf("r%02d", r), "pooled_fit.rds")), 1:5)
rep1d <- vapply(ns, function(n) {
  cand <- has_long(n); if (!length(cand)) stop("No finished 1D long fit for n = ", n)
  representative("1d_skn", n, cand)
}, integer(1))
read_panel <- function(n, r) {
  pd <- file.path(REPSTUDY, "1d_skn", sprintf("n%04d", n), sprintf("r%02d", r))
  lf <- file.path(LONG, "1d_skn", sprintf("n%04d", n), sprintf("r%02d", r))
  p <- read.csv(file.path(pd, "points.csv")); s <- read.csv(file.path(pd, "summary.csv"))
  pooled <- readRDS(file.path(lf, "pooled_fit.rds"))
  stopifnot(s$n == n, s$rep == r, identical(as.character(s$data_id), pooled$data_id),
            max(abs(p$z - pooled$z[, 1])) < 1e-12, max(abs(p$truth - pooled$truth)) < 1e-10)
  p$GP <- pooled$post_mean; p$lower <- pooled$lower; p$upper <- pooled$upper   # long-fit summaries
  p$quartic <- p$Log_poly4
  p
}
panels <- Map(read_panel, ns, rep1d)

col_gp <- "#1F6FB4"; col_band <- adjustcolor(col_gp, alpha.f = 0.25)
col_ker <- "#C0392B"; col_pol <- "#2E8B57"
ymax <- max(sapply(panels, function(p) max(p$upper, p$GP, p$Kernel_CV, p$Kernel_oracle, p$quartic, p$truth)))
ylim <- c(0, 1.03 * ymax)

draw <- function() {
  K <- length(ns)
  layout(matrix(c(seq_len(K), rep(K + 1, K)), 2, byrow = TRUE), heights = c(1, 0.16))
  par(mar = c(4, 4.2, 2.2, 0.8), mgp = c(2.4, 0.7, 0), las = 1,
      cex.axis = 1.3, cex.lab = 1.45, cex.main = 1.9)
  for (k in seq_along(ns)) {
    p <- panels[[k]]
    plot(NA, xlim = c(0, 1), ylim = ylim, xaxs = "i", xlab = "z",
         ylab = if (k == 1) expression(rho(z)) else "", main = bquote(n == .(ns[k])))
    abline(h = pretty(ylim), v = seq(0, 1, 0.25), col = "grey90")
    polygon(c(p$z, rev(p$z)), c(p$lower, rev(p$upper)), col = col_band, border = NA)
    lines(p$z, p$truth, lwd = 2, col = "black")
    lines(p$z, p$Kernel_CV, lwd = 1.5, col = col_ker)
    lines(p$z, p$Kernel_oracle, lwd = 1.5, lty = 2, col = col_ker)
    lines(p$z, p$quartic, lwd = 1.5, col = col_pol)
    lines(p$z, p$GP, lwd = 2, col = col_gp)
    box()
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", cex = 1.65, x.intersp = 0.6, text.width = NA,
    legend = c("Ground truth", "Posterior mean", "95% credible band", "Kernel (CV)",
               "Kernel (oracle)", "Log-polynomial (k = 4)"),
    col = c("black", col_gp, col_band, col_ker, col_ker, col_pol),
    lwd = c(2, 2, 10, 2, 2, 2), lty = c(1, 1, 1, 1, 2, 1), seg.len = 1.8)
}
pdf(file.path(fig_dir, "sim1d_fits.pdf"), width = 12, height = 4.4); draw(); dev.off()
png(file.path(prev_dir, "sim1d_fits.png"), width = 2400, height = 880, res = 200); draw(); dev.off()
cat("Wrote sim1d_fits.pdf (datasets", paste(ns, rep1d, sep = ":", collapse = ", "), ")\n")

## ---------------------------------------------------------------- Figure sim2d_maps
# Bivariate scenario, posterior means from the replicate-study fits (60 x 60 reporting grid).
# Top row: posterior mean at n = 50, 250, 1000 and the ground truth.
# Bottom row (n = 1000, same dataset): kernel (CV), kernel (oracle), log-polynomial of degree 4.
reps2d <- setNames(vapply(ns, function(n) representative("2d_hole", n, 1:50), integer(1)), ns)
read_2d <- function(n, r) {
  pd <- file.path(REPSTUDY, "2d_hole", sprintf("n%04d", n), sprintf("r%02d", r))
  p <- read.csv(file.path(pd, "points.csv")); s <- read.csv(file.path(pd, "summary.csv"))
  stopifnot(s$n == n, s$rep == r)
  p$quartic <- p$Log_poly4
  p
}
maps <- Map(read_2d, as.integer(names(reps2d)), reps2d)
to_grid <- function(p, v) {
  x <- sort(unique(p$z1)); y <- sort(unique(p$z2))
  M <- matrix(NA_real_, length(x), length(y)); M[cbind(match(p$z1, x), match(p$z2, y))] <- v
  stopifnot(!anyNA(M)); list(x = x, y = y, z = M)
}
p1000 <- maps[[3]]
panels2d <- list(
  list(g = to_grid(maps[[1]], maps[[1]]$GP), main = "n = 50"),
  list(g = to_grid(maps[[2]], maps[[2]]$GP), main = "n = 250"),
  list(g = to_grid(p1000, p1000$GP), main = "n = 1000"),
  list(g = to_grid(p1000, p1000$truth), main = "Ground truth"),
  list(g = to_grid(p1000, p1000$Kernel_CV), main = "Kernel (CV)"),
  list(g = to_grid(p1000, p1000$Kernel_oracle), main = "Kernel (oracle)"),
  list(g = to_grid(p1000, p1000$quartic), main = "Log-polynomial (k = 4)"))
zlim2d <- c(0, max(sapply(panels2d, function(a) max(a$g$z))))
pal2d <- hcl.colors(200, "Plasma")

smooth2d <- TRUE   # display: bilinear interpolation of the 60 x 60 grid (FALSE = raw pixels)
draw_map <- function(g, main) {
  if (smooth2d) {
    idx <- 1 + floor((length(pal2d) - 1) * (g$z - zlim2d[1]) / diff(zlim2d))
    cols <- matrix(pal2d[pmin(pmax(idx, 1), length(pal2d))], nrow(g$z))   # rows: z1, cols: z2
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
         xlab = expression(z[1]), ylab = expression(z[2]), main = main)
    rasterImage(as.raster(t(cols)[ncol(cols):1, ]), 0, 0, 1, 1, interpolate = TRUE)
  } else {
    image(g$x, g$y, g$z, zlim = zlim2d, col = pal2d, xlim = c(0, 1), ylim = c(0, 1),
          xaxs = "i", yaxs = "i", useRaster = TRUE,
          xlab = expression(z[1]), ylab = expression(z[2]), main = main)
  }
  box()
}
draw2d <- function() {
  wl <- 0.14                                            # width of the row-label column
  layout(matrix(c(9, 1, 2, 3, 4, 10, 5, 6, 7, 8), 2, byrow = TRUE), widths = c(wl, 1, 1, 1, 1))
  par(mar = c(4, 4.2, 2.4, 1), mgp = c(2.4, 0.7, 0), las = 1,
      cex.axis = 1.3, cex.lab = 1.45, cex.main = 1.9)
  for (a in panels2d) draw_map(a$g, a$main)
  yy <- seq(zlim2d[1], zlim2d[2], length.out = 200)    # colour bar
  par(mar = c(4, 11, 2.4, 9))
  image(1, yy, matrix(yy, 1), col = pal2d, axes = FALSE, xlab = "", ylab = "",
        main = expression(rho(z)), cex.main = 1.9, useRaster = TRUE)
  axis(4, cex.axis = 1.6); box()
  par(mar = c(4, 0, 2.4, 0))                            # row labels
  for (lab in c("Posterior mean", "Other methods, n = 1000")) {
    plot.new(); text(0.6, 0.5, lab, srt = 90, cex = 1.8, font = 2)
  }
  # dotted separators: between the rows, and between the estimates and the last column
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  xs <- (wl + 3) / (wl + 4)
  segments(wl / (wl + 4), 0.5, 1, 0.5, lty = 3, lwd = 1.5, col = "grey35")
  segments(xs, 0.01, xs, 0.99, lty = 3, lwd = 1.5, col = "grey35")
}
pdf(file.path(fig_dir, "sim2d_maps.pdf"), width = 12.4, height = 6.6); draw2d(); dev.off()
png(file.path(prev_dir, "sim2d_maps.png"), width = 2480, height = 1320, res = 200); draw2d(); dev.off()
cat("Wrote sim2d_maps.pdf (datasets", paste(names(reps2d), reps2d, sep = ":", collapse = ", "), ")\n")
