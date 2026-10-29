# Multiple-chain diagnostics use the maintained posterior package.
# All inputs are post-burn-in draws, with iterations x chains x variables.

pilot_probe_indices <- function(z) {
  d <- ncol(z)
  points <- if (d == 1) matrix(seq(0, 1, length.out = 21), ncol = 1) else
    rbind(as.matrix(expand.grid(seq(0.05, 0.95, length.out = 5),
                               seq(0.05, 0.95, length.out = 5))),
          c(0.3, 0.8), c(0.8, 0.3), c(0.5, 0.5),
          c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  unique(apply(points, 1, function(p) which.min(rowSums(sweep(z, 2, p)^2))))
}

summarise_pilot <- function(paths, out_dir, basic_file = NULL) {
  first <- readRDS(paths[1])
  z <- first$z; truth <- first$truth; d <- ncol(z)
  probe <- pilot_probe_indices(z)
  n_draw <- nrow(first$gp$draws); n_chain <- length(paths)
  variable_names <- c("rho_star", paste0("ell", seq_len(d)), "poisson_loglik",
                      "mean_intensity", sprintf("rho_%04d", probe))
  diagnostics_draws <- array(NA_real_, c(n_draw, n_chain, length(variable_names)),
    dimnames = list(NULL, paste0("chain", seq_len(n_chain)), variable_names))
  # Read one chain at a time; do not keep four complete fitted objects in memory.
  pool <- matrix(NA_real_, n_draw * n_chain, nrow(z))
  means <- matrix(NA_real_, nrow(z), n_chain)
  chain_rows <- vector("list", n_chain)
  reference_id <- first$dataset$data_id
  rm(first); gc(FALSE)
  for (i in seq_along(paths)) {
    f <- readRDS(paths[i]); gp <- f$gp
    stopifnot(identical(f$dataset$data_id, reference_id), identical(f$z, z),
              identical(f$truth, truth), nrow(gp$draws) == n_draw,
              all(is.finite(gp$draws)))
    keep <- gp$keep
    diagnostics_draws[, i, ] <- cbind(gp$trace$rho_star[keep],
      gp$trace$ell[keep, , drop = FALSE],
      gp$trace$loglik[keep] + f$events * log(gp$trace$rho_star[keep]),
      rowMeans(gp$draws), gp$draws[, probe, drop = FALSE])
    pool[((i - 1) * n_draw + 1):(i * n_draw), ] <- gp$draws
    means[, i] <- gp$post_mean
    chain_rows[[i]] <- data.frame(chain = f$chain, data_id = reference_id,
      data_seed = f$dataset$specification$seed, mcmc_seed = f$mcmc_seed,
      relative_L2 = rel_L2(gp$post_mean, truth), seconds = gp$time,
      acceptance_w = gp$acc_w, acceptance_ell = paste(gp$acc_ell, collapse = ";"),
      ell_start = paste(f$initial$ell, collapse = ";"), rho_star_start = f$initial$rho_star)
    rm(f, gp); gc(FALSE)
  }
  diagnostics <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(diagnostics_draws), "mean", "sd", "rhat",
    "ess_bulk", "ess_tail", "mcse_mean"))
  diagnostics$passes_targets <- with(diagnostics,
    is.finite(rhat) & is.finite(ess_bulk) & is.finite(ess_tail) &
    rhat < 1.01 & ess_bulk >= 200 & ess_tail >= 200)   # protocol rule (as run_long_fits.R)
  write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)
  write.csv(do.call(rbind, chain_rows), file.path(out_dir, "chains.csv"), row.names = FALSE)
  probes <- data.frame(variable = sprintf("rho_%04d", probe), evaluation_index = probe, z[probe, , drop = FALSE])
  names(probes)[-(1:2)] <- paste0("z", seq_len(d))
  write.csv(probes, file.path(out_dir, "monitored_points.csv"), row.names = FALSE)
  mean_rho <- colMeans(pool)
  # Column-wise quantiles avoid a full transposed copy of the pooled draw matrix.
  intervals <- vapply(seq_len(ncol(pool)), function(j)
    quantile(pool[, j], c(0.025, 0.975), names = FALSE), numeric(2))
  pooled <- list(data_id = reference_id, z = z, truth = truth,
    post_mean = mean_rho, lower = intervals[1, ], upper = intervals[2, ],
    chain_means = means, draws_per_chain = n_draw, n_chains = n_chain)
  saveRDS(pooled, file.path(out_dir, "pooled_fit.rds"))
  saveRDS(diagnostics_draws, file.path(out_dir, "monitored_draws.rds"))
  rm(pool); gc(FALSE)
  basic <- if (!is.null(basic_file) && file.exists(basic_file)) readRDS(basic_file) else NULL
  if (!is.null(basic)) stopifnot(identical(basic$dataset$data_id, reference_id), identical(basic$z, z))
  summary <- data.frame(scenario = basename(out_dir), data_id = reference_id,
    chains = n_chain, draws_per_chain = n_draw,
    pooled_relative_L2 = rel_L2(mean_rho, truth),
    basic_relative_L2 = if (is.null(basic)) NA_real_ else rel_L2(basic$gp$post_mean, truth),
    max_rhat = max(diagnostics$rhat), min_bulk_ess = min(diagnostics$ess_bulk),
    min_tail_ess = min(diagnostics$ess_tail), all_pass_targets = all(diagnostics$passes_targets),
    interval_truth_fraction = mean(truth >= pooled$lower & truth <= pooled$upper))
  write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)

  cols <- c("#3274A1", "#C44E52", "#55A868", "#8172B3")
  show_traces <- function() {
    selected <- seq_len(d + 2)
    par(mfrow = c(length(selected), 1), mar = c(3, 4, 2, 1))
    for (j in selected) {
      matplot(diagnostics_draws[, , j], type = "l", lty = 1, col = cols,
        xlab = "Retained post-burn-in draw", ylab = variable_names[j], main = variable_names[j])
      if (j == 1) legend("topright", paste("Chain", seq_len(n_chain)), col = cols,
                         lty = 1, bty = "n", cex = 0.8)
    }
  }
  png(file.path(out_dir, "traces.png"), 1600, 500 * (d + 2), res = 160)
  show_traces(); dev.off()
  pdf(file.path(out_dir, "traces.pdf"), 10, 3 * (d + 2))
  show_traces(); dev.off()
  show_fit <- function() {
    if (d == 1) {
      par(mar = c(4, 4, 3, 1))
      plot(z[, 1], truth, type = "n", ylim = range(0, truth, pooled$upper),
           xlab = "Covariate z", ylab = "Intensity", main = "Four-chain posterior fit")
      polygon(c(z[, 1], rev(z[, 1])), c(pooled$lower, rev(pooled$upper)),
              col = adjustcolor(cols[1], 0.15), border = NA)
      matlines(z[, 1], means, col = cols, lty = 2, lwd = 1)
      lines(z[, 1], truth, lwd = 2)
      lines(z[, 1], mean_rho, col = cols[1], lwd = 2)
      legend("topleft", c("Truth", "Pooled mean", "Individual chain means"),
             col = c("black", cols[1], "grey50"), lty = c(1, 1, 2), bty = "n")
    } else {
      values <- c(list(truth, mean_rho), lapply(seq_len(n_chain), function(i) means[, i]))
      titles <- c("Truth", "Pooled posterior mean", paste("Chain", seq_len(n_chain)))
      limits <- range(0, unlist(values)); palette <- hcl.colors(100, "YlOrRd", rev = TRUE)
      xx <- sort(unique(z[, 1])); yy <- sort(unique(z[, 2]))
      par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
      for (i in seq_along(values)) {
        M <- matrix(NA_real_, length(xx), length(yy))
        M[cbind(match(z[, 1], xx), match(z[, 2], yy))] <- values[[i]]
        image(xx, yy, M, zlim = limits, col = palette, xlab = "z1", ylab = "z2",
              main = titles[i], asp = 1, useRaster = TRUE)
      }
    }
  }
  png(file.path(out_dir, "fit_comparison.png"), if (d == 1) 1600 else 1800,
      if (d == 1) 1000 else 1250, res = 160)
  show_fit(); dev.off()
  pdf(file.path(out_dir, "fit_comparison.pdf"), if (d == 1) 10 else 12,
      if (d == 1) 6.25 else 8.33)
  show_fit(); dev.off()
  summary
}
