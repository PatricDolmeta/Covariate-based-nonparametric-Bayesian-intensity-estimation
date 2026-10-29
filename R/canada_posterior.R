###############################################################################
# canada_posterior.R
# Convergence diagnostics for a fit (several independent chains) and posterior summaries of
# rho computed from the stored node draws (w, rho*): for any set of points z,
#   rho^(s)(z) = rho*^(s) * sigmoid( A(z) w^(s) ),  A(z) = the sampler's interpolation matrix,
# so posterior means / quantiles at arbitrary z (curves, slices, every pixel of a year) are
# computed exactly for the discretised model, without any smoothing.
# Requires gp_sampler.R (interp_matrix) and the 'posterior' package.
###############################################################################

chain_file <- function(fit_dir, chain) file.path(fit_dir, sprintf("chain%d.rds", chain))
chain_target_file <- function(fit_dir, chain) file.path(fit_dir, sprintf("chain%d_target.txt", chain))
chain_target <- function(fit_dir, chain) {
  f <- chain_target_file(fit_dir, chain)
  if (!file.exists(f) || !file.exists(chain_file(fit_dir, chain))) return(0L)
  as.integer(readLines(f, n = 1))
}

# ------------------------------------------------------------------ pooled node draws
load_pool <- function(fit_dir, chains = NULL) {
  if (is.null(chains)) chains <- as.integer(sub("^chain([0-9]+)\\.rds$", "\\1",
    list.files(fit_dir, pattern = "^chain[0-9]+\\.rds$")))
  if (!length(chains)) stop("No chain results in ", fit_dir)
  W <- list(); rs <- list(); first <- NULL
  for (k in sort(chains)) {
    f <- readRDS(chain_file(fit_dir, k))
    if (is.null(first)) first <- f[c("fit_id", "region", "model", "covs", "d", "nd", "target")] else
      stopifnot(identical(f$nd, first$nd), identical(f$covs, first$covs), f$target == first$target)
    W[[length(W) + 1]] <- f$gp$node_w; rs[[length(rs) + 1]] <- f$gp$node_rs
    rm(f); gc(FALSE)
  }
  W <- do.call(rbind, W); rs <- unlist(rs)
  if (!nrow(W)) stop("No stored node draws in ", fit_dir)
  c(first, list(Wt = t(W), rs = rs, chains = sort(chains)))     # nodes x draws
}

posterior_rho <- function(pool, z, probs = NULL, max_cells = 2e7) {
  Wt <- if (!is.null(pool$Wt)) pool$Wt else t(pool$W)        # nodes x draws
  z <- as.matrix(z); nz <- nrow(z); S <- ncol(Wt)
  m <- numeric(nz); Q <- if (is.null(probs)) NULL else matrix(NA_real_, nz, length(probs))
  chunk <- max(1L, floor(max_cells / S))
  for (st in seq(1L, nz, by = chunk)) {
    idx <- st:min(nz, st + chunk - 1L)
    A <- interp_matrix(pool$nd, z[idx, , drop = FALSE])
    R <- plogis(as.matrix(A %*% Wt))                       # |idx| x S
    R <- R * rep(pool$rs, each = length(idx))
    m[idx] <- rowMeans(R)
    if (!is.null(probs))
      Q[idx, ] <- matrix(t(apply(R, 1, quantile, probs = probs, names = FALSE)), length(idx))
  }
  list(mean = m, quantiles = Q)
}

# Plug-in posterior mean of the spatial intensity, rho_hat(Z^(i)(x)), at the valid pixels
plugin_intensity <- function(pool, R, year) {
  i <- match(as.character(year), as.character(R$years))
  if (is.na(i)) stop("Year not in data: ", year)
  z <- do.call(cbind, lapply(pool$covs, function(cv) R$U[[cv]][, i]))
  posterior_rho(pool, z)$mean
}

# ------------------------------------------------------------------ diagnostics
# Monitored quantities: rho*, ell_j, full Poisson log-likelihood, rho at the probe points.
# Criteria (defaults): rank-normalised R-hat < 1.01, combined bulk and tail ESS >= 200 for all
# monitored quantities, and relative RMS Monte Carlo noise of the posterior mean at the probes
# sqrt(mean(MCSE^2)) / sqrt(mean(mean^2)) < 1%.
assess_fit <- function(fit_dir, chains, crit = list(rhat = 1.01, ess = 200, mc_noise = 0.01)) {
  fits <- lapply(chains, function(k) {
    f <- readRDS(chain_file(fit_dir, k)); f$gp$node_w <- NULL; f
  })
  f1 <- fits[[1]]; d <- f1$d
  stopifnot(all(vapply(fits, function(f) f$target == f1$target, logical(1))),
            all(vapply(fits, function(f) nrow(f$gp$draws) == nrow(f1$gp$draws), logical(1))))
  probe_names <- sprintf("rho_probe%02d", seq_len(nrow(f1$z_eval)))
  vars <- c("rho_star", paste0("ell", seq_len(d)), "loglik", probe_names)
  arr <- array(NA_real_, c(nrow(f1$gp$draws), length(fits), length(vars)),
               dimnames = list(NULL, paste0("chain", chains), vars))
  for (k in seq_along(fits)) {
    gp <- fits[[k]]$gp; keep <- gp$keep
    arr[, k, ] <- cbind(gp$trace$rho_star[keep], gp$trace$ell[keep, , drop = FALSE],
                        gp$trace$loglik[keep] + fits[[k]]$N * log(gp$trace$rho_star[keep]),
                        gp$draws)
  }
  diag <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(arr),
    "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  pooled <- apply(arr, 3, function(x) quantile(as.vector(x), c(0.025, 0.975), names = FALSE))
  diag$q2.5 <- pooled[1, ]; diag$q97.5 <- pooled[2, ]
  pr <- diag$variable %in% probe_names
  mc_noise <- sqrt(mean(diag$mcse_mean[pr]^2)) / sqrt(mean(diag$mean[pr]^2))
  diag$passes <- is.finite(diag$rhat) & diag$rhat < crit$rhat &
    is.finite(diag$ess_bulk) & diag$ess_bulk >= crit$ess &
    is.finite(diag$ess_tail) & diag$ess_tail >= crit$ess
  passed <- all(diag$passes) && is.finite(mc_noise) && mc_noise < crit$mc_noise
  write.csv(diag, file.path(fit_dir, "diagnostics.csv"), row.names = FALSE)

  chains_df <- do.call(rbind, lapply(fits, function(f) data.frame(
    chain = f$chain, seed = f$seed, sweeps = f$target, seconds = f$gp$time,
    acc_w = f$gp$acc_w, acc_ell = paste(round(f$gp$acc_ell, 3), collapse = ";"),
    beta = f$gp$beta, s_ell = paste(signif(f$gp$s_ell, 3), collapse = ";"))))
  write.csv(chains_df, file.path(fit_dir, "chains.csv"), row.names = FALSE)

  cols <- rep_len(c("#3274A1", "#C44E52", "#55A868", "#8172B3"), length(fits))
  traces <- function() {
    sel <- seq_len(d + 2)
    par(mfrow = c(length(sel), 1), mar = c(3, 4.5, 2, 1))
    for (j in sel) {
      matplot(arr[, , j], type = "l", lty = 1, col = cols,
              xlab = "Retained post-burn-in draw", ylab = vars[j],
              main = sprintf("%s  (R-hat %.3f, bulk ESS %.0f)", vars[j], diag$rhat[j], diag$ess_bulk[j]))
    }
  }
  pdf(file.path(fit_dir, "traces.pdf"), 9, 2.4 * (d + 2)); tryCatch(traces(), finally = dev.off())

  fmt_ci <- function(v) with(diag[diag$variable == v, ],
                             sprintf("%.4g [%.4g, %.4g]", mean, q2.5, q97.5))
  data.frame(fit_id = f1$fit_id, region = f1$region, model = f1$model, d = d,
             nodes = nrow(f1$nd$nodes), G = f1$G, bins_dom = f1$bins[["dom"]], bins_obs = f1$bins[["obs"]],
             n = f1$n, events = f1$N, chains = length(fits), sweeps = f1$target,
             burn = f1$burn, thin = f1$thin,
             median_minutes_per_chain = median(chains_df$seconds) / 60,
             seconds_per_1000_sweeps = median(chains_df$seconds) / f1$target * 1000,
             acc_w = mean(chains_df$acc_w),
             acc_ell = paste(round(colMeans(do.call(rbind, lapply(fits, function(f) f$gp$acc_ell))), 3), collapse = ";"),
             max_rhat = max(diag$rhat), min_ess_bulk = min(diag$ess_bulk),
             min_ess_tail = min(diag$ess_tail), mc_noise = mc_noise,
             rho_star = fmt_ci("rho_star"),
             ell = paste(vapply(paste0("ell", seq_len(d)), fmt_ci, ""), collapse = "; "),
             passed = passed, stringsAsFactors = FALSE)
}
