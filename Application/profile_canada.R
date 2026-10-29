###############################################################################
# profile_canada.R -- data summary, binning accuracy and timing before the final runs.
# Run from this folder:  Rscript profile_canada.R            (Ontario; about 10-30 minutes)
#                        Rscript profile_canada.R all        (all three provinces)
# Short runs only: not inference. Output: Results/Profile/ (CSV tables + printed summary).
#
# 1. Data summary per province (pixels, events, snapped border events).
# 2. Binning accuracy: for 20 prior draws of (ell, w), with rho* set so that the expected number
#    of events equals the observed one, the binned log-likelihood is compared with the
#    unbinned one (G = NULL: every pixel of every year, every event). What matters for the
#    MH ratios is the variation of the error across states, reported as its SD (and the
#    maximum pairwise difference); a constant offset cancels.
# 3. Timing: one log-likelihood evaluation per G, and a 300-sweep chain (100 burn-in) per
#    covariate dimension, extrapolated to the default budget (60000 sweeps).
###############################################################################

# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("Matrix", "coda")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
suppressPackageStartupMessages({ library(Matrix); library(coda) })
source("../R/gp_sampler.R"); source("../R/canada_data.R")
args <- commandArgs(TRUE)
regions <- if (length(args) && args[1] == "all") names(CANADA_REGIONS) else "ontario"
out <- file.path("Results", "Profile"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261010)
prior <- list(a = 1, b = 2.5, lambda = 0.1)      # final prior (run_canada.R defaults)

loglik_fun <- function(binned, nd) {
  A_dom <- interp_matrix(nd, binned$z_dom); A_obs <- interp_matrix(nd, binned$z_obs)
  list(ll = function(w, rs) sum(binned$h_obs * plogis(as.vector(A_obs %*% w), log.p = TRUE)) -
         rs * sum(binned$h_dom * plogis(as.vector(A_dom %*% w))),
       mass = function(w) sum(binned$h_dom * plogis(as.vector(A_dom %*% w))))
}
nodes_for <- function(d) if (d == 1) make_nodes(1, V = 200) else make_nodes(d, type = "grid", grid_m = c(NA, 20, 20)[d])
G_values <- list(c(500, 2000), c(100, 200), c(50, 100, 200))

data_rows <- list(); bin_rows <- list(); time_rows <- list()
for (region in regions) {
  R <- load_canada_region(region); data_rows[[region]] <- canada_data_summary(R)
  cat("\n==== ", region, "\n"); print(data_rows[[region]], row.names = FALSE)
  for (model in c("temp", "temp_prec", "full")) {
    covs <- CANADA_MODELS[[model]]; d <- length(covs); nd <- nodes_for(d)
    reps <- canada_reps(R, covs)
    exact <- loglik_fun(bin_data(reps, d, G = NULL), nd)
    draws <- lapply(1:20, function(s) {
      ell <- r_ell_prior(d, prior$a, prior$b, prior$lambda)
      w <- fmul(make_factor(nd, NULL, ell), rnorm(nrow(nd$nodes)))
      list(w = w, rs = sum(R$n_events) / exact$mass(w))
    })
    ll_exact <- vapply(draws, function(x) exact$ll(x$w, x$rs), 0)
    t_exact <- system.time(for (k in 1:5) exact$ll(draws[[1]]$w, draws[[1]]$rs))[["elapsed"]] / 5
    for (G in G_values[[d]]) {
      b <- bin_data(reps, d, G = G); L <- loglik_fun(b, nd)
      err <- vapply(draws, function(x) L$ll(x$w, x$rs), 0) - ll_exact
      tG <- system.time(for (k in 1:20) L$ll(draws[[1]]$w, draws[[1]]$rs))[["elapsed"]] / 20
      bin_rows[[length(bin_rows) + 1]] <- data.frame(region = region, model = model, d = d,
        nodes = nrow(nd$nodes), G = G, bins_dom = nrow(b$z_dom), bins_obs = nrow(b$z_obs),
        quadrature_points = sum(vapply(reps, function(r) nrow(r$Zq), 0L)),
        err_sd = sd(err), err_max_pairwise = diff(range(err)), ll_sd_across_draws = sd(ll_exact),
        ms_per_loglik = 1000 * tG, ms_per_loglik_unbinned = 1000 * t_exact)
      cat(sprintf("  %-9s G=%-4d bins %7d  error SD %.3g (pairwise max %.3g; SD of loglik across draws %.3g)  %.2f ms/eval (unbinned %.1f ms)\n",
                  model, G, nrow(b$z_dom), sd(err), diff(range(err)), sd(ll_exact), 1000 * tG, 1000 * t_exact))
    }
  }
}
write.csv(do.call(rbind, data_rows), file.path(out, "data_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, bin_rows), file.path(out, "binning_accuracy.csv"), row.names = FALSE)

# timing of short chains (Ontario; the default configuration of the final runs)
R <- load_canada_region("ontario")
for (cfg in list(list(model = "temp", G = 2000), list(model = "temp_prec", G = 200),
                 list(model = "full", G = 100), list(model = "full", G = 50))) {
  covs <- CANADA_MODELS[[cfg$model]]; d <- length(covs); nd <- nodes_for(d)
  b <- bin_data(canada_reps(R, covs), d, G = cfg$G)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(1)
  fit <- gp_pcn_mcmc(b, nd, canada_probes(d), n_iter = 300, burn = 100, thin = 10,
                     a_ell = prior$a, b_ell = prior$b, lambda = prior$lambda,
                     save_nodes_every = c(10, 10, 50)[d], verbose = FALSE)
  time_rows[[length(time_rows) + 1]] <- data.frame(region = "ontario", model = cfg$model, d = d,
    nodes = nrow(nd$nodes), G = cfg$G, bins_dom = nrow(b$z_dom), seconds_300 = fit$time,
    projected_minutes_60000 = fit$time / 300 * 60000 / 60,
    acc_w_short = fit$acc_w, rho_star_last = tail(fit$trace$rho_star, 1),
    ell_last = paste(round(tail(fit$trace$ell, 1), 2), collapse = ";"))
  cat(sprintf("Ontario %-9s nodes %5d G %4d: %.1f s for 300 sweeps -> about %.0f min per 60000-sweep chain (one core)\n",
              cfg$model, nrow(nd$nodes), cfg$G, fit$time, fit$time / 300 * 60000 / 60))
}
write.csv(do.call(rbind, time_rows), file.path(out, "timing.csv"), row.names = FALSE)
writeLines(c(R.version.string, paste("BLAS:", extSoftVersion()[["BLAS"]]), format(Sys.time())),
           file.path(out, "session.txt"))
cat("\nWritten to", out, "\n")
