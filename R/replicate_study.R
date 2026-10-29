###############################################################################
# replicate_study.R -- fixed-cap one-chain fits for the replicated error study
# Protocol: Simulation Study/SIMULATION_PROTOCOL.md, Section 4 (version 2.2, 12 October 2026).
#
# One dataset = one worker process:
#   saved observations -> competitors (kernel with oracle and CV bandwidths, log-polynomials of
#   degree 1-4) -> one GP chain (60000 sweeps, 20000 burn-in/adaptation, thinning 10) ->
#   single-chain summaries and flags -> verified compact export -> removal of that dataset's
#   temporary chain/checkpoint only.
#
# Sourced together with data_generation.R, gp_sampler.R, competitors.R and chain_diagnostics.R
# (probe locations). Version 2.2 changes the length-scale hyper-prior to
# (a_ell, b_ell, lambda) = (1, 2.5, 0.1), the kernel bandwidth grid (competitors.R) and the
# parametric competitors (full log-polynomials of degree 1-4, previously log-linear and
# additive log-cubic, with full degree 2-4 as post-hoc add-ons) and the 2D GP grid (20 x 20
# nodes, previously 25 x 25).
###############################################################################

repstudy_science <- function(smoke = FALSE) list(
  version = "replicate_study",
  n_iter = if (smoke) 600L else 60000L,
  burn = if (smoke) 200L else 20000L,
  thin = 10L,
  n_w_steps = 10L,
  ell_update = "noncentred",
  target_w = 0.30, target_ell = 0.44,
  checkpoint_every = if (smoke) 200L else 10000L,
  prior = list(a_ell = 1, b_ell = 2.5, lambda = 0.1, a_rho = 1, b_rho = 0.1, c_rho = 100),
  poly_degrees = 1:4,
  V_1d = if (smoke) 20L else 200L,
  grid_m = if (smoke) 8L else 20L,
  m = 50L,
  seed_offset = 70000007L,          # mcmc_seed = data_seed + 70000007 (distinct from all earlier streams)
  noise_flag = 0.01,                # relative RMS posterior-mean MC noise
  split_rhat_flag = 1.05)           # within-chain split R-hat on rho(z) probes and mean intensity

# Methods on the common reporting grid, in table order.
REPSTUDY_METHODS <- c("GP", "Kernel_oracle", "Kernel_CV", paste0("Log_poly", 1:4))
REPSTUDY_LABELS <- c(GP = "GP", Kernel_oracle = "Kernel oracle", Kernel_CV = "Kernel CV",
                 Log_poly1 = "Log-poly 1", Log_poly2 = "Log-poly 2", Log_poly3 = "Log-poly 3",
                 Log_poly4 = "Log-poly 4")

repstudy_tables <- c("summary", "points", "diagnostics", "parameters", "errors", "chain",
                 "timings", "concentration", "kernel_scores")

repstudy_dir <- function(out, scenario, n, rep)
  file.path(out, scenario, sprintf("n%04d", as.integer(n)), sprintf("r%02d", as.integer(rep)))

repstudy_data_seed <- function(d, n, rep) as.integer(1000000L * d + 1000L * n + rep)

repstudy_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  p <- tempfile("pending_", tmpdir = dirname(path)); on.exit(unlink(p))
  saveRDS(x, p)
  if (!file.rename(p, path)) stop("Atomic save failed: ", path)
}

repstudy_json <- function(x, path) {
  p <- tempfile("pending_", tmpdir = dirname(path)); on.exit(unlink(p))
  jsonlite::write_json(x, p, pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 16)
  if (!file.rename(p, path)) stop("Atomic save failed: ", path)
}

# ---------------------------------------------------------------- data and competitors
repstudy_load_data <- function(cfg, scenario, n, rep) {
  sc <- experiment_scenarios()[[scenario]]
  stopifnot(!is.null(sc))
  seed <- repstudy_data_seed(sc$d, n, rep)
  key <- sprintf("%s_n%d_r%d", scenario, as.integer(n), as.integer(rep))
  reps <- get_replicates(key, n, sc$rho0, sc$ls, seed = seed, m = cfg$science$m,
                         data_dir = cfg$data_dir)
  list(reps = reps, dataset = attr(reps, "dataset"), sc = sc, seed = seed)
}

repstudy_competitors <- function(reps, sc, z, degrees = 1:4) {
  d <- sc$d
  sec <- c(kernel = system.time(kc <- kernel_competitor(reps, d, sc$rho0))[["elapsed"]])
  stopifnot(isTRUE(all.equal(kc$truth, sc$rho0(kc$centres))))
  common <- kernel_on_grid(kc, z)                     # linear (1D) / bilinear (2D) from bin centres
  colnames(common) <- c("Kernel_oracle", "Kernel_CV")
  polys <- lapply(degrees, function(k) {
    t0 <- proc.time()[["elapsed"]]; f <- poly_competitor(reps, z, k)
    f$seconds <- proc.time()[["elapsed"]] - t0; f
  })
  names(polys) <- paste0("Log_poly", degrees)
  poly_curves <- sapply(polys, `[[`, "curve")
  if (!is.matrix(poly_curves)) poly_curves <- matrix(poly_curves, ncol = length(degrees),
                                                     dimnames = list(NULL, names(polys)))
  sec <- c(sec, vapply(polys, `[[`, numeric(1), "seconds"))
  stopifnot(nrow(common) == nrow(z), all(is.finite(common)), all(is.finite(poly_curves)))
  list(kernel_common = common, poly = poly_curves, seconds = sec,
       poly_info = data.frame(method = names(polys),
         converged = vapply(polys, `[[`, logical(1), "converged"),
         iterations = vapply(polys, `[[`, integer(1), "iterations"),
         n_coef = vapply(polys, function(f) length(f$coef), integer(1)), row.names = NULL),
       bandwidths = list(oracle = kc$oracle$h, cv = kc$cv$h,
                         oracle_at_edge = kc$oracle$at_edge, cv_at_edge = kc$cv$at_edge),
       kernel_scores = as.data.frame(kc$scores),
       kernel_reporting = sprintf(paste("bandwidth selection on %d^%d bins; %s interpolation from the bin",
         "centres to the common reporting grid"), kc$G, d, if (d == 1) "linear" else "bilinear"))
}

# ---------------------------------------------------------------- one GP chain
repstudy_sample <- function(cfg, dat, nd, z, work) {
  s <- cfg$science; p <- s$prior; d <- dat$sc$d; n <- length(dat$reps)
  mseed <- as.integer(dat$seed + s$seed_offset)
  # Record the cold start (a draw from the prior); gp_pcn_mcmc draws the same values from
  # the same seed in the same order.
  set_mcmc_seed(mseed)
  ell0 <- r_ell_prior(d, p$a_ell, p$b_ell, p$lambda)
  initial <- list(ell = ell0, w = fmul(make_factor(nd, NULL, ell0), rnorm(nrow(nd$nodes))),
                  rho_star = rtgamma(p$a_rho, p$b_rho, p$c_rho + log(n)))
  set_mcmc_seed(mseed)
  gp <- do.call(gp_pcn_mcmc, c(list(binned = bin_data(dat$reps, d), nd = nd, z_eval = z,
    n_iter = s$n_iter, burn = s$burn, thin = s$thin, ell_update = s$ell_update,
    n_w_steps = s$n_w_steps, target_w = s$target_w, target_ell = s$target_ell,
    checkpoint_path = file.path(work, "checkpoint.rds"), checkpoint_every = s$checkpoint_every,
    checkpoint_guard = cfg$guard, verbose = TRUE), p))
  stopifnot(nrow(gp$draws) == (s$n_iter - s$burn) %/% s$thin, all(is.finite(gp$draws)))
  list(gp = gp, mcmc_seed = mseed, initial = initial)
}

# ---------------------------------------------------------------- single-chain summaries
repstudy_summarise <- function(cfg, scenario, n, rep, dat, z, truth, smp, comp, data_seconds) {
  t_post <- proc.time()[["elapsed"]]
  s <- cfg$science; g <- smp$gp; D <- g$draws; keep <- g$keep
  d <- ncol(z); nd <- nrow(D); p <- ncol(D); events <- n_events(dat$reps)
  stopifnot(p == nrow(z), length(truth) == p, all(is.finite(D)), all(D >= 0),
            nd == (s$n_iter - s$burn) %/% s$thin, sum(keep) == nd,
            max(abs(colMeans(D) - g$post_mean)) < 1e-10)

  # Monitored scalar/probe traces (same quantities as the validated pilots, one chain).
  probes <- pilot_probe_indices(z)
  vars <- c("rho_star", paste0("ell", seq_len(d)), "poisson_loglik", "mean_intensity",
            sprintf("rho_%04d", probes))
  a <- array(NA_real_, c(nd, 1L, length(vars)), dimnames = list(NULL, "chain1", vars))
  a[, 1, ] <- cbind(g$trace$rho_star[keep], g$trace$ell[keep, , drop = FALSE],
                    g$trace$loglik[keep] + events * log(g$trace$rho_star[keep]),
                    rowMeans(D), D[, probes, drop = FALSE])
  diagnostics <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a),
    "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  # With one chain, posterior::rhat is the rank-normalised split R-hat of the two halves.
  names(diagnostics)[names(diagnostics) == "rhat"] <- "split_rhat_within_chain"
  target <- diagnostics$variable == "mean_intensity" | grepl("^rho_[0-9]{4}$", diagnostics$variable)

  # Posterior mean and its single-chain MCSE at every reporting location.
  mu <- colMeans(D)
  mean_mcse <- vapply(seq_len(p), function(j) posterior::mcse_mean(matrix(D[, j], ncol = 1L)), numeric(1))
  variance <- apply(D, 2, var)
  interval <- vapply(seq_len(p), function(j) quantile(D[, j], c(.025, .975), names = FALSE), numeric(2))
  qse <- matrix(NA_real_, p, 2)
  for (j in probes) qse[j, ] <- posterior::mcse_quantile(matrix(D[, j], ncol = 1L), probs = c(.025, .975))

  noise <- sqrt(mean(mean_mcse^2)) / sqrt(mean(mu^2))
  flag_noise <- !is.finite(noise) || noise > s$noise_flag
  sr <- diagnostics$split_rhat_within_chain[target]
  flag_rhat <- any(!is.finite(sr) | sr > s$split_rhat_flag)

  # All methods on the common reporting grid.
  curves <- cbind(GP = mu, comp$kernel_common, comp$poly)
  stopifnot(identical(colnames(curves), REPSTUDY_METHODS),
            all(is.finite(curves)))
  errors <- data.frame(method = colnames(curves),
                       relative_L2 = unname(apply(curves, 2, rel_L2, truth = truth)))

  # Descriptive hyperparameter summaries (validated estimates come from long multichain runs).
  hn <- c("rho_star", paste0("ell", seq_len(d)), paste0("lengthscale", seq_len(d)))
  sa <- array(NA_real_, c(nd, 1L, length(hn)), dimnames = list(NULL, "chain1", hn))
  sa[, 1, 1:(1 + d)] <- a[, 1, 1:(1 + d)]
  sa[, 1, (2 + d):(1 + 2 * d)] <- 1 / (sqrt(2) * a[, 1, 2:(1 + d)])   # exp(-(z-z')^2/(2 L^2))
  parameters <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(sa),
    "mean", "sd", "ess_bulk", "ess_tail", "mcse_mean",
    ~posterior::quantile2(.x, probs = c(.025, .25, .5, .75, .975))))
  parameters$scope <- "descriptive one-chain capped summary; not a validated hyperparameter estimate"

  # Optional posterior-concentration summaries (REPORTING_PLAN_1D.md).
  den <- sqrt(mean(truth^2)); l2sq <- dzsq <- numeric(nd); esqrt <- numeric(p)
  for (block in split(seq_len(p), ceiling(seq_len(p) / 50))) {
    X <- D[, block, drop = FALSE]
    l2sq <- l2sq + rowSums(sweep(X, 2, truth[block])^2) / p
    SX <- sqrt(X); esqrt[block] <- colMeans(SX)
    dzsq <- dzsq + rowSums(sweep(SX, 2, sqrt(truth[block]))^2) / p
  }
  stopifnot(isTRUE(all.equal(mean(l2sq), mean((mu - truth)^2) + mean(variance) * (nd - 1) / nd,
                             tolerance = 1e-8)))
  concentration <- data.frame(quantity = c("relative_posterior_RMS_L2", "posterior_RMS_dZ"),
    estimate = c(sqrt(mean(l2sq)) / den, sqrt(mean(dzsq))),
    mean_square = c(mean(l2sq) / den^2, mean(dzsq)),
    mcse_delta = c(posterior::mcse_mean(matrix(l2sq, ncol = 1L)) / (2 * sqrt(mean(l2sq)) * den),
                   posterior::mcse_mean(matrix(dzsq, ncol = 1L)) / (2 * sqrt(mean(dzsq)))),
    scope = "one-chain capped run; uncertainty precision not validated")

  zz <- as.data.frame(matrix(z, ncol = d)); names(zz) <- if (d == 1) "z" else c("z1", "z2")
  points <- cbind(zz, data.frame(truth = truth, curves, posterior_sd = sqrt(variance),
    lower = interval[1, ], upper = interval[2, ], mean_mcse = mean_mcse,
    lower_mcse = qse[, 1], upper_mcse = qse[, 2], posterior_mean_sqrt = esqrt, check.names = FALSE))

  chain <- data.frame(chain = 1L, mcmc_seed = smp$mcmc_seed, n_iter = s$n_iter, burn = s$burn,
    thin = s$thin, draws = nd, sampling_seconds = g$time, acceptance_w = g$acc_w,
    beta = g$beta, rho_star_start = smp$initial$rho_star)
  for (j in seq_len(d)) {
    chain[[paste0("acceptance_ell", j)]] <- g$acc_ell[j]
    chain[[paste0("s_ell", j)]] <- g$s_ell[j]
    chain[[paste0("ell_start", j)]] <- smp$initial$ell[j]
  }

  ell_rows <- grepl("^ell[0-9]+$", diagnostics$variable)
  summary <- data.frame(scenario = scenario, n = as.integer(n), rep = as.integer(rep),
    data_id = dat$dataset$data_id, data_seed = dat$seed, events = events, chains = 1L,
    n_iter = s$n_iter, burn = s$burn, thin = s$thin, draws = nd, mcmc_seed = smp$mcmc_seed,
    gp_relative_L2 = errors$relative_L2[1], relative_RMS_MC_noise = noise,
    max_split_rhat_target = max(sr), min_bulk_ess_target = min(diagnostics$ess_bulk[target]),
    min_tail_ess_target = min(diagnostics$ess_tail[target]),
    min_bulk_ess_ell = min(diagnostics$ess_bulk[ell_rows]),
    bulk_ess_rho_star = diagnostics$ess_bulk[diagnostics$variable == "rho_star"],
    between_chain_rhat = "unavailable: one chain",
    flag_noise = flag_noise, flag_split_rhat = flag_rhat,
    flag_for_inspection = flag_noise || flag_rhat,
    sampling_minutes = g$time / 60, competitor_minutes = sum(comp$seconds) / 60,
    smoke = cfg$smoke, protocol = s$version)

  timings <- data.frame(
    component = c("data_load", "kernel_oracle_and_CV", tolower(names(comp$seconds)[-1]),
                  "gp_sampling", "gp_postprocessing"),
    seconds = c(data_seconds, unname(comp$seconds), g$time, proc.time()[["elapsed"]] - t_post),
    scope = c("load or generate saved observations",
              "one kernel computation selects both bandwidths over the full grid",
              rep("exact pixel-quadrature Poisson ML", length(comp$seconds) - 1),
              "cumulative chain time including resumed sessions",
              "single-chain summaries before export"))

  list(schema = "replicate_study",
       meta = list(scenario = scenario, n = as.integer(n), rep = as.integer(rep), d = d,
                   dataset = dat$dataset, settings = s, smoke = cfg$smoke,
                   sources = cfg$guard$sources, R = R.version.string,
                   BLAS = extSoftVersion()[["BLAS"]], kernel_reporting = comp$kernel_reporting),
       summary = summary, points = points, diagnostics = diagnostics, parameters = parameters,
       errors = errors, chain = chain, timings = timings, concentration = concentration,
       kernel_scores = comp$kernel_scores, bandwidths = comp$bandwidths, poly_info = comp$poly_info,
       monitored_draws = a, parameter_draws = sa, monitored_indices = probes,
       retained_iterations = s$burn + s$thin * seq_len(nd), initial = smp$initial)
}

# ---------------------------------------------------------------- figures and export
repstudy_figures <- function(x, dir) {
  cols <- c("#3274A1", "#DD8452", "#DD8452", "#8C8C8C", "#8172B3", "#55A868", "#C44E52")
  a <- x$points; d <- x$meta$d; methods <- x$errors$method
  fit <- function() {
    if (d == 1) {
      par(mar = c(4, 4, 3, 1))
      plot(a$z, a$truth, type = "n", ylim = range(0, as.matrix(a[, c("truth", methods)])),
           xlab = "Covariate z", ylab = "Intensity",
           main = sprintf("%s n=%d, dataset %d (one chain, %d sweeps)", x$summary$scenario,
                          x$summary$n, x$summary$rep, x$summary$n_iter))
      lines(a$z, a$truth, lwd = 2)
      matlines(a$z, as.matrix(a[, methods]), col = cols, lty = c(1, 1, 2, 3, 3, 3, 3), lwd = 1.5)
      legend("topleft", c("Truth", methods), col = c("black", cols), lty = c(1, 1, 1, 2, 3, 3, 3, 3),
             bty = "n", cex = .75)
    } else {
      xx <- sort(unique(a$z1)); yy <- sort(unique(a$z2))
      vals <- as.matrix(a[, c("truth", methods)]); lim <- range(0, vals)
      pal <- hcl.colors(100, "YlOrRd", rev = TRUE)
      par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
      for (k in seq_len(ncol(vals))) {
        M <- matrix(NA_real_, length(xx), length(yy))
        M[cbind(match(a$z1, xx), match(a$z2, yy))] <- vals[, k]
        image(xx, yy, M, zlim = lim, col = pal, xlab = "z1", ylab = "z2", asp = 1,
              useRaster = TRUE, main = gsub("_", " ", colnames(vals)[k]))
      }
    }
  }
  trace <- function() {
    v <- dimnames(x$monitored_draws)[[3]]
    show <- c("rho_star", grep("^ell[0-9]+$", v, value = TRUE), "poisson_loglik", "mean_intensity")
    par(mfrow = c(length(show), 1), mar = c(3, 4, 2, 1))
    for (nm in show) plot(x$retained_iterations, x$monitored_draws[, 1, nm], type = "l",
                          col = cols[1], xlab = "Sweep (post burn-in)", ylab = nm, main = nm)
  }
  for (kind in c("fit", "trace")) {
    draw <- if (kind == "fit") fit else trace
    h <- if (kind == "fit") (if (d == 1) 6 else 8) else 2.4 * (3 + d)
    w <- if (kind == "fit" && d == 2) 15 else 9
    pdf(file.path(dir, paste0(kind, ".pdf")), width = w, height = h); draw(); dev.off()
    png(file.path(dir, paste0(kind, ".png")), width = 150 * w, height = 150 * h, res = 150)
    draw(); dev.off()
  }
}

# CSV round trip: read everything as text, then compare by column type.
repstudy_same_table <- function(a, path, tol = 1e-10) {
  b <- read.csv(path, check.names = FALSE, colClasses = "character", na.strings = "NA")
  if (nrow(a) != nrow(b) || !identical(names(a), names(b))) return(FALSE)
  for (v in names(a)) {
    x <- a[[v]]; y <- b[[v]]
    if (is.logical(x)) {
      if (!identical(x, as.logical(y))) return(FALSE)
    } else if (is.numeric(x)) {
      x <- as.numeric(x); y <- suppressWarnings(as.numeric(y))
      if (!identical(is.na(x), is.na(y))) return(FALSE)
      ok <- !is.na(x)
      same <- (x[ok] == y[ok]) |
        (is.finite(x[ok]) & abs(x[ok] - y[ok]) <= tol * pmax(1, abs(x[ok])))
      if (!all(same %in% TRUE)) return(FALSE)
    } else if (!identical(as.character(x), y)) return(FALSE)
  }
  TRUE
}

repstudy_verify <- function(dir, expected = NULL) {
  x <- readRDS(file.path(dir, "compact_fit.rds"))
  if (!is.null(expected) && !identical(x, expected)) stop("Compact RDS round trip failed: ", dir)
  for (nm in repstudy_tables)
    if (!repstudy_same_table(x[[nm]], file.path(dir, paste0(nm, ".csv"))))
      stop("CSV round trip failed: ", file.path(dir, paste0(nm, ".csv")))
  figs <- file.path(dir, c("fit.pdf", "fit.png", "trace.pdf", "trace.png"))
  stopifnot(nrow(x$points) == nrow(eval_grid(x$meta$d)), x$summary$chains == 1L,
            dim(x$monitored_draws)[1] == x$summary$draws,
            isTRUE(all.equal(x$summary$gp_relative_L2, rel_L2(x$points$GP, x$points$truth),
                             tolerance = 1e-12)),
            all(file.exists(figs)), all(file.info(figs)$size > 100))
  invisible(x)
}

repstudy_write <- function(x, dir) {
  repstudy_atomic(x, file.path(dir, "compact_fit.rds"))
  for (nm in repstudy_tables) write.csv(x[[nm]], file.path(dir, paste0(nm, ".csv")), row.names = FALSE)
  repstudy_figures(x, dir)
  repstudy_verify(dir, x)
  files <- file.path(dir, c("compact_fit.rds", paste0(repstudy_tables, ".csv"),
                            "fit.pdf", "fit.png", "trace.pdf", "trace.png"))
  repstudy_json(list(verified = TRUE, verified_at = as.character(Sys.time()),
                 data_id = x$summary$data_id,
                 files = as.list(setNames(unname(tools::md5sum(files)), basename(files)))),
            file.path(dir, "verification.json"))
  invisible(x)
}

repstudy_check_manifest <- function(dir) {
  v <- jsonlite::read_json(file.path(dir, "verification.json"), simplifyVector = TRUE)
  expected <- unlist(v$files)
  actual <- unname(tools::md5sum(file.path(dir, names(expected))))
  if (!isTRUE(v$verified) || !identical(unname(expected), actual))
    stop("Compact record changed after verification: ", dir)
  repstudy_verify(dir)
}

# Remove only this dataset's registered temporary files, after verified export.
repstudy_cleanup <- function(dir) {
  x <- repstudy_check_manifest(dir)
  work <- file.path(dir, "Work")
  disposable <- file.path(work, c("chain.rds", "checkpoint.rds"))
  disposable <- disposable[file.exists(disposable)]
  if (length(disposable)) stopifnot(all(file.remove(disposable)))
  files <- list.files(dir, full.names = TRUE); files <- files[!file.info(files)$isdir]
  repstudy_json(list(status = "complete", data_id = x$summary$data_id, verified = TRUE,
                 flag_for_inspection = x$summary$flag_for_inspection,
                 full_grid_draws_retained = FALSE, checkpoint_retained = FALSE,
                 compact_bytes = sum(file.info(files)$size, na.rm = TRUE),
                 completed_at = as.character(Sys.time())),
            file.path(dir, "final_status.json"))
  invisible(x)
}

# ---------------------------------------------------------------- worker
repstudy_worker <- function(cfg, scenario, n, rep) {
  dir <- repstudy_dir(cfg$out, scenario, n, rep); work <- file.path(dir, "Work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  # Per-dataset lock: a worker orphaned by a killed controller must not share this directory.
  lock <- file.path(work, ".worker_lock")
  if (dir.exists(lock)) {
    pf <- file.path(lock, "pid")
    pid <- if (file.exists(pf)) suppressWarnings(as.integer(readLines(pf, warn = FALSE))) else NA_integer_
    if (length(pid) == 1 && !is.na(pid) && isTRUE(tryCatch(ps::ps_is_running(ps::ps_handle(as.integer(pid))), error = function(e) FALSE)))
      stop("Another worker (pid ", pid, ") is still running on ", dir)
    unlink(lock, recursive = TRUE)
  }
  if (!dir.create(lock)) stop("Could not create worker lock: ", lock)
  writeLines(as.character(Sys.getpid()), file.path(lock, "pid"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(file.path(dir, "final_status.json"))) {
    repstudy_check_manifest(dir); cat("Already complete:", dir, "\n"); return(invisible("complete"))
  }
  if (file.exists(file.path(dir, "verification.json"))) {    # interrupted during cleanup
    repstudy_cleanup(dir); cat("Finalised:", dir, "\n"); return(invisible("complete"))
  }
  s <- cfg$science
  data_seconds <- system.time(dat <- repstudy_load_data(cfg, scenario, n, rep))[["elapsed"]]
  d <- dat$sc$d; z <- eval_grid(d); truth <- dat$sc$rho0(z)
  nd <- if (d == 1) make_nodes(1, V = s$V_1d) else make_nodes(2, type = "grid", grid_m = s$grid_m)
  cat(sprintf("%s n=%d r=%d: data_id %s, %d events\n", scenario, as.integer(n), as.integer(rep),
              dat$dataset$data_id, n_events(dat$reps))); flush.console()

  comp_path <- file.path(work, "competitors.rds")
  if (!file.exists(comp_path))
    repstudy_atomic(c(list(data_id = dat$dataset$data_id), repstudy_competitors(dat$reps, dat$sc, z, s$poly_degrees)), comp_path)
  comp <- readRDS(comp_path)
  stopifnot(identical(comp$data_id, dat$dataset$data_id))

  chain_path <- file.path(work, "chain.rds")
  if (!file.exists(chain_path)) {
    smp <- repstudy_sample(cfg, dat, nd, z, work)
    repstudy_atomic(c(smp, list(data_id = dat$dataset$data_id, settings = s)), chain_path)
  }
  smp <- readRDS(chain_path)
  stopifnot(identical(smp$data_id, dat$dataset$data_id), identical(smp$settings, s))
  cp_path <- file.path(work, "checkpoint.rds")
  if (file.exists(cp_path)) {
    cp <- readRDS(cp_path)
    stopifnot(cp$iteration == s$n_iter, identical(cp$draws, smp$gp$draws),
              identical(cp$contract$external_guard, cfg$guard))
    rm(cp); gc(FALSE)
  }
  x <- repstudy_summarise(cfg, scenario, n, rep, dat, z, truth, smp, comp, data_seconds)
  rm(smp); gc(FALSE)
  repstudy_write(x, dir)
  repstudy_cleanup(dir)
  cat(sprintf("Complete %s n=%d r=%d: GP relative L2 %.5f, MC noise %.4f%%, flagged %s\n",
              scenario, as.integer(n), as.integer(rep), x$summary$gp_relative_L2,
              100 * x$summary$relative_RMS_MC_noise, x$summary$flag_for_inspection))
  invisible("complete")
}

# ---------------------------------------------------------------- aggregate tables (no sampling)
repstudy_report <- function(cfg) {
  status_files <- list.files(cfg$out, pattern = "^final_status\\.json$", recursive = TRUE,
                             full.names = TRUE)
  if (!length(status_files)) { cat("No completed datasets yet.\n"); return(invisible(NULL)) }
  xs <- lapply(dirname(status_files), repstudy_check_manifest)
  tab <- file.path(cfg$out, "Tables"); dir.create(tab, showWarnings = FALSE)
  by_cell <- function(df, f) {
    parts <- split(df, list(df$scenario, df$n), drop = TRUE)
    out <- do.call(rbind, lapply(parts, f)); rownames(out) <- NULL
    out[order(out$scenario, out$n), , drop = FALSE]
  }
  msd <- function(v) c(mean = mean(v), sd = if (length(v) > 1) sd(v) else NA_real_,
                       mcse = if (length(v) > 1) sd(v) / sqrt(length(v)) else NA_real_)

  summ <- do.call(rbind, lapply(xs, `[[`, "summary"))
  summ <- summ[order(summ$scenario, summ$n, summ$rep), ]; rownames(summ) <- NULL
  write.csv(summ, file.path(tab, "dataset_summary.csv"), row.names = FALSE)

  err <- do.call(rbind, lapply(xs, function(x)
    cbind(x$summary[rep(1L, nrow(x$errors)), c("scenario", "n", "rep", "data_id")], x$errors)))
  rownames(err) <- NULL
  write.csv(err, file.path(tab, "method_errors.csv"), row.names = FALSE)

  methods <- REPSTUDY_METHODS
  error_table <- by_cell(err, function(e) do.call(rbind, lapply(methods, function(m) {
    v <- e$relative_L2[e$method == m]; s <- msd(v)
    data.frame(scenario = e$scenario[1], n = e$n[1], method = m, datasets = length(v),
               mean_relative_L2 = s[["mean"]], sd_relative_L2 = s[["sd"]], simulation_MCSE = s[["mcse"]])
  })))
  write.csv(error_table, file.path(tab, "error_table.csv"), row.names = FALSE)

  gp <- err[err$method == "GP", c("scenario", "n", "rep", "relative_L2")]
  names(gp)[4] <- "gp_relative_L2"
  pd <- merge(err[err$method != "GP", ], gp, by = c("scenario", "n", "rep"))
  pd$difference <- pd$relative_L2 - pd$gp_relative_L2
  paired <- by_cell(pd, function(e) do.call(rbind, lapply(setdiff(methods, "GP"), function(m) {
    v <- e$difference[e$method == m]; s <- msd(v)
    data.frame(scenario = e$scenario[1], n = e$n[1], method = m, datasets = length(v),
               mean_difference_minus_GP = s[["mean"]], sd_difference = s[["sd"]],
               simulation_MCSE = s[["mcse"]], fraction_GP_smaller_error = mean(v > 0))
  })))
  write.csv(paired, file.path(tab, "paired_differences_vs_GP.csv"), row.names = FALSE)

  computation <- by_cell(summ, function(e) data.frame(scenario = e$scenario[1], n = e$n[1],
    datasets = nrow(e), events_median = median(e$events),
    sampling_minutes_median = median(e$sampling_minutes), sampling_minutes_max = max(e$sampling_minutes),
    competitor_minutes_median = median(e$competitor_minutes),
    MC_noise_median = median(e$relative_RMS_MC_noise), MC_noise_max = max(e$relative_RMS_MC_noise),
    max_split_rhat_target = max(e$max_split_rhat_target),
    min_bulk_ess_target = min(e$min_bulk_ess_target),
    bulk_ess_ell_median = median(e$min_bulk_ess_ell), bulk_ess_ell_min = min(e$min_bulk_ess_ell),
    flagged_noise = sum(e$flag_noise), flagged_split_rhat = sum(e$flag_split_rhat),
    flagged_any = sum(e$flag_for_inspection)))
  chains <- do.call(rbind, lapply(xs, function(x) data.frame(x$summary[c("scenario", "n", "rep")],
    acceptance_w = x$chain$acceptance_w, acceptance_ell1 = x$chain$acceptance_ell1)))
  acc <- by_cell(chains, function(e) data.frame(scenario = e$scenario[1], n = e$n[1],
    acceptance_w_median = median(e$acceptance_w), acceptance_ell1_median = median(e$acceptance_ell1)))
  computation <- merge(computation, acc, by = c("scenario", "n"))
  write.csv(computation, file.path(tab, "computation_table.csv"), row.names = FALSE)

  par_rows <- do.call(rbind, lapply(xs, function(x)
    cbind(x$summary[rep(1L, nrow(x$parameters)), c("scenario", "n", "rep")], x$parameters)))
  hyper <- by_cell(par_rows, function(e) do.call(rbind, lapply(split(e, e$variable), function(v)
    data.frame(scenario = v$scenario[1], n = v$n[1], parameter = v$variable[1], datasets = nrow(v),
               posterior_mean_median = median(v$mean),
               posterior_mean_q25 = unname(quantile(v$mean, .25)),
               posterior_mean_q75 = unname(quantile(v$mean, .75)),
               bulk_ess_median = median(v$ess_bulk), bulk_ess_min = min(v$ess_bulk),
               scope = "descriptive one-chain capped summaries"))))
  write.csv(hyper, file.path(tab, "hyperparameter_descriptive.csv"), row.names = FALSE)

  log_path <- file.path(cfg$out, "job_log.csv")
  failed <- if (file.exists(log_path)) {
    jl <- read.csv(log_path)
    jl <- jl[!duplicated(paste(jl$scenario, jl$n, jl$rep), fromLast = TRUE), ]
    jl <- jl[jl$status != "complete", ]
    jl[!file.exists(file.path(repstudy_dir(cfg$out, jl$scenario, jl$n, jl$rep), "final_status.json")), ]
  } else NULL
  if (!is.null(failed) && nrow(failed))
    write.csv(failed, file.path(tab, "unresolved_technical_failures.csv"), row.names = FALSE)

  lines <- c("# Replicate study: fixed-cap one-chain fits", "",
    sprintf("Updated %s. Protocol %s: one chain, %d sweeps, %d burn-in, thinning %d. %s",
            format(Sys.time(), "%Y-%m-%d %H:%M"), cfg$science$version, cfg$science$n_iter,
            cfg$science$burn, cfg$science$thin,
            if (cfg$smoke) "SMOKE TEST: tiny settings, not scientific output." else ""), "",
    paste0("| Scenario | n | Datasets | Flagged | ", paste(REPSTUDY_LABELS[methods], collapse = " | "), " |"),
    paste0("|---|---:|---:|---:|", strrep("---:|", length(methods))))
  for (k in seq_len(nrow(computation))) {
    cc <- computation[k, ]
    e <- error_table[error_table$scenario == cc$scenario & error_table$n == cc$n, ]
    cell <- vapply(methods, function(m) {
      r <- e[e$method == m, ]
      sprintf("%.4f (%s)", r$mean_relative_L2, if (is.na(r$sd_relative_L2)) "-" else sprintf("%.4f", r$sd_relative_L2))
    }, character(1))
    lines <- c(lines, sprintf("| %s | %d | %d | %d | %s |", cc$scenario, cc$n, cc$datasets,
                              cc$flagged_any, paste(cell, collapse = " | ")))
  }
  lines <- c(lines, "", "Cells: mean relative L2 error over datasets (between-dataset SD).",
    "Flags (relative RMS posterior-mean MC noise > 1%, or within-chain split R-hat > 1.05 on rho(z) probes/mean intensity) mark datasets for inspection; they are never excluded.",
    sprintf("Unresolved technical failures: %d%s", if (is.null(failed)) 0L else nrow(failed),
            if (!is.null(failed) && nrow(failed)) " (see Tables/unresolved_technical_failures.csv)" else ""),
    "", "Tables: error_table.csv, paired_differences_vs_GP.csv, computation_table.csv, hyperparameter_descriptive.csv, dataset_summary.csv, method_errors.csv.")
  writeLines(lines, file.path(cfg$out, "STATUS.md"))
  cat(lines, sep = "\n"); cat("\n")
  invisible(list(summary = summ, errors = err, error_table = error_table))
}
