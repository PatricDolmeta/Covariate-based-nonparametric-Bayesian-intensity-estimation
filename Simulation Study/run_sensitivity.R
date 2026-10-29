###############################################################################
# run_sensitivity.R -- sensitivity study (Supplement S4.5)
# Specification: SENSITIVITY_SPEC.md (version 3, 12 October 2026: new baseline prior
# (a_ell, b_ell, lambda) = (1, 2.5, 0.1); blocks S5 and S6 re-centred on it).
#
# One factor at a time around the replicate-study baseline: one chain, 60000 sweeps, 20000 burn-in,
# thinning 10, mcmc_seed = data_seed + 70000007 (common random numbers with the replicate study).
# The baseline of S1, S2, S4-S6 is the replicate-study fit of the same dataset (read, never rerun);
# S3 fits its own baseline (q = 50) on new data simulated at 200 x 200 pixels.
#
# Run from 'Data Analysis/Simulation Study':
#   Rscript run_sensitivity.R --list                        # blocks and settings
#   Rscript run_sensitivity.R --smoke                       # workflow check, tiny budget
#   Rscript run_sensitivity.R --blocks=all --workers=3      # everything (280 fits)
#   Rscript run_sensitivity.R --blocks=S4a                  # one block
#   Rscript run_sensitivity.R --blocks=S1,S2 --scenarios=1d_skn --n=50
#   Rscript run_sensitivity.R --settings=S6_a1b4 --reps=1:2  # single settings
#   Rscript run_sensitivity.R --report                      # tables only, no sampling
#
# Options (defaults in brackets)
#   --blocks=all | S1,S2,S3,S4a,S4b,S5,S6          [all]
#   --settings=ID[,ID...]    run exactly these settings (overrides --blocks), see --list
#   --scenarios=1d_skn,2d_hole                     [both]
#   --n=50,1000                                    [50,1000]
#   --reps=1:5 | 1,3,7:9                           [1:5]
#   --workers=2                                    [2]  fits run in parallel, one R process each
#   --out=DIR      [Results/Sensitivity; smoke: Results/Sensitivity_Smoke]
#
# Interrupting is safe: rerun the same command. Completed fits are skipped and running chains
# resume from their last 10000-sweep checkpoint. All blocks share one output directory, bound
# to the sources and settings recorded in config.rds; after a code or settings change use a
# new --out directory.
#
# From RStudio (working directory = this folder):
#   SENS_ARGS <- c("--blocks=S4a", "--workers=2"); source("run_sensitivity.R")
# (the fits still run as separate Rscript processes; the console stays busy until the batch ends).
###############################################################################
# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("Matrix", "coda", "sn", "jsonlite", "posterior", "processx", "ps")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })

SENS_SOURCES <- c("../R/data_generation.R", "../R/gp_sampler.R", "../R/replicate_study.R", "run_sensitivity.R")
for (f in setdiff(SENS_SOURCES, "run_sensitivity.R")) source(f)

SENS_BLOCKS <- c("S1", "S2", "S3", "S4a", "S4b", "S5", "S6")
REPSTUDY_DIR <- "Results/Replicate_Study"
or_default <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------- settings (one factor each)
sens_settings <- function() list(
  S0_base    = list(block = "S0",  label = "baseline refit (smoke test only)"),
  S1_coarse  = list(block = "S1",  label = "GP nodes: V = 100 (1D), 15x15 grid (2D)", V_1d = 100L, grid_m = 15L),
  S1_fine    = list(block = "S1",  label = "GP nodes: V = 400 (1D), 35x35 grid (2D)", V_1d = 400L, grid_m = 35L),
  S2_coarse  = list(block = "S2",  label = "bins: 500 (1D), 100x100 (2D)", G1 = 500L, G2 = 100L),
  S2_fine    = list(block = "S2",  label = "bins: 8000 (1D), 400x400 (2D)", G1 = 8000L, G2 = 400L),
  S3_q25     = list(block = "S3",  label = "quadrature 25x25, data at 200x200", m = 200L, quad_m = 25L),
  S3_q50     = list(block = "S3",  label = "quadrature 50x50, data at 200x200 (S3 baseline)", m = 200L, quad_m = 50L),
  S3_q100    = list(block = "S3",  label = "quadrature 100x100, data at 200x200", m = 200L, quad_m = 100L),
  S4a_cInf   = list(block = "S4a", label = "rho* not truncated (c = Inf)", prior = list(c_rho = Inf)),
  S4b_b0p01  = list(block = "S4b", label = "rho* prior rate b = 0.01", prior = list(b_rho = 0.01)),
  S4b_b2     = list(block = "S4b", label = "rho* prior rate b = 2", prior = list(b_rho = 2)),
  S5_l0p05   = list(block = "S5",  label = "tilt lambda = 0.05", prior = list(lambda = 0.05)),
  S5_l0p25   = list(block = "S5",  label = "tilt lambda = 0.25", prior = list(lambda = 0.25)),
  S6_a1b1p5  = list(block = "S6",  label = "(a_ell, b_ell) = (1, 1.5): heavier tail", prior = list(a_ell = 1, b_ell = 1.5)),
  S6_a1b4    = list(block = "S6",  label = "(a_ell, b_ell) = (1, 4): lighter tail", prior = list(a_ell = 1, b_ell = 4)))

sens_science <- function(smoke) {
  s <- repstudy_science(smoke)             # budget, sampler, baseline prior and grids of the replicate study
  s$version <- "sensitivity"
  s$settings <- sens_settings()
  s$fine_seed_offset <- 5000000L       # data seed of the m = 200 datasets (S3)
  s
}

# Numerical and prior settings of one fit: replicate-study baseline with one factor changed.
sens_resolve <- function(science, id, d) {
  def <- science$settings[[id]]
  if (is.null(def)) stop("Unknown setting: ", id)
  prior <- science$prior
  if (!is.null(def$prior)) prior[names(def$prior)] <- def$prior
  list(id = id, block = def$block,
       V_1d = or_default(def$V_1d, science$V_1d), grid_m = or_default(def$grid_m, science$grid_m),
       G = if (d == 1) or_default(def$G1, 2000L) else or_default(def$G2, 200L),
       m = or_default(def$m, science$m), quad_m = or_default(def$quad_m, science$m), prior = prior)
}

sens_job_dir <- function(out, block, setting, scenario, n, rep)
  file.path(out, block, setting, scenario, sprintf("n%04d", as.integer(n)), sprintf("r%02d", as.integer(rep)))

# ---------------------------------------------------------------- data
sens_data_key <- function(cfg, scenario, n, rep, m) {
  sc <- experiment_scenarios()[[scenario]]
  seed <- repstudy_data_seed(sc$d, n, rep)
  if (m == cfg$science$m) return(list(key = sprintf("%s_n%d_r%d", scenario, as.integer(n), as.integer(rep)),
                                      seed = seed, sc = sc))
  list(key = sprintf("sens_%s_m%d_n%d_r%d", scenario, as.integer(m), as.integer(n), as.integer(rep)),
       seed = as.integer(cfg$science$fine_seed_offset + seed), sc = sc)
}

sens_load_data <- function(cfg, scenario, n, rep, r) {
  k <- sens_data_key(cfg, scenario, n, rep, r$m)
  if (!cfg$smoke && r$m == cfg$science$m && !file.exists(file.path(cfg$data_dir, paste0(k$key, ".rds"))))
    stop("Saved replicate-study dataset missing: ", k$key)
  reps <- get_replicates(k$key, n, k$sc$rho0, k$sc$ls, seed = k$seed, m = r$m, quad_m = r$quad_m,
                         data_dir = cfg$data_dir)
  dataset <- attr(reps, "dataset")
  for (i in seq_along(reps)) reps[[i]]$Zpix <- NULL      # not used by the GP sampler
  list(reps = reps, dataset = dataset, sc = k$sc, seed = k$seed)
}

# ---------------------------------------------------------------- one GP chain
sens_fit <- function(cfg, dat, r, z, work, guard) {
  s <- cfg$science; p <- r$prior; d <- dat$sc$d; n <- length(dat$reps)
  nd <- if (d == 1) make_nodes(1, V = r$V_1d) else make_nodes(2, type = "grid", grid_m = r$grid_m)
  mseed <- as.integer(dat$seed + s$seed_offset)
  # Record the cold start (a draw from the prior), as in the replicate study; gp_pcn_mcmc draws the
  # same values from the same seed in the same order.
  set_mcmc_seed(mseed)
  ell0 <- r_ell_prior(d, p$a_ell, p$b_ell, p$lambda)
  initial <- list(ell = ell0, w = fmul(make_factor(nd, NULL, ell0), rnorm(nrow(nd$nodes))),
                  rho_star = rtgamma(p$a_rho, p$b_rho, p$c_rho + log(n)))
  set_mcmc_seed(mseed)
  gp <- do.call(gp_pcn_mcmc, c(list(binned = bin_data(dat$reps, d, r$G), nd = nd, z_eval = z,
    n_iter = s$n_iter, burn = s$burn, thin = s$thin, ell_update = s$ell_update,
    n_w_steps = s$n_w_steps, target_w = s$target_w, target_ell = s$target_ell,
    checkpoint_path = file.path(work, "checkpoint.rds"), checkpoint_every = s$checkpoint_every,
    checkpoint_guard = guard, verbose = TRUE), p))
  stopifnot(nrow(gp$draws) == (s$n_iter - s$burn) %/% s$thin, all(is.finite(gp$draws)))
  list(gp = gp, mcmc_seed = mseed, initial = initial, nodes = nrow(nd$nodes))
}

sens_summarise <- function(cfg, scenario, rep, r, dat, z, truth, smp) {
  g <- smp$gp; D <- g$draws; keep <- g$keep; d <- dat$sc$d; n <- length(dat$reps)
  mu <- colMeans(D)
  mcse <- vapply(seq_len(ncol(D)), function(j) posterior::mcse_mean(matrix(D[, j], ncol = 1L)), numeric(1))
  noise <- sqrt(mean(mcse^2)) / sqrt(mean(mu^2))
  hn <- c("rho_star", paste0("ell", seq_len(d)))
  H <- cbind(g$trace$rho_star[keep], g$trace$ell[keep, , drop = FALSE])
  a <- array(H, c(nrow(H), 1L, ncol(H)), dimnames = list(NULL, "chain1", hn))
  par <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a), "mean", "sd",
    ~posterior::quantile2(.x, probs = c(.025, .5, .975)), "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  names(par)[names(par) == "rhat"] <- "split_rhat_within_chain"   # one chain: split R-hat
  pv <- function(v, col) par[[col]][par$variable == v]
  upper <- r$prior$c_rho + log(n)
  summary <- data.frame(block = r$block, setting = r$id, scenario = scenario, n = as.integer(n),
    rep = as.integer(rep), data_key = dat$dataset$key, data_id = dat$dataset$data_id,
    data_seed = dat$seed, mcmc_seed = smp$mcmc_seed, events = n_events(dat$reps),
    gp_nodes = smp$nodes, bins_per_axis = r$G, data_m = r$m, quad_m = r$quad_m,
    a_ell = r$prior$a_ell, b_ell = r$prior$b_ell, lambda = r$prior$lambda,
    a_rho = r$prior$a_rho, b_rho = r$prior$b_rho, c_rho = r$prior$c_rho, rho_star_upper = upper,
    gp_relative_L2 = rel_L2(mu, truth), relative_RMS_MC_noise = noise,
    rho_star_mean = pv("rho_star", "mean"), rho_star_q2.5 = pv("rho_star", "q2.5"),
    rho_star_q97.5 = pv("rho_star", "q97.5"),
    rho_star_frac_near_bound = if (is.finite(upper)) mean(H[, 1] >= 0.99 * upper) else 0,
    ell1_mean = pv("ell1", "mean"), ell2_mean = if (d == 2) pv("ell2", "mean") else NA_real_,
    min_bulk_ess_hyper = min(par$ess_bulk), max_split_rhat_hyper = max(par$split_rhat_within_chain),
    acceptance_w = g$acc_w, sampling_minutes = g$time / 60,
    n_iter = cfg$science$n_iter, burn = cfg$science$burn, thin = cfg$science$thin,
    smoke = cfg$smoke, protocol = cfg$science$version, stringsAsFactors = FALSE)
  list(summary = summary, z = z, truth = truth, mean = mu, mean_mcse = mcse, parameters = par,
       resolved = r, initial = smp$initial, sources = cfg$guard$sources)
}

# ---------------------------------------------------------------- worker (one fit)
sens_worker <- function(cfg, setting, scenario, n, rep) {
  sc <- experiment_scenarios()[[scenario]]
  r <- sens_resolve(cfg$science, setting, sc$d)
  dir <- sens_job_dir(cfg$out, r$block, setting, scenario, n, rep); work <- file.path(dir, "Work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(work, ".worker_lock")          # an orphaned worker must not share this folder
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
  if (file.exists(file.path(dir, "summary.csv"))) { cat("Already complete:", dir, "\n"); return(invisible()) }

  dat <- sens_load_data(cfg, scenario, n, rep, r)
  z <- eval_grid(sc$d); truth <- sc$rho0(z)
  cat(sprintf("%s %s n=%d r=%d: data %s (%s), %d events\n", setting, scenario, as.integer(n),
              as.integer(rep), dat$dataset$key, dat$dataset$data_id, n_events(dat$reps))); flush.console()
  guard <- list(config = cfg$guard, setting = setting, scenario = scenario,
                n = as.integer(n), rep = as.integer(rep))
  chain_path <- file.path(work, "chain.rds")
  if (!file.exists(chain_path)) {
    smp <- sens_fit(cfg, dat, r, z, work, guard)
    repstudy_atomic(c(smp, list(data_id = dat$dataset$data_id, resolved = r)), chain_path)
  }
  smp <- readRDS(chain_path)
  stopifnot(identical(smp$data_id, dat$dataset$data_id), identical(smp$resolved, r))
  x <- sens_summarise(cfg, scenario, rep, r, dat, z, truth, smp)
  rm(smp); gc(FALSE)
  repstudy_atomic(x, file.path(dir, "result.rds"))
  write.csv(x$parameters, file.path(dir, "parameters.csv"), row.names = FALSE)
  tmp <- tempfile("pending_", tmpdir = dir)                # summary.csv last: marks completion
  write.csv(x$summary, tmp, row.names = FALSE)
  if (!file.rename(tmp, file.path(dir, "summary.csv"))) stop("Could not save summary.csv in ", dir)
  unlink(file.path(work, c("chain.rds", "checkpoint.rds")))
  cat(sprintf("Complete %s %s n=%d r=%d: GP relative L2 %.5f, MC noise %.3f%%, %.1f min\n", setting,
              scenario, as.integer(n), as.integer(rep), x$summary$gp_relative_L2,
              100 * x$summary$relative_RMS_MC_noise, x$summary$sampling_minutes))
  invisible()
}

# ---------------------------------------------------------------- report (no sampling)
sens_baseline <- function(cfg, x) {
  s <- x$summary
  bid <- if (s$block == "S3") "S3_q50" else if (cfg$smoke) "S0_base" else NA_character_
  if (!is.na(bid)) {
    f <- file.path(sens_job_dir(cfg$out, cfg$science$settings[[bid]]$block, bid, s$scenario, s$n, s$rep),
                   "result.rds")
    if (!file.exists(f)) return(NULL)
    b <- readRDS(f)$summary; m <- readRDS(f)$mean
    return(list(source = bid, mean = m, err = b$gp_relative_L2, noise = b$relative_RMS_MC_noise,
                data_id = b$data_id, rho_star = b$rho_star_mean, ell1 = b$ell1_mean, ell2 = b$ell2_mean))
  }
  pd <- repstudy_dir(cfg$repstudy_dir, s$scenario, s$n, s$rep)
  if (!file.exists(file.path(pd, "final_status.json"))) return(NULL)
  ps <- read.csv(file.path(pd, "summary.csv"), colClasses = c(data_id = "character")); pp <- read.csv(file.path(pd, "points.csv"))
  pa <- read.csv(file.path(pd, "parameters.csv"))
  zc <- unname(as.matrix(pp[, if (ncol(x$z) == 1) "z" else c("z1", "z2"), drop = FALSE]))
  stopifnot(nrow(zc) == nrow(x$z), max(abs(zc - unname(x$z))) < 1e-9)
  pm <- function(v) if (v %in% pa$variable) pa$mean[pa$variable == v] else NA_real_
  list(source = "replicate_study", mean = pp$GP, err = ps$gp_relative_L2, noise = ps$relative_RMS_MC_noise,
       data_id = ps$data_id, rho_star = pm("rho_star"), ell1 = pm("ell1"), ell2 = pm("ell2"))
}

sens_report <- function(cfg) {
  files <- list.files(cfg$out, pattern = "^summary\\.csv$", recursive = TRUE, full.names = TRUE)
  files <- files[basename(dirname(dirname(files))) != "Tables"]
  if (!length(files)) { cat("No completed fits yet.\n"); return(invisible(NULL)) }
  baselines <- c("S0_base", "S3_q50")
  rows <- list(); missing <- character(0)
  for (f in files) {
    x <- readRDS(file.path(dirname(f), "result.rds")); s <- x$summary
    if (s$setting %in% baselines) next
    b <- sens_baseline(cfg, x)
    if (is.null(b)) { missing <- c(missing, dirname(f)); next }
    if (!identical(b$data_id, s$data_id)) stop("Data mismatch with the baseline for ", dirname(f))
    den <- sqrt(mean(x$truth^2))
    rows[[length(rows) + 1]] <- data.frame(block = s$block, setting = s$setting,
      label = cfg$science$settings[[s$setting]]$label, scenario = s$scenario, n = s$n, rep = s$rep,
      baseline = b$source, err_base = b$err, err = s$gp_relative_L2,
      rel_change = (s$gp_relative_L2 - b$err) / b$err,
      delta_mean = sqrt(mean((x$mean - b$mean)^2)) / den,
      # Monte Carlo envelope of delta_mean, on the same scale (truth norm): each fit's RMS MCSE is
      # its relative noise times the RMS of its own posterior mean. The two chains are treated as
      # independent (their covariance, e.g. from shared random numbers, is not available), so this
      # is a rough yardstick for delta_mean, not a test.
      mc_envelope = sqrt((s$relative_RMS_MC_noise * sqrt(mean(x$mean^2)))^2 +
                         (b$noise * sqrt(mean(b$mean^2)))^2) / den,
      rho_star_base = b$rho_star, rho_star = s$rho_star_mean,
      rho_star_q2.5 = s$rho_star_q2.5, rho_star_q97.5 = s$rho_star_q97.5,
      rho_star_frac_near_bound = s$rho_star_frac_near_bound,
      ell1_base = b$ell1, ell1 = s$ell1_mean, ell2_base = b$ell2, ell2 = s$ell2_mean,
      min_bulk_ess_hyper = s$min_bulk_ess_hyper, sampling_minutes = s$sampling_minutes,
      stringsAsFactors = FALSE)
  }
  if (length(missing)) cat(length(missing), "fit(s) without a completed baseline yet (skipped).\n")
  if (!length(rows)) { cat("No fit with a completed baseline yet.\n"); return(invisible(NULL)) }
  res <- do.call(rbind, rows)
  res <- res[order(res$scenario, match(res$block, SENS_BLOCKS), res$setting, res$n, res$rep), ]
  tab <- file.path(cfg$out, "Tables"); dir.create(tab, showWarnings = FALSE)
  write.csv(res, file.path(tab, "sensitivity_by_dataset.csv"), row.names = FALSE)
  key <- paste(res$scenario, res$setting, res$n, sep = "|")
  summ <- do.call(rbind, lapply(split(res, factor(key, levels = unique(key))), function(a) {
    v <- 100 * a$rel_change
    data.frame(scenario = a$scenario[1], block = a$block[1], setting = a$setting[1], label = a$label[1],
      n = a$n[1], datasets = nrow(a), reps = paste(a$rep, collapse = ","),
      rel_change_median_pct = median(v), rel_change_min_pct = min(v), rel_change_max_pct = max(v),
      cell = sprintf("%+.1f [%+.1f, %+.1f]", median(v), min(v), max(v)),
      err_base_median = median(a$err_base), err_median = median(a$err),
      delta_mean_median_pct = 100 * median(a$delta_mean), delta_mean_max_pct = 100 * max(a$delta_mean),
      mc_envelope_median_pct = 100 * median(a$mc_envelope),
      rho_star_base_median = median(a$rho_star_base), rho_star_median = median(a$rho_star),
      rho_star_frac_near_bound_max = max(a$rho_star_frac_near_bound),
      ell1_base_median = median(a$ell1_base), ell1_median = median(a$ell1),
      ell2_base_median = median(a$ell2_base), ell2_median = median(a$ell2),
      sampling_minutes_median = median(a$sampling_minutes), stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL
  write.csv(summ, file.path(tab, "sensitivity_summary.csv"), row.names = FALSE)
  print(summ[, c("scenario", "setting", "n", "datasets", "cell", "delta_mean_median_pct",
                 "mc_envelope_median_pct")], row.names = FALSE, digits = 3)
  cat("Tables written to", tab, "\n")
  invisible(summ)
}

# ---------------------------------------------------------------- controller
sens_args <- function(args) {
  known <- c("--blocks", "--settings", "--scenarios", "--n", "--reps", "--workers", "--out",
             "--smoke", "--report", "--list")
  bad <- args[!sub("=.*$", "", args) %in% known]
  if (length(bad)) stop("Unknown option(s): ", paste(bad, collapse = " "))
  value <- function(key, default) {
    a <- args[startsWith(args, paste0(key, "="))]
    if (length(a) > 1) stop("Duplicate option ", key)
    if (length(a)) sub(paste0("^", key, "="), "", a) else default
  }
  smoke <- "--smoke" %in% args
  list(smoke = smoke, report = "--report" %in% args, list = "--list" %in% args,
       blocks = value("--blocks", "all"), settings = value("--settings", NA_character_),
       scenarios = value("--scenarios", "1d_skn,2d_hole"),
       n = value("--n", if (smoke) "6" else "50,1000"), reps = value("--reps", if (smoke) "98" else "1:5"),
       workers = as.integer(value("--workers", "2")),
       out = value("--out", if (smoke) "Results/Sensitivity_Smoke" else "Results/Sensitivity"))
}

sens_ints <- function(x, what) {
  if (!grepl("^[0-9:,]+$", x)) stop(what, " must look like 1:5 or 1,3,7:9")
  sort(unique(as.integer(eval(parse(text = paste0("c(", x, ")"))))))
}

sens_main <- function(args) {
  o <- sens_args(args)
  settings <- sens_settings()
  if (o$list) {
    cat("Baseline: replicate-study fit of the same dataset (S3: S3_q50 on the m = 200 data).\n")
    for (id in names(settings)) cat(sprintf("  %-4s %-11s %s\n", settings[[id]]$block, id, settings[[id]]$label))
    return(invisible(NULL))
  }
  if (is.na(o$workers) || o$workers < 1L) stop("--workers must be a positive integer")
  dir.create(o$out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(o$out)
  science <- sens_science(o$smoke)
  cfg <- list(science = science, smoke = o$smoke, out = out,
              data_dir = if (o$smoke) file.path(out, "Smoke_Data") else DATA_DIR,
              repstudy_dir = REPSTUDY_DIR,
              guard = list(science = science, sources = as.list(tools::md5sum(SENS_SOURCES))))
  if (o$report) { sens_report(cfg); return(invisible(NULL)) }
  config_path <- file.path(out, "config.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path), cfg))
      stop("Sources or settings differ from those recorded in ", config_path,
           ". Use a new --out directory; existing results were preserved.")
  } else repstudy_atomic(cfg, config_path)

  # Selection: settings x scenarios x n x datasets.
  if (!is.na(o$settings)) {
    ids <- strsplit(o$settings, ",")[[1]]
    if (!all(ids %in% names(settings))) stop("Unknown setting(s); see --list")
  } else {
    blocks <- if (identical(o$blocks, "all")) SENS_BLOCKS else strsplit(o$blocks, ",")[[1]]
    if (!all(blocks %in% SENS_BLOCKS)) stop("Unknown block(s); use ", paste(SENS_BLOCKS, collapse = ","))
    ids <- names(settings)[vapply(settings, function(x) x$block %in% blocks, logical(1))]
  }
  if (any(grepl("^S3_", ids)) && !"S3_q50" %in% ids) {
    ids <- c(ids, "S3_q50"); cat("Adding S3_q50, the baseline of block S3.\n")
  }
  if (o$smoke) ids <- union("S0_base", ids) else ids <- setdiff(ids, "S0_base")
  scen <- strsplit(o$scenarios, ",")[[1]]
  if (!all(scen %in% names(experiment_scenarios()))) stop("Unknown scenario(s)")
  ns <- sens_ints(o$n, "--n"); reps <- sens_ints(o$reps, "--reps")
  if (!o$smoke && (!all(ns %in% c(50L, 250L, 500L, 1000L)) || any(reps < 1L | reps > 50L)))
    stop("n must be a replicate-study sample size and datasets must lie in 1..50")
  jobs <- expand.grid(setting = ids, scenario = scen, n = ns, rep = reps, stringsAsFactors = FALSE)
  jobs$block <- vapply(jobs$setting, function(id) settings[[id]]$block, "")
  jobs$dir <- sens_job_dir(out, jobs$block, jobs$setting, jobs$scenario, jobs$n, jobs$rep)
  # Longest fits first (2D, large n), for load balancing.
  jobs <- jobs[order(jobs$scenario != "2d_hole", -jobs$n, jobs$rep, jobs$setting), ]

  if (!o$smoke) {                          # baselines of S1, S2, S4-S6 must exist
    need <- unique(jobs[!jobs$block %in% c("S0", "S3"), c("scenario", "n", "rep")])
    pd <- repstudy_dir(REPSTUDY_DIR, need$scenario, need$n, need$rep)
    miss <- pd[!file.exists(file.path(pd, "final_status.json"))]
    if (length(miss)) stop("Replicate-study fits missing (baselines):\n", paste(miss, collapse = "\n"))
  }

  lock <- file.path(out, ".controller_lock")
  if (dir.exists(lock)) {
    pf <- file.path(lock, "pid")
    pid <- if (file.exists(pf)) suppressWarnings(as.integer(readLines(pf, warn = FALSE))) else NA_integer_
    if (length(pid) == 1 && !is.na(pid) && isTRUE(tryCatch(ps::ps_is_running(ps::ps_handle(as.integer(pid))), error = function(e) FALSE)))
      stop("Another controller is already running on ", out)
    unlink(lock, recursive = TRUE)
  }
  stopifnot(dir.create(lock)); writeLines(as.character(Sys.getpid()), file.path(lock, "pid"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  for (src in SENS_SOURCES) {
    dst <- file.path(out, "Source_Snapshot", basename(src))
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(dst)) stopifnot(file.copy(src, dst))
  }
  capture.output(sessionInfo(), file = file.path(out, paste0("sessionInfo_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".txt")))

  done <- file.exists(file.path(jobs$dir, "summary.csv"))
  # Generate the m = 200 datasets of S3 once, serially, before the workers start.
  s3 <- unique(jobs[!done & jobs$block == "S3", c("scenario", "n", "rep")])
  for (i in seq_len(nrow(s3))) {
    k <- sens_data_key(cfg, s3$scenario[i], s3$n[i], s3$rep[i], 200L)
    if (!file.exists(file.path(cfg$data_dir, paste0(k$key, ".rds")))) {
      cat("Simulating", k$key, "...\n"); flush.console()
      invisible(get_replicates(k$key, s3$n[i], k$sc$rho0, k$sc$ls, seed = k$seed, m = 200L,
                               data_dir = cfg$data_dir)); gc(FALSE)
    }
  }
  cat(sprintf("%d fits requested (%s), %d already complete, %d to run, %d worker(s). Output: %s\n",
              nrow(jobs), paste(ids, collapse = ","), sum(done), sum(!done), o$workers, out)); flush.console()

  queue <- which(!done); active <- list(); failures <- 0L
  log_path <- file.path(out, "job_log.csv")
  on.exit(for (a in active) if (a$proc$is_alive()) a$proc$kill_tree(), add = TRUE)
  t0 <- Sys.time()
  while (length(queue) || length(active)) {
    while (length(queue) && length(active) < o$workers) {
      j <- jobs[queue[1], ]; queue <- queue[-1]
      work <- file.path(j$dir, "Work"); dir.create(work, recursive = TRUE, showWarnings = FALSE)
      stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
      proc <- processx::process$new(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
        c("run_sensitivity.R", "--worker", config_path, j$setting, j$scenario,
          as.character(j$n), as.character(j$rep)),
        wd = getwd(), stdout = file.path(work, paste0("worker_", stamp, ".log")),
        stderr = file.path(work, paste0("worker_", stamp, "_stderr.log")),
        env = c("current", VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"),
        cleanup = TRUE)
      active[[length(active) + 1L]] <- list(proc = proc, job = j, started = Sys.time())
      cat(sprintf("[%s] start  %-11s %s n=%d r=%d\n", format(Sys.time(), "%H:%M:%S"), j$setting,
                  j$scenario, j$n, j$rep)); flush.console()
    }
    Sys.sleep(2)
    alive <- vapply(active, function(a) a$proc$is_alive(), logical(1))
    for (a in active[!alive]) {
      code <- a$proc$get_exit_status()
      ok <- identical(as.integer(code), 0L) && file.exists(file.path(a$job$dir, "summary.csv"))
      if (!ok) failures <- failures + 1L
      mins <- as.numeric(difftime(Sys.time(), a$started, units = "mins"))
      row <- data.frame(setting = a$job$setting, scenario = a$job$scenario, n = a$job$n, rep = a$job$rep,
                        status = if (ok) "complete" else "technical_failure",
                        exit_status = if (is.null(code)) NA_integer_ else as.integer(code),
                        started = format(a$started), finished = format(Sys.time()), wall_minutes = round(mins, 3))
      write.table(row, log_path, sep = ",", row.names = FALSE, qmethod = "double",
                  col.names = !file.exists(log_path), append = file.exists(log_path))
      cat(sprintf("[%s] %-6s %-11s %s n=%d r=%d (%.1f min)%s\n", format(Sys.time(), "%H:%M:%S"),
                  if (ok) "done" else "FAILED", a$job$setting, a$job$scenario, a$job$n, a$job$rep, mins,
                  if (ok) "" else paste0(" -- see ", file.path(a$job$dir, "Work"), "/worker_*_stderr.log")))
      flush.console()
    }
    active <- active[alive]
  }
  cat(sprintf("Batch finished in %.1f minutes; %d technical failure(s).\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), failures))
  if (failures) cat("Rerun the same command to retry failed fits; they resume from their checkpoints.\n")
  sens_report(cfg)
}

args <- if (exists("SENS_ARGS", envir = globalenv(), inherits = FALSE))
  get("SENS_ARGS", envir = globalenv()) else commandArgs(TRUE)
if (length(args) && args[1] == "--worker") {
  options(warn = 1)   # write each warning to the worker's stderr log as it occurs
  stopifnot(length(args) == 6L)
  cfg <- readRDS(args[2])
  if (!identical(cfg$guard$sources, as.list(tools::md5sum(SENS_SOURCES))))
    stop("Sources changed since this output directory was created; use a new --out directory.")
  sens_worker(cfg, args[3], args[4], as.integer(args[5]), as.integer(args[6]))
} else sens_main(args)
