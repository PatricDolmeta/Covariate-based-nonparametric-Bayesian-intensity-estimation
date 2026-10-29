# run_long_fits.R
# Long two-chain fits on saved simulation datasets, for both scenarios (1d_skn, 2d_hole).
# Purpose: (i) posterior summaries of the hyper-parameters (paper Table Tab:Hyper) and MCMC
# diagnostics; (ii) reference fits for validating the one-chain replicate-study budget
# (validate_replicate_study_cap.R). Consolidated 12 October 2026 from run_long_fits_2d.R; same prior,
# sampler and numerical settings as the replicate study (R/replicate_study.R, repstudy_science()).
#
# Run from 'Data Analysis/Simulation Study':
#   Rscript run_long_fits.R --smoke                                  # workflow check on a tiny problem
#   Rscript run_long_fits.R --scenario=1d_skn  --n=50 --reps=1:5 [--workers=2] [--compact]
#   Rscript run_long_fits.R --scenario=2d_hole --n=1000 --reps=1:5 --workers=4 --compact
#   Rscript run_long_fits.R --report                                 # tables only, no sampling
#
# Settings: two chains per dataset from independent prior draws, burn-in 50000, thinning 10;
# data_seed = 1e6 d + 1000 n + r (saved datasets only, never regenerated); mcmc_seed =
# data_seed + 7 + 10000 * chain. Initial target 100000 sweeps in 1D, 200000 in 2D (250000 at
# n = 1000); a dataset failing the convergence rule is extended by 50000 sweeps, up to 400000.
# Convergence rule (agreed 9 October 2026): rank-normalised R-hat < 1.01 and bulk and tail
# ESS >= 200 for all monitored quantities (rho_star, ell_j, Poisson log-likelihood, mean
# intensity, rho(z) at fixed probe points). Chains checkpoint every 10000 sweeps; rerunning the
# same command resumes an interrupted batch.
# --compact: once a dataset passes and its summaries are written, delete its full chain and
# checkpoint files. Kept: diagnostics.csv, chains.csv, monitored_draws.rds and pooled_fit.rds
# (posterior mean and pointwise 95% bands). Without --compact nothing is deleted.
# Output: Results/Long_Fits/<scenario>/n<NNNN>/r<RR>/ and Results/Long_Fits/Tables/.
# From RStudio: LONG_ARGS <- c("--scenario=1d_skn", "--n=50", "--reps=1:5"); source("run_long_fits.R")

# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("Matrix", "coda", "sn", "jsonlite", "posterior", "processx", "ps")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
source("../R/data_generation.R"); source("../R/gp_sampler.R")
source("../R/chain_diagnostics.R"); source("../R/replicate_study.R")

SOURCES <- c("../R/data_generation.R", "../R/gp_sampler.R", "../R/chain_diagnostics.R",
             "../R/replicate_study.R")
LONG_ROOT <- file.path("Results", "Long_Fits")
SCENARIOS <- c("1d_skn", "2d_hole")

atomic_rds <- function(object, path) {
  tmp <- tempfile(pattern = "pending_", tmpdir = dirname(path)); on.exit(unlink(tmp))
  saveRDS(object, tmp); if (!file.rename(tmp, path)) stop("Could not save ", path)
}
root_dir <- function(smoke) if (smoke) file.path("Results", "Long_Fits_Smoke") else LONG_ROOT
settings_for <- function(scenario, n, smoke) {
  d <- experiment_scenarios()[[scenario]]$d; ps <- repstudy_science(FALSE)
  list(scenario = scenario, d = d, n = n,
    V_1d = if (smoke) 20L else ps$V_1d, grid_m = if (smoke) 7L else ps$grid_m,
    burn = if (smoke) 250L else 50000L, thin = if (smoke) 5L else 10L,
    initial = if (smoke) 1000L else if (d == 1) 100000L else if (n >= 1000L) 250000L else 200000L,
    increment = if (smoke) 500L else 50000L, max_iter = if (smoke) 1500L else 400000L,
    rhat_max = 1.01, ess_min = 200, prior = ps$prior)
}
parse_reps <- function(x) {
  if (grepl("^[0-9]+:[0-9]+$", x)) { b <- as.integer(strsplit(x, ":")[[1]]); seq(b[1], b[2]) }
  else as.integer(strsplit(x, ",")[[1]])
}

load_data <- function(cfg) {
  s <- cfg$settings; sc <- experiment_scenarios()[[s$scenario]]
  key <- sprintf("%s_n%d_r%d", s$scenario, s$n, cfg$rep)
  if (!cfg$smoke) stopifnot(file.exists(file.path(DATA_DIR, paste0(key, ".rds"))))
  get_replicates(key, s$n, sc$rho0, sc$ls, seed = cfg$data_seed, m = 50,
    data_dir = if (cfg$smoke) file.path(cfg$root, "Data") else DATA_DIR)
}

# ---------------------------------------------------------------- one chain (separate R process)
run_worker <- function(config_path, chain, target) {
  cfg <- readRDS(config_path); s <- cfg$settings; sc <- experiment_scenarios()[[s$scenario]]
  reps <- load_data(cfg); dataset <- attr(reps, "dataset")
  stopifnot(identical(dataset$data_id, cfg$data_id))
  mcmc_seed <- as.integer(cfg$data_seed + 7 + 10000 * chain)
  nd <- if (s$d == 1) make_nodes(1, V = s$V_1d) else make_nodes(2, type = "grid", grid_m = s$grid_m)
  z <- eval_grid(s$d)
  # Record the sampler's random cold start, then reset its stream (as in the existing long fits).
  set_mcmc_seed(mcmc_seed)
  ell <- r_ell_prior(s$d, s$prior$a_ell, s$prior$b_ell, s$prior$lambda)
  initial <- list(ell = ell, w = fmul(make_factor(nd, NULL, ell), rnorm(nrow(nd$nodes))),
    rho_star = rtgamma(s$prior$a_rho, s$prior$b_rho, s$prior$c_rho + log(s$n)))
  set_mcmc_seed(mcmc_seed)
  cat(s$scenario, "n:", s$n, "dataset:", cfg$rep, "chain:", chain, "target:", target, "data_id:", dataset$data_id,
      "mcmc_seed:", mcmc_seed, "\n")
  gp <- do.call(gp_pcn_mcmc, c(list(binned = bin_data(reps, s$d), nd = nd, z_eval = z,
    n_iter = target, burn = s$burn, thin = s$thin, ell_update = "noncentred", n_w_steps = 10,
    verbose = TRUE, checkpoint_path = file.path(cfg$dir, sprintf("chain_%02d_checkpoint.rds", chain)),
    checkpoint_every = if (cfg$smoke) 100L else 10000L, checkpoint_guard = cfg$guard), s$prior))
  stopifnot(all(is.finite(gp$draws)))
  atomic_rds(list(dataset = dataset, n = s$n, events = n_events(reps), chain = chain, mcmc_seed = mcmc_seed,
    initial = initial, settings = c(s, list(n_iter = target)), z = z, truth = sc$rho0(z), gp = gp),
    file.path(cfg$dir, sprintf("chain_%02d.rds", chain)))
  writeLines(as.character(target), file.path(cfg$dir, sprintf("chain_%02d.done", chain)))
  cat("Saved chain", chain, "at", target, "sweeps\n")
}

# ---------------------------------------------------------------- report (no sampling)
run_report <- function(smoke) {
  root <- root_dir(smoke); tab_dir <- file.path(root, "Tables")
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (f in list.files(root, pattern = "^final_status\\.csv$", recursive = TRUE, full.names = TRUE)) {
    st <- read.csv(f, colClasses = c(data_id = "character")); dir <- dirname(f)
    dg <- read.csv(file.path(dir, "diagnostics.csv"))
    ch <- read.csv(file.path(dir, "chains.csv"), colClasses = c(data_id = "character"))
    repfit <- sprintf("Results/Replicate_Study/%s/n%04d/r%02d/summary.csv", st$scenario, st$n, st$rep)
    if (!smoke && file.exists(repfit))
      stopifnot(identical(read.csv(repfit, colClasses = c(data_id = "character"))$data_id, ch$data_id[1]))
    g <- function(v) if (v %in% dg$variable) dg$mean[dg$variable == v] else NA_real_
    ell <- grepl("^ell[0-9]+$", dg$variable)
    rows[[length(rows) + 1]] <- data.frame(scenario = st$scenario, n = st$n, rep = st$rep,
      data_id = ch$data_id[1], chains = nrow(ch), sweeps = st$n_iter,
      ell1 = g("ell1"), ell2 = g("ell2"), rho_star = g("rho_star"),
      ell_min_ess_bulk = min(dg$ess_bulk[ell]), ell_max_rhat = max(dg$rhat[ell]),
      all_max_rhat = max(dg$rhat), all_min_ess_bulk = min(dg$ess_bulk), all_min_ess_tail = min(dg$ess_tail),
      passes_rule = all(is.finite(dg$rhat) & dg$rhat < 1.01 & dg$ess_bulk >= 200 & dg$ess_tail >= 200),
      chain_minutes_max = max(ch$seconds) / 60)
  }
  if (!length(rows)) { cat("No completed long fits yet.\n"); return(invisible(NULL)) }
  res <- do.call(rbind, rows); res <- res[order(res$scenario, res$n, res$rep), ]
  write.csv(res, file.path(tab_dir, "long_fit_results.csv"), row.names = FALSE)
  tab <- do.call(rbind, lapply(split(res, list(res$scenario, res$n), drop = TRUE), function(a)
    data.frame(scenario = a$scenario[1], n = a$n[1], datasets = nrow(a),
      reps = paste(a$rep, collapse = ","), ell1_median = median(a$ell1), ell2_median = median(a$ell2),
      rho_star_median = median(a$rho_star), ell_min_ess_bulk = min(a$ell_min_ess_bulk),
      ell_max_rhat = max(a$ell_max_rhat), sweeps_max = max(a$sweeps), all_pass_rule = all(a$passes_rule))))
  tab <- tab[order(tab$scenario, tab$n), ]
  write.csv(tab, file.path(tab_dir, "long_fit_table.csv"), row.names = FALSE)
  print(tab, row.names = FALSE)
  invisible(tab)
}

# ---------------------------------------------------------------- controller
run_main <- function(args) {
  value <- function(key, default) {
    a <- args[startsWith(args, paste0(key, "="))]
    if (length(a) > 1) stop("Duplicate option ", key)
    if (length(a)) sub(paste0("^", key, "="), "", a) else default
  }
  smoke <- "--smoke" %in% args; compact <- "--compact" %in% args
  scenario <- if (smoke) value("--scenario", "2d_hole") else value("--scenario", NA)
  if (!scenario %in% SCENARIOS) stop("Give --scenario=1d_skn or --scenario=2d_hole")
  dim_z <- experiment_scenarios()[[scenario]]$d
  n <- if (smoke) 6L else as.integer(value("--n", NA))
  reps <- if (smoke) 98L else parse_reps(value("--reps", ""))
  workers <- as.integer(value("--workers", "2"))
  stopifnot(n %in% c(6L, 50L, 250L, 500L, 1000L), length(reps) >= 1, all(is.finite(reps)),
            all(reps >= 1), is.finite(workers), workers >= 1)
  root <- file.path(root_dir(smoke), scenario); s <- settings_for(scenario, n, smoke)
  guard <- list(settings = s, sources = as.list(tools::md5sum(SOURCES)))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = file.path(root, sprintf("sessionInfo_n%04d.txt", n)))
  ds <- lapply(reps, function(r) {
    dir <- file.path(root, sprintf("n%04d", n), sprintf("r%02d", r))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cfg <- list(settings = s, rep = r, data_seed = repstudy_data_seed(dim_z, n, r), smoke = smoke,
                root = root, dir = normalizePath(dir), guard = guard)
    x <- load_data(cfg); cfg$data_id <- attr(x, "dataset")$data_id; rm(x)
    path <- file.path(dir, "config.rds")
    if (file.exists(path)) {
      if (!identical(readRDS(path), cfg)) stop("Settings or sources changed for ", dir, "; existing output preserved.")
    } else atomic_rds(cfg, path)
    list(rep = r, dir = cfg$dir, config = normalizePath(path), data_id = cfg$data_id,
         stage = NA_integer_, done = file.exists(file.path(dir, "final_status.csv")))
  })
  queue <- list(); active <- list(); processes <- list()
  on.exit(for (p in processes) if (p$is_alive()) p$kill_tree(), add = TRUE)
  chain_done <- function(dir, ch, stage) {
    f <- file.path(dir, sprintf("chain_%02d.done", ch))
    file.exists(f) && as.integer(readLines(f, warn = FALSE)) == stage
  }
  enqueue_stage <- function(k) {
    writeLines(as.character(ds[[k]]$stage), file.path(ds[[k]]$dir, "stage.txt"))
    for (ch in 1:2) if (!chain_done(ds[[k]]$dir, ch, ds[[k]]$stage)) queue[[length(queue) + 1]] <<- c(k, ch)
  }
  waiting <- function(k) any(vapply(queue, function(j) j[1] == k, logical(1))) ||
    any(vapply(active, function(a) a$k == k, logical(1)))
  check_dataset <- function(k) {
    d <- ds[[k]]; paths <- file.path(d$dir, sprintf("chain_%02d.rds", 1:2))
    summarise_pilot(paths, d$dir, NULL)
    dg <- read.csv(file.path(d$dir, "diagnostics.csv"))
    a <- list(max_rhat = max(dg$rhat), min_ess_bulk = min(dg$ess_bulk), min_ess_tail = min(dg$ess_tail),
      pass = isTRUE(all(is.finite(dg$rhat) & is.finite(dg$ess_bulk) & is.finite(dg$ess_tail) &
        dg$rhat < s$rhat_max & dg$ess_bulk >= s$ess_min & dg$ess_tail >= s$ess_min)))
    cat(sprintf("%s n = %d, dataset %d, %d sweeps: max R-hat %.4f, min bulk ESS %.0f, min tail ESS %.0f -> %s\n",
      scenario, n, d$rep, d$stage, a$max_rhat, a$min_ess_bulk, a$min_ess_tail, if (a$pass) "pass" else "fail"))
    flush.console()
    if (a$pass || (smoke && d$stage >= s$max_iter)) {
      write.csv(data.frame(scenario = scenario, n = n, rep = d$rep, data_id = d$data_id, n_iter = d$stage, a, smoke = smoke,
        finished = as.character(Sys.time())), file.path(d$dir, "final_status.csv"), row.names = FALSE)
      ds[[k]]$done <<- TRUE
      if (compact) unlink(c(paths, file.path(d$dir, sprintf("chain_%02d_checkpoint.rds", 1:2))))
    } else if (d$stage + s$increment <= s$max_iter) {
      ds[[k]]$stage <<- d$stage + s$increment; enqueue_stage(k)
    } else {
      write.csv(data.frame(scenario = scenario, n = n, rep = d$rep, n_iter = d$stage, a), file.path(d$dir, "needs_review.csv"),
        row.names = FALSE)
      ds[[k]]$done <<- TRUE
      warning("Dataset ", d$rep, " fails the rule at ", d$stage, " sweeps; files kept for review.")
    }
  }
  for (k in seq_along(ds)) if (!ds[[k]]$done) {
    sf <- file.path(ds[[k]]$dir, "stage.txt")
    ds[[k]]$stage <- if (file.exists(sf)) as.integer(readLines(sf, warn = FALSE)) else s$initial
    enqueue_stage(k)
  }
  for (k in seq_along(ds)) if (!ds[[k]]$done && !waiting(k)) check_dataset(k)
  last_report <- Sys.time()
  while (length(queue) || length(active)) {
    while (length(queue) && length(active) < workers) {
      j <- queue[[1]]; queue <- queue[-1]; d <- ds[[j[1]]]
      log <- file.path(d$dir, sprintf("chain_%02d_iter%d.log", j[2], d$stage))
      p <- processx::process$new(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
        c("run_long_fits.R", "--worker", d$config, as.character(j[2]), as.character(d$stage)),
        wd = getwd(), stdout = log, stderr = sub("\\.log$", "_stderr.log", log), cleanup = TRUE,
        env = c("current", VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"))
      active[[length(active) + 1]] <- list(k = j[1], chain = j[2], p = p, log = log)
      processes[[length(processes) + 1]] <- p
      cat(format(Sys.time()), "started n =", n, "dataset", d$rep, "chain", j[2], "target", d$stage, "\n")
      flush.console()
    }
    Sys.sleep(2)
    alive <- vapply(active, function(a) a$p$is_alive(), logical(1))
    for (i in rev(which(!alive))) {
      a <- active[[i]]; status <- a$p$get_exit_status(); active <- active[-i]
      if (is.null(status) || status != 0)
        stop("Chain failed: n = ", n, ", dataset ", ds[[a$k]]$rep, ", chain ", a$chain, "; see ", a$log)
      if (!waiting(a$k)) check_dataset(a$k)
    }
    if (as.numeric(difftime(Sys.time(), last_report, units = "secs")) >= 60 && length(active)) {
      for (a in active) if (file.exists(a$log)) {
        parts <- unlist(strsplit(paste(readLines(a$log, warn = FALSE), collapse = "\n"), "[\r\n]"))
        prog <- parts[grepl("^it [0-9]+", parts)]
        cat(" dataset", ds[[a$k]]$rep, "chain", a$chain, ":", if (length(prog)) tail(prog, 1) else "initialising", "\n")
      }
      flush.console(); last_report <- Sys.time()
    }
  }
  cat("Batch finished: n =", n, "datasets", paste(reps, collapse = ","), "\n")
  run_report(smoke)
}

# From RStudio: set the working directory to this folder, then e.g.
#   LONG_ARGS <- c("--smoke"); source("run_long_fits.R")
#   LONG_ARGS <- c("--scenario=1d_skn", "--n=50", "--reps=1:5", "--workers=2"); source("run_long_fits.R")
# (the chains still run as separate Rscript processes; the console stays busy until the batch ends).
args <- if (exists("LONG_ARGS", envir = globalenv(), inherits = FALSE))
  get("LONG_ARGS", envir = globalenv()) else commandArgs(TRUE)
if (length(args) && args[1] == "--worker") {
  options(warn = 1)   # write each warning to the worker's stderr log as it occurs
  stopifnot(length(args) == 4)
  run_worker(args[2], as.integer(args[3]), as.integer(args[4]))
} else if ("--report" %in% args) {
  run_report("--smoke" %in% args)
} else run_main(args)
