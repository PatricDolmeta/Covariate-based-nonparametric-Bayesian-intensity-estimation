###############################################################################
# run_canada.R -- revised Bayesian analysis of the Canadian wildfire data (BA revision)
#
# Run from 'Data Analysis/Application', with Rscript:
#   Rscript run_canada.R --smoke                       # short workflow check (Results/Smoke)
#   Rscript run_canada.R                               # the 13 fits of the paper (Results/Canada)
#   Rscript run_canada.R --fits=ontario:full           # selected fits, comma separated
#   Rscript run_canada.R --report                      # tables only, no sampling
# Options (defaults): --chains=2 --workers=2 --sweeps=60000 --burn=20000 --thin=10
#   --extend=50000 --max_sweeps=260000 --V1=200 --m2=20 --m3=20 --G1=2000 --G2=200 --G3=50
#   --node_thin1=10 --node_thin2=10 --node_thin3=50 --rhat=1.01 --ess=200 --mc_noise=0.01
#   --a_ell=1 --b_ell=2.5 --lambda=0.1 (gamma part and tilt of the length-scale hyper-prior;
#   chosen on 11 October 2026 after the Ontario prior comparison, Prior Tuning/)
#   --out=Results/Canada
# Other priors, e.g. (as in the prior comparison):
#   Rscript run_canada.R --fits=ontario:temp,ontario:prec,ontario:wind --a_ell=3 --b_ell=3 --max_sweeps=60000 --out="Prior Tuning/Ontario_a3_b3"
#
# Each chain runs in its own R process (single-threaded BLAS), checkpointing every 10000
# sweeps. When all chains of a fit reach the current target, the fit is assessed
# (R/canada_posterior.R::assess_fit); if it fails the criteria, all its chains are extended by
# --extend sweeps (same burn-in, frozen tuning, exact continuation) up to --max_sweeps.
# Rerunning the same command resumes from the checkpoints. Changing any setting or source file
# requires a new --out directory.
###############################################################################

# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("Matrix", "coda", "posterior", "processx", "ps")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
source("../R/gp_sampler.R")
source("../R/canada_data.R")
source("../R/canada_posterior.R")
CANADA_SOURCES <- c("run_canada.R", "../R/gp_sampler.R", "../R/canada_data.R", "../R/canada_posterior.R")

DEFAULTS <- list(fits = "paper", chains = 2L, workers = 2L, sweeps = 60000L, burn = 20000L,
  thin = 10L, extend = 50000L, max_sweeps = 260000L, V1 = 200L, m2 = 20L, m3 = 20L,
  G1 = 2000L, G2 = 200L, G3 = 50L, node_thin1 = 10L, node_thin2 = 10L, node_thin3 = 50L,
  rhat = 1.01, ess = 200, mc_noise = 0.01, checkpoint_every = 10000L, seed_base = 20261010L,
  a_ell = 1, b_ell = 2.5, lambda = 0.1, out = file.path("Results", "Canada"))
# Settings that define the statistical/numerical configuration of an output directory
SETTING_KEYS <- c("chains", "sweeps", "burn", "thin", "extend", "max_sweeps", "V1", "m2", "m3",
  "G1", "G2", "G3", "node_thin1", "node_thin2", "node_thin3", "rhat", "ess", "mc_noise",
  "checkpoint_every", "seed_base")
# Prior: length-scale hyper-prior from the options (defaults = paper, as in the simulations), rho* fixed.
# The prior is recorded in config.rds and guarded like the settings.
prior_from <- function(o) list(a_ell = o$a_ell, b_ell = o$b_ell, lambda = o$lambda,
                               a_rho = 1, b_rho = 0.1, c_rho = 100)

parse_options <- function(args) {
  o <- DEFAULTS; o$smoke <- "--smoke" %in% args; o$report <- "--report" %in% args
  if (o$smoke) {
    o$fits <- "ontario:temp,ontario:temp_prec,ontario:full"; o$sweeps <- 300L; o$burn <- 100L
    o$extend <- 0L; o$max_sweeps <- 300L; o$checkpoint_every <- 100L
    o$out <- file.path("Results", "Smoke")
  }
  for (a in setdiff(args, c("--smoke", "--report"))) {
    if (!grepl("^--[A-Za-z0-9_]+=", a)) stop("Unrecognised argument: ", a)
    k <- sub("^--([A-Za-z0-9_]+)=.*$", "\\1", a); v <- sub("^--[A-Za-z0-9_]+=", "", a)
    if (!k %in% names(DEFAULTS)) stop("Unknown option: --", k)
    o[[k]] <- if (is.character(DEFAULTS[[k]])) v else if (is.integer(DEFAULTS[[k]])) as.integer(v) else as.numeric(v)
  }
  stopifnot(o$sweeps > o$burn, o$max_sweeps >= o$sweeps, o$chains >= 2, o$workers >= 1,
            o$a_ell > 0, o$b_ell > 0, o$lambda > 0)
  o
}

all_fit_ids <- function() as.vector(t(outer(names(CANADA_REGIONS), names(CANADA_MODELS), paste, sep = "__")))
resolve_fits <- function(spec) {
  if (identical(spec, "paper")) return(CANADA_PAPER_FITS$fit_id)
  ids <- sub(":", "__", strsplit(spec, ",")[[1]], fixed = TRUE)
  bad <- setdiff(ids, all_fit_ids())
  if (length(bad)) stop("Unknown fit(s): ", paste(bad, collapse = ", "),
                        ". Use region:model with region in {", paste(names(CANADA_REGIONS), collapse = ","),
                        "} and model in {", paste(names(CANADA_MODELS), collapse = ","), "}")
  ids
}
fit_seed <- function(cfg, fit_id, chain)
  as.integer(cfg$settings$seed_base + 1000L * match(fit_id, all_fit_ids()) + chain)

make_config <- function(o) {
  data_md5 <- vapply(names(CANADA_REGIONS), function(r) unname(tools::md5sum(canada_data_file(r))), "")
  list(settings = o[SETTING_KEYS], prior = prior_from(o),
       guard = list(sources = as.list(tools::md5sum(CANADA_SOURCES)), data = as.list(data_md5)),
       created = format(Sys.time()), R = R.version.string, BLAS = extSoftVersion()[["BLAS"]])
}

fit_setup <- function(cfg, fit_id) {
  parts <- strsplit(fit_id, "__", fixed = TRUE)[[1]]
  s <- cfg$settings; covs <- CANADA_MODELS[[parts[2]]]; d <- length(covs)
  R <- load_canada_region(parts[1])
  G <- c(s$G1, s$G2, s$G3)[d]
  binned <- bin_data(canada_reps(R, covs), d, G = G)
  nd <- if (d == 1) make_nodes(1, V = s$V1) else make_nodes(d, type = "grid", grid_m = c(NA, s$m2, s$m3)[d])
  list(region = parts[1], model = parts[2], covs = covs, d = d, R = R, G = G, binned = binned,
       nd = nd, z_eval = canada_probes(d), node_thin = c(s$node_thin1, s$node_thin2, s$node_thin3)[d])
}

# ------------------------------------------------------------------ worker: one chain
run_worker <- function(out, fit_id, chain, target) {
  cfg <- readRDS(file.path(out, "config.rds"))
  if (!identical(cfg$guard$sources, as.list(tools::md5sum(CANADA_SOURCES))))
    stop("Source files changed since ", out, " was created; use a new --out directory.")
  s <- cfg$settings; S <- fit_setup(cfg, fit_id)
  fdir <- file.path(out, "Fits", fit_folder(fit_id)); dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  seed <- fit_seed(cfg, fit_id, chain)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  cat(sprintf("%s chain %d: d=%d, nodes=%d, bins=%d/%d, events=%d, target=%d sweeps, seed=%d\n",
              fit_id, chain, S$d, nrow(S$nd$nodes), nrow(S$binned$z_dom), nrow(S$binned$z_obs),
              S$binned$N, target, seed))
  gp <- do.call(gp_pcn_mcmc, c(list(binned = S$binned, nd = S$nd, z_eval = S$z_eval,
    n_iter = target, burn = s$burn, thin = s$thin, save_nodes_every = S$node_thin,
    checkpoint_path = file.path(fdir, sprintf("chain%d_checkpoint.rds", chain)),
    checkpoint_every = s$checkpoint_every, checkpoint_guard = cfg$guard, verbose = TRUE), cfg$prior))
  res <- list(fit_id = fit_id, region = S$region, model = S$model, covs = S$covs, d = S$d,
              chain = chain, seed = seed, target = target, burn = s$burn, thin = s$thin,
              G = S$G, bins = c(dom = nrow(S$binned$z_dom), obs = nrow(S$binned$z_obs)),
              n = S$binned$n, N = S$binned$N, nd = S$nd, z_eval = S$z_eval,
              data = canada_data_summary(S$R), prior = cfg$prior, gp = gp,
              finished = format(Sys.time()), R = R.version.string, BLAS = extSoftVersion()[["BLAS"]])
  tmp <- paste0(chain_file(fdir, chain), ".pending")
  saveRDS(res, tmp)
  if (!file.rename(tmp, chain_file(fdir, chain))) stop("Could not save ", chain_file(fdir, chain))
  writeLines(as.character(target), chain_target_file(fdir, chain))
  cat("\nDone:", fit_id, "chain", chain, "at", target, "sweeps,", round(gp$time / 60, 1), "min\n")
}

# ------------------------------------------------------------------ master
status_file <- function(out, fit_id) file.path(out, "Fits", fit_folder(fit_id), "status.rds")
read_status <- function(out, fit_id, cfg) {
  f <- status_file(out, fit_id)
  if (file.exists(f)) readRDS(f) else list(target = cfg$settings$sweeps, final = FALSE)
}
write_status <- function(out, fit_id, st) {
  dir.create(dirname(status_file(out, fit_id)), recursive = TRUE, showWarnings = FALSE)
  saveRDS(st, status_file(out, fit_id))
}

evaluate_fit <- function(out, cfg, fit_id, target) {
  s <- cfg$settings; fdir <- file.path(out, "Fits", fit_folder(fit_id))
  summ <- assess_fit(fdir, seq_len(s$chains), list(rhat = s$rhat, ess = s$ess, mc_noise = s$mc_noise))
  cat(sprintf("[%s] assessed %s at %d sweeps: max R-hat %.4f, min bulk/tail ESS %.0f/%.0f, MC noise %.3f%% -> %s\n",
              format(Sys.time(), "%H:%M:%S"), fit_id, target, summ$max_rhat, summ$min_ess_bulk,
              summ$min_ess_tail, 100 * summ$mc_noise, if (summ$passed) "PASS" else "not passed"))
  next_target <- target + s$extend
  if (summ$passed || s$extend <= 0 || next_target > s$max_sweeps) {
    write_status(out, fit_id, list(target = target, final = TRUE, passed = summ$passed, summary = summ))
    if (!summ$passed) cat("  ", fit_id, ": criteria not met at the sweep limit; inspect", file.path(fdir, "traces.pdf"), "\n")
    return(NULL)
  }
  write_status(out, fit_id, list(target = next_target, final = FALSE, passed = FALSE, summary = summ))
  data.frame(fit_id = fit_id, chain = seq_len(s$chains), target = next_target, stringsAsFactors = FALSE)
}

write_report <- function(out) {
  ids <- all_fit_ids(); ids <- ids[file.exists(status_file(out, ids))]
  rows <- lapply(ids, function(id) {
    f <- status_file(out, id); if (!file.exists(f)) return(NULL)
    st <- readRDS(f); if (is.null(st$summary)) return(NULL)
    cbind(st$summary, final = st$final)
  })
  rows <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(rows)) { cat("No assessed fits in", out, "\n"); return(invisible(NULL)) }
  dir.create(file.path(out, "Tables"), showWarnings = FALSE)
  write.csv(rows, file.path(out, "Tables", "computation_summary.csv"), row.names = FALSE)
  print(rows[, c("fit_id", "sweeps", "median_minutes_per_chain", "acc_w", "max_rhat",
                 "min_ess_bulk", "min_ess_tail", "mc_noise", "passed", "final")], row.names = FALSE)
  cat("Written:", file.path(out, "Tables", "computation_summary.csv"), "\n")
  invisible(rows)
}

run_master <- function(o) {
  out <- o$out; dir.create(file.path(out, "Fits"), recursive = TRUE, showWarnings = FALSE)
  if (o$report) return(write_report(out))
  cfg <- make_config(o); cfg_path <- file.path(out, "config.rds")
  if (file.exists(cfg_path)) {
    old <- readRDS(cfg_path)
    if (!identical(old$settings, cfg$settings) || !identical(old$guard, cfg$guard) || !identical(old$prior, cfg$prior))
      stop("Settings, sources or data differ from those of ", out, "; use a new --out directory.")
    cfg <- old
  } else saveRDS(cfg, cfg_path)
  fits <- resolve_fits(o$fits); s <- cfg$settings
  queue <- NULL
  for (id in fits) {
    st <- read_status(out, id, cfg); if (isTRUE(st$final)) next
    fdir <- file.path(out, "Fits", fit_folder(id))
    need <- which(vapply(seq_len(s$chains), function(k) chain_target(fdir, k), 0L) < st$target)
    if (length(need)) queue <- rbind(queue, data.frame(fit_id = id, chain = need, target = st$target, stringsAsFactors = FALSE))
    else queue <- rbind(queue, evaluate_fit(out, cfg, id, st$target))
  }
  cat(sprintf("%d fit(s) requested; %d chain job(s) to run with %d worker(s). Output: %s\n",
              length(fits), if (is.null(queue)) 0L else nrow(queue), o$workers, out)); flush.console()
  active <- list(); failures <- 0L; t0 <- Sys.time()
  log_path <- file.path(out, "job_log.csv")
  on.exit(for (a in active) if (a$proc$is_alive()) a$proc$kill_tree(), add = TRUE)
  while ((!is.null(queue) && nrow(queue)) || length(active)) {
    while (!is.null(queue) && nrow(queue) && length(active) < o$workers) {
      j <- queue[1, ]; queue <- queue[-1, , drop = FALSE]
      wdir <- file.path(out, "Fits", fit_folder(j$fit_id), "Logs"); dir.create(wdir, recursive = TRUE, showWarnings = FALSE)
      stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
      proc <- processx::process$new(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
        c("run_canada.R", "--worker", out, j$fit_id, as.character(j$chain), as.character(j$target)),
        wd = getwd(), stdout = file.path(wdir, sprintf("chain%d_%s.log", j$chain, stamp)),
        stderr = file.path(wdir, sprintf("chain%d_%s_stderr.log", j$chain, stamp)),
        env = c("current", VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"),
        cleanup = TRUE)
      active[[length(active) + 1L]] <- list(proc = proc, job = j, started = Sys.time())
      cat(sprintf("[%s] start  %s chain %d -> %d sweeps\n", format(Sys.time(), "%H:%M:%S"),
                  j$fit_id, j$chain, j$target)); flush.console()
    }
    Sys.sleep(2)
    alive <- vapply(active, function(a) a$proc$is_alive(), logical(1))
    for (a in active[!alive]) {
      code <- a$proc$get_exit_status(); fdir <- file.path(out, "Fits", fit_folder(a$job$fit_id))
      ok <- identical(as.integer(code), 0L) && chain_target(fdir, a$job$chain) >= a$job$target
      if (!ok) failures <- failures + 1L
      mins <- as.numeric(difftime(Sys.time(), a$started, units = "mins"))
      row <- data.frame(fit_id = a$job$fit_id, chain = a$job$chain, target = a$job$target,
                        status = if (ok) "complete" else "failed", started = format(a$started),
                        finished = format(Sys.time()), wall_minutes = round(mins, 2))
      write.table(row, log_path, sep = ",", row.names = FALSE, col.names = !file.exists(log_path),
                  append = file.exists(log_path))
      cat(sprintf("[%s] %-6s %s chain %d (%.1f min)%s\n", format(Sys.time(), "%H:%M:%S"),
                  if (ok) "done" else "FAILED", a$job$fit_id, a$job$chain, mins,
                  if (ok) "" else paste0(" -- see ", file.path(fdir, "Logs"))))
      flush.console()
    }
    active <- active[alive]
    # assess fits whose chains have all reached their target and have no pending jobs
    for (id in fits) {
      st <- read_status(out, id, cfg); if (isTRUE(st$final)) next
      if (any(vapply(active, function(a) a$job$fit_id == id, logical(1)))) next
      if (!is.null(queue) && any(queue$fit_id == id)) next
      fdir <- file.path(out, "Fits", fit_folder(id))
      if (all(vapply(seq_len(s$chains), function(k) chain_target(fdir, k), 0L) >= st$target))
        queue <- rbind(queue, evaluate_fit(out, cfg, id, st$target))
    }
  }
  cat(sprintf("Finished in %.1f minutes; %d failed job(s).\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), failures))
  if (failures) cat("Rerun the same command to retry; chains resume from their checkpoints.\n")
  write_report(out)
}

args <- commandArgs(TRUE)
if (length(args) && args[1] == "--worker") {
  stopifnot(length(args) == 5L)
  run_worker(args[2], args[3], as.integer(args[4]), as.integer(args[5]))
} else if (Sys.getenv("CANADA_SOURCE_ONLY") != "1") run_master(parse_options(args))
