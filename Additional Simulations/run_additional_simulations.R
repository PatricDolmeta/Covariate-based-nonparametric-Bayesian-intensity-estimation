###############################################################################
# run_additional_simulations.R -- additional simulation experiments, Supplement S4.4-S4.6
#
# New data-generating settings (additional_scenarios.R) run through the unchanged engine of the
# main replicate study (../R/replicate_study.R): for each dataset, the kernel estimator with oracle
# and CV bandwidths, the log-polynomials of degree 1-4 and one GP chain (60000 sweeps, 20000
# burn-in, thinning 10), with the prior and all numerical settings of the main study.
#
# Run from 'Data Analysis/Additional Simulations' (Terminal or RStudio, see README):
#   Rscript run_additional_simulations.R --smoke                        # tiny non-scientific check
#   Rscript run_additional_simulations.R --cells=all --reps=1:5         # stage 1
#   Rscript run_additional_simulations.R --cells=all --reps=6:50        # stage 2
#   Rscript run_additional_simulations.R --report                       # rebuild tables, no sampling
#
# Options
#   --cells=all | scenario[:n][,scenario[:n]...]
#        e.g. --cells=1d_negexp,1d_minmax (all four n) or --cells=2d_slope:250,2d_detcov:1000
#        scenarios: 1d_negexp 1d_minmax 2d_slope 2d_twobump 2d_detcov 2d_overparam
#        default sample sizes n = 50, 250, 500, 1000 (any n in 2..1000 may be given explicitly)
#   --reps=1:5 | 6:50 | 1,3,7:9                dataset indices, 1..100
#   --workers=2                                datasets run in parallel, one R process each
#   --out=DIR                                  default Results/Additional_Simulations
#                                              (smoke: Results/Additional_Simulations_Smoke)
#
# Interrupting is safe: rerun the same command. Completed datasets are skipped after their
# verified record is re-checked; a running chain resumes from its last checkpoint. The output
# directory is bound to the exact sources and settings recorded in config.rds (including the
# files in ../R/); after any code change use a new --out directory.
#
# From the RStudio console (working directory 'Data Analysis/Additional Simulations'):
#   ADDSIM_ARGS <- c("--cells=all", "--reps=1:5", "--workers=3"); source("run_additional_simulations.R")
###############################################################################
local({ need <- c("Matrix", "coda", "sn", "jsonlite", "posterior", "processx", "ps")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })

ADDSIM_SOURCES <- c("../R/data_generation.R", "../R/gp_sampler.R", "../R/competitors.R",
                    "../R/chain_diagnostics.R", "../R/replicate_study.R",
                    "additional_scenarios.R", "run_additional_simulations.R")
if (!all(file.exists(ADDSIM_SOURCES)))
  stop("Run this script from 'Data Analysis/Additional Simulations' (missing: ",
       paste(ADDSIM_SOURCES[!file.exists(ADDSIM_SOURCES)], collapse = ", "), ")")
for (f in setdiff(ADDSIM_SOURCES, "run_additional_simulations.R")) source(f)

ADDSIM_MAIN_DATA <- "../Simulation Study/Data"                                # 2d_overparam input
ADDSIM_MAIN_ERRORS <- "../Simulation Study/Results/Replicate_Study/Tables/method_errors.csv"

addsim_args <- function(args) {
  known <- c("--cells", "--reps", "--workers", "--out", "--smoke", "--report")
  bad <- args[!sub("=.*$", "", args) %in% known]
  if (length(bad)) stop("Unknown option(s): ", paste(bad, collapse = " "))
  value <- function(key, default) {
    a <- args[startsWith(args, paste0(key, "="))]
    if (length(a) > 1) stop("Duplicate option ", key)
    if (length(a)) sub(paste0("^", key, "="), "", a) else default
  }
  smoke <- "--smoke" %in% args
  list(smoke = smoke, report = "--report" %in% args,
       out = value("--out", if (smoke) "Results/Additional_Simulations_Smoke" else "Results/Additional_Simulations"),
       cells = value("--cells", if (smoke) paste0(names(ADDSIM_SCENARIOS), ":6", collapse = ",") else NA_character_),
       reps = value("--reps", if (smoke) "99" else NA_character_),
       workers = as.integer(value("--workers", "2")))
}

addsim_parse_cells <- function(x, smoke) {
  if (identical(x, "all")) x <- paste(names(ADDSIM_SCENARIOS), collapse = ",")
  parts <- strsplit(strsplit(x, ",")[[1]], ":")
  if (any(lengths(parts) > 2L)) stop("--cells must look like 1d_negexp or 2d_slope:250")
  out <- do.call(rbind, lapply(parts, function(p) {
    n <- if (length(p) == 2L) suppressWarnings(as.integer(p[2])) else ADDSIM_N
    data.frame(scenario = p[1], n = n, stringsAsFactors = FALSE)
  }))
  if (!all(out$scenario %in% names(ADDSIM_SCENARIOS)))
    stop("Unknown scenario(s): ", paste(setdiff(out$scenario, names(ADDSIM_SCENARIOS)), collapse = ", "),
         ". Available: ", paste(names(ADDSIM_SCENARIOS), collapse = ", "))
  if (anyNA(out$n) || any(out$n < 2L) || any(out$n > 1000L)) stop("Sample sizes must lie in 2..1000")
  unique(out)
}

addsim_parse_reps <- function(x) {
  if (!grepl("^[0-9:,]+$", x)) stop("--reps must look like 1:5 or 1,3,7:9")
  r <- sort(unique(as.integer(eval(parse(text = paste0("c(", x, ")"))))))
  if (!length(r) || any(r < 1L) || any(r > 100L)) stop("Repetitions must lie in 1..100")
  r
}

# ---------------------------------------------------------------- extra tables (no sampling)
# Supplement-style tables: truth norms, mean (SD) of the relative L2 errors per scenario and n,
# and, for 2d_overparam, the paired comparison with the main-study 1d_skn fits of the same data.
addsim_report_extras <- function(cfg) {
  tab <- file.path(cfg$out, "Tables"); ef <- file.path(tab, "method_errors.csv")
  if (!file.exists(ef)) return(invisible(NULL))
  err <- read.csv(ef, colClasses = c(data_id = "character"))
  scen <- intersect(names(ADDSIM_SCENARIOS), unique(err$scenario))
  norms <- do.call(rbind, lapply(scen, function(s) {
    sc <- ADDSIM_SCENARIOS[[s]]; z <- eval_grid(sc$d); v <- sc$rho0(z)
    data.frame(scenario = s, section = sc$section, equation = sc$equation, formula = sc$formula,
               L2_norm = sqrt(mean(v^2)), minimum = min(v), maximum = max(v))
  }))
  write.csv(norms, file.path(tab, "truth_norms.csv"), row.names = FALSE)

  lines <- c("# Additional simulation experiments (Supplement S4.4-S4.6)", "",
    sprintf("Updated %s. Engine and settings of the main replicate study (%d sweeps, %d burn-in, thinning %d).%s",
            format(Sys.time(), "%Y-%m-%d %H:%M"), cfg$science$n_iter, cfg$science$burn, cfg$science$thin,
            if (cfg$smoke) " SMOKE TEST: tiny settings, not scientific output." else ""),
    "Cells: mean relative L2 error over datasets (between-dataset SD); L2 norms of rho0 on [0,1]^d (reporting grid).", "")
  for (s in scen) {
    e <- err[err$scenario == s, ]; ns <- sort(unique(e$n)); nr <- norms[norms$scenario == s, ]
    lines <- c(lines, sprintf("## %s (%s, Eq. %s): ||rho0||_L2 = %.2f", s, nr$section, nr$equation, nr$L2_norm), "",
      paste0("| Method | ", paste(sprintf("n = %d", ns), collapse = " | "), " |"),
      paste0("|---|", strrep("---:|", length(ns))),
      paste0("| Datasets | ", paste(vapply(ns, function(n) as.character(sum(e$n == n & e$method == "GP")), ""), collapse = " | "), " |"))
    for (m in REPSTUDY_METHODS) {
      cell <- vapply(ns, function(n) {
        v <- e$relative_L2[e$n == n & e$method == m]
        if (!length(v)) "-" else if (length(v) == 1) sprintf("%.3f", v) else sprintf("%.3f (%.3f)", mean(v), sd(v))
      }, "")
      lines <- c(lines, sprintf("| %s | %s |", REPSTUDY_LABELS[[m]], paste(cell, collapse = " | ")))
    }
    lines <- c(lines, "")
  }

  if ("2d_overparam" %in% scen && !cfg$smoke && file.exists(ADDSIM_MAIN_ERRORS)) {
    main <- read.csv(ADDSIM_MAIN_ERRORS, colClasses = c(data_id = "character"))
    main <- main[main$scenario == "1d_skn", c("n", "rep", "data_id", "method", "relative_L2")]
    names(main)[5] <- "main_relative_L2"
    op <- merge(err[err$scenario == "2d_overparam", c("n", "rep", "data_id", "method", "relative_L2")],
                main, by = c("n", "rep", "data_id", "method"))
    if (nrow(op)) {
      op$difference <- op$relative_L2 - op$main_relative_L2
      cmp <- do.call(rbind, lapply(split(op, list(op$n, op$method), drop = TRUE), function(v) data.frame(
        n = v$n[1], method = v$method[1], datasets = nrow(v),
        mean_overparam = mean(v$relative_L2), mean_main_1d = mean(v$main_relative_L2),
        mean_difference = mean(v$difference),
        MCSE_difference = if (nrow(v) > 1) sd(v$difference) / sqrt(nrow(v)) else NA_real_,
        fraction_overparam_larger = mean(v$difference > 0))))
      cmp <- cmp[order(match(cmp$method, REPSTUDY_METHODS), cmp$n), ]; rownames(cmp) <- NULL
      write.csv(cmp, file.path(tab, "overparam_vs_main_1d.csv"), row.names = FALSE)
      g <- cmp[cmp$method == "GP", ]
      lines <- c(lines, "## 2d_overparam vs main-study 1d_skn fits of the same datasets (GP)", "",
        "| n | Datasets | Over-parametrised | Main study (1D) | Difference (MCSE) | Over-param. larger |",
        "|---:|---:|---:|---:|---:|---:|",
        sprintf("| %d | %d | %.4f | %.4f | %+.4f (%s) | %.0f%% |", g$n, g$datasets, g$mean_overparam,
                g$mean_main_1d, g$mean_difference,
                ifelse(is.na(g$MCSE_difference), "-", sprintf("%.4f", g$MCSE_difference)),
                100 * g$fraction_overparam_larger),
        "", "All methods: Tables/overparam_vs_main_1d.csv (paired by data_id).", "")
    }
  }
  lines <- c(lines, "Tables: truth_norms.csv, error_table.csv, paired_differences_vs_GP.csv, computation_table.csv, method_errors.csv, dataset_summary.csv.")
  writeLines(lines, file.path(cfg$out, "ADDITIONAL_SUMMARY.md"))
  cat(lines, sep = "\n"); cat("\n")
  invisible(norms)
}

# ---------------------------------------------------------------- controller
addsim_main <- function(args) {
  o <- addsim_args(args)
  if (is.na(o$workers) || o$workers < 1L) stop("--workers must be a positive integer")
  dir.create(o$out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(o$out)
  science <- repstudy_science(o$smoke)
  cfg <- list(science = science, smoke = o$smoke, out = out,
              data_dir = if (o$smoke) file.path(out, "Smoke_Data") else DATA_DIR,
              main_data_dir = if (o$smoke) NULL else ADDSIM_MAIN_DATA,
              guard = list(science = science, sources = as.list(tools::md5sum(ADDSIM_SOURCES))))
  if (o$report) { repstudy_report(cfg); addsim_report_extras(cfg); return(invisible(NULL)) }
  if (is.na(o$cells) || is.na(o$reps)) stop("Give --cells and --reps (or use --report or --smoke).")
  cells <- addsim_parse_cells(o$cells, o$smoke); reps <- addsim_parse_reps(o$reps)

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

  config_path <- file.path(out, "config.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path), cfg))
      stop("Sources or settings differ from those recorded in ", config_path,
           ". Use a new --out directory; existing results were preserved.")
  } else repstudy_atomic(cfg, config_path)
  for (src in ADDSIM_SOURCES) {
    dst <- file.path(out, "Source_Snapshot", basename(src))
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(dst)) stopifnot(file.copy(src, dst))
    stopifnot(identical(unname(tools::md5sum(dst)), unname(tools::md5sum(src))))
  }
  capture.output(sessionInfo(), file = file.path(out, paste0("sessionInfo_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".txt")))

  # Repetition-major order: every requested cell gets dataset r before any cell gets r+1.
  jobs <- do.call(rbind, lapply(reps, function(r) data.frame(cells, rep = r)))
  jobs$dir <- repstudy_dir(out, jobs$scenario, jobs$n, jobs$rep)
  done <- file.exists(file.path(jobs$dir, "final_status.json"))
  cat(sprintf("%d datasets requested, %d already complete, %d to run, %d worker(s). Output: %s\n",
              nrow(jobs), sum(done), sum(!done), o$workers, out)); flush.console()

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
        c("run_additional_simulations.R", "--worker", config_path, j$scenario, as.character(j$n), as.character(j$rep)),
        wd = getwd(), stdout = file.path(work, paste0("worker_", stamp, ".log")),
        stderr = file.path(work, paste0("worker_", stamp, "_stderr.log")),
        env = c("current", VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                OMP_NUM_THREADS = "1"),
        cleanup = TRUE)
      active[[length(active) + 1L]] <- list(proc = proc, job = j, started = Sys.time())
      cat(sprintf("[%s] start  %s n=%d r=%d\n", format(Sys.time(), "%H:%M:%S"),
                  j$scenario, j$n, j$rep)); flush.console()
    }
    Sys.sleep(2)
    alive <- vapply(active, function(a) a$proc$is_alive(), logical(1))
    for (a in active[!alive]) {
      code <- a$proc$get_exit_status()
      ok <- identical(as.integer(code), 0L) && file.exists(file.path(a$job$dir, "final_status.json"))
      if (!ok) failures <- failures + 1L
      mins <- as.numeric(difftime(Sys.time(), a$started, units = "mins"))
      row <- data.frame(scenario = a$job$scenario, n = a$job$n, rep = a$job$rep,
                        status = if (ok) "complete" else "technical_failure",
                        exit_status = if (is.null(code)) NA_integer_ else as.integer(code),
                        started = format(a$started), finished = format(Sys.time()),
                        wall_minutes = round(mins, 3))
      write.table(row, log_path, sep = ",", row.names = FALSE, qmethod = "double",
                  col.names = !file.exists(log_path), append = file.exists(log_path))
      cat(sprintf("[%s] %-6s %s n=%d r=%d (%.1f min)%s\n", format(Sys.time(), "%H:%M:%S"),
                  if (ok) "done" else "FAILED", a$job$scenario, a$job$n, a$job$rep, mins,
                  if (ok) "" else paste0(" -- see ", file.path(a$job$dir, "Work"), "/worker_*_stderr.log")))
      flush.console()
    }
    active <- active[alive]
  }
  cat(sprintf("Batch finished in %.1f minutes; %d technical failure(s).\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), failures))
  if (failures) cat("Rerun the same command to retry failed datasets; they resume from their checkpoints.\n")
  repstudy_report(cfg)
  addsim_report_extras(cfg)
}

a <- if (exists("ADDSIM_ARGS", envir = globalenv(), inherits = FALSE))
  get("ADDSIM_ARGS", envir = globalenv()) else commandArgs(TRUE)
if (length(a) && a[1] == "--worker") {
  options(warn = 1)   # write each warning to the worker's stderr log as it occurs
  stopifnot(length(a) == 5L)
  cfg <- readRDS(a[2])
  if (!identical(cfg$guard$sources, as.list(tools::md5sum(ADDSIM_SOURCES))))
    stop("Sources changed since this output directory was created; use a new --out directory.")
  repstudy_worker(cfg, a[3], as.integer(a[4]), as.integer(a[5]))
} else if (Sys.getenv("ADDSIM_SOURCE_ONLY") != "1") addsim_main(a)
