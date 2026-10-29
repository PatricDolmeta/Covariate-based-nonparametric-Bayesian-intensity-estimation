###############################################################################
# run_replicate_study.R -- replicated estimation-error study, fixed-cap one-chain fits
# Protocol: SIMULATION_PROTOCOL.md, Section 4 (version 2.2, 12 October 2026).
#
# Run from 'Data Analysis/Simulation Study' (Terminal or RStudio's terminal):
#   Rscript run_replicate_study.R --smoke                          # tiny non-scientific check (~minutes)
#   Rscript run_replicate_study.R --cells=all --reps=1:5           # stage 1: also the cap validation
#   Rscript validate_replicate_study_cap.R                         # paired check against long fits
#   Rscript run_replicate_study.R --cells=all --reps=6:50          # stage 2
#   Rscript run_replicate_study.R --report                         # rebuild tables, no sampling
#
# Options
#   --cells=all | scenario:n[,scenario:n...]   e.g. --cells=1d_skn:250,2d_hole:1000
#   --reps=1:5 | 6:50 | 1,3,7:9                dataset (repetition) indices, 1..50
#   --workers=2                                datasets run in parallel, one R process each
#   --out=DIR                                  default Results/Replicate_Study
#                                              (smoke: Results/Replicate_Study_Smoke)
#
# Interrupting is safe: rerun the same command. Completed datasets are skipped after their
# verified record is re-checked; a running chain resumes from its last 10000-sweep checkpoint.
# The output directory is bound to the exact sources and settings recorded in config.rds; after
# any code change use a new --out directory.
#
# From the RStudio console (working directory 'Data Analysis/Simulation Study'):
#   REPSTUDY_ARGS <- c("--cells=all", "--reps=1:5", "--workers=4"); source("run_replicate_study.R")
###############################################################################
# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("Matrix", "coda", "sn", "jsonlite", "posterior", "processx", "ps")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })

REPSTUDY_SOURCES <- c("../R/data_generation.R", "../R/gp_sampler.R", "../R/competitors.R",
                  "../R/chain_diagnostics.R", "../R/replicate_study.R", "run_replicate_study.R")
for (f in setdiff(REPSTUDY_SOURCES, "run_replicate_study.R")) source(f)

repstudy_args <- function(args) {
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
       out = value("--out", if (smoke) "Results/Replicate_Study_Smoke" else "Results/Replicate_Study"),
       cells = value("--cells", if (smoke) "1d_skn:6,2d_hole:6" else NA_character_),
       reps = value("--reps", if (smoke) "98:99" else NA_character_),
       workers = as.integer(value("--workers", "2")))
}

repstudy_parse_cells <- function(x, smoke) {
  if (identical(x, "all"))
    return(data.frame(scenario = rep(c("1d_skn", "2d_hole"), each = 4),
                      n = rep(c(50L, 250L, 500L, 1000L), 2), stringsAsFactors = FALSE))
  parts <- strsplit(strsplit(x, ",")[[1]], ":")
  if (any(lengths(parts) != 2L)) stop("--cells must look like 1d_skn:250,2d_hole:1000")
  out <- data.frame(scenario = vapply(parts, `[`, "", 1),
                    n = suppressWarnings(as.integer(vapply(parts, `[`, "", 2))),
                    stringsAsFactors = FALSE)
  if (!all(out$scenario %in% names(experiment_scenarios())) || anyNA(out$n) || any(out$n < 2L))
    stop("Unknown scenario or invalid n in --cells")
  if (!smoke && !all(out$n %in% c(50L, 250L, 500L, 1000L)))
    stop("Replicate-study sample sizes are 50, 250, 500, 1000")
  unique(out)
}

repstudy_parse_reps <- function(x, smoke) {
  if (!grepl("^[0-9:,]+$", x)) stop("--reps must look like 1:5 or 1,3,7:9")
  r <- sort(unique(as.integer(eval(parse(text = paste0("c(", x, ")"))))))
  if (!length(r) || any(r < 1L) || (!smoke && any(r > 50L))) stop("Repetitions must lie in 1..50")
  r
}

repstudy_main <- function(args) {
  o <- repstudy_args(args)
  if (is.na(o$workers) || o$workers < 1L) stop("--workers must be a positive integer")
  dir.create(o$out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(o$out)
  science <- repstudy_science(o$smoke)
  cfg <- list(science = science, smoke = o$smoke, out = out,
              data_dir = if (o$smoke) file.path(out, "Smoke_Data") else DATA_DIR,
              guard = list(science = science, sources = as.list(tools::md5sum(REPSTUDY_SOURCES))))
  if (o$report) { repstudy_report(cfg); return(invisible(NULL)) }
  if (is.na(o$cells) || is.na(o$reps)) stop("Give --cells and --reps (or use --report or --smoke).")
  cells <- repstudy_parse_cells(o$cells, o$smoke); reps <- repstudy_parse_reps(o$reps, o$smoke)

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
  for (src in REPSTUDY_SOURCES) {
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
        c("run_replicate_study.R", "--worker", config_path, j$scenario, as.character(j$n), as.character(j$rep)),
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
}

a <- if (exists("REPSTUDY_ARGS", envir = globalenv(), inherits = FALSE))
  get("REPSTUDY_ARGS", envir = globalenv()) else commandArgs(TRUE)
if (length(a) && a[1] == "--worker") {
  options(warn = 1)   # write each warning to the worker's stderr log as it occurs
  stopifnot(length(a) == 5L)
  cfg <- readRDS(a[2])
  if (!identical(cfg$guard$sources, as.list(tools::md5sum(REPSTUDY_SOURCES))))
    stop("Sources changed since this output directory was created; use a new --out directory.")
  repstudy_worker(cfg, a[3], as.integer(a[4]), as.integer(a[5]))
} else if (Sys.getenv("REPSTUDY_SOURCE_ONLY") != "1") repstudy_main(a)
