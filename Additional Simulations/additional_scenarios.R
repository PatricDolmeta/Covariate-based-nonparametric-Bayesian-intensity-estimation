###############################################################################
# additional_scenarios.R -- data-generating settings of the additional simulation experiments
# (Supplement Sections S4.4-S4.6), plugged into the replicate-study engine of ../R/.
#
# Sourced AFTER ../R/data_generation.R, gp_sampler.R, competitors.R, chain_diagnostics.R and
# replicate_study.R. It defines the new ground truths and data loaders, and redefines
# repstudy_load_data() and experiment_scenarios() so that the unchanged engine of the main
# study (competitors -> one GP chain -> single-chain summaries -> verified compact export) runs
# on the new datasets. Nothing in ../R/ is modified: prior, sampler settings (60000 sweeps,
# 20000 burn-in, thinning 10), GP discretisation (V = 200 in 1D, 20 x 20 grid in 2D), binning,
# kernel bandwidth grids and log-polynomial competitors are exactly those of the main study.
#
# Scenarios. Window W = [0,1]^2 (the paper's [-1/2,1/2]^2 shifted), m = 50 pixels per axis; the
# random covariates are the GP fields of the main study (length-scales 0.005 and 0.05).
#   1d_negexp     S4.4, Eq. (1DTruth2): rho0(z) = 2 exp(3 (1 - z) - 1)
#   1d_minmax     S4.4, Eq. (1DTruth3): rho0(z) = 2 + 2 P(z; 3/4, 3/8) - 2 P(z; 1/4, 3/8),
#                 P(z; c, a) = 1 - p(2 delta(z, c) / a), p(t) = 6t^5 - 15t^4 + 10t^3 clamped to
#                 [0, 1] (the "0.6 t^5" in the Supplement is a typo), delta = periodic distance
#   2d_slope      S4.5, Eq. (2DTruth3): rho0(z) = max{0, 30 - 90 f_SN(z; (0.3,0.3), 0.5 I, (-1,-1))}
#   2d_twobump    S4.5, Eq. (2DTruth2): rho0(z) = 6 f_SN(z; (0.3,0.8), 0.03 I, (-1,-1))
#                                               + 14 f_SN(z; (0.7,0.2), 0.05 I, (3,-2))
#   2d_detcov     S4.6, Eq. (2DTruthCoord): rho0(z1, z2) = 15 z2 f_SN(z1; 0.8, 0.3, -5), with Z1
#                 the random field of the 1D study and the deterministic Z2(x) = x1 (paper:
#                 1/2 + x1 on [-1/2,1/2]^2)
#   2d_overparam  S4.6, over-parametrisation: the SAME observations as the main study's 1d_skn
#                 datasets (read from ../Simulation Study/Data), fitted with the extra
#                 deterministic covariate Z2(x) = x1; truth rho0(z1, z2) = 5 f_SN(z1; 0.8, 0.3, -5)
# f_SN is the skew-normal density of package sn (dsn, dmsn with tau = 0); in 2D the second
# argument is the scale matrix Omega.
#
# The deterministic covariate is evaluated at the pixel centres and 16-bit coded like the
# random fields, so it enters the likelihood, the quadrature and the event simulation exactly
# as a random covariate does.
#
# Seeds: data_seed = seed_base + 1000 n + rep, with a separate seed_base per scenario (no overlap
# with the main study, whose seeds are 1e6 d + 1000 n + rep); 2d_overparam reuses the 1d_skn
# seeds of the main study because it reuses its datasets. MCMC seed = data seed + 70000007, as in
# the main study.
###############################################################################

# ---------------------------------------------------------------- ground truths (z: N x d matrix)
# Each truth is self-contained (its definition is stored in the dataset specification).
rho0_1d_negexp <- function(z) 2 * exp(3 * (1 - z[, 1]) - 1)

rho0_1d_minmax <- function(z) {
  smootherstep <- function(t) { t <- pmin(pmax(t, 0), 1); 6 * t^5 - 15 * t^4 + 10 * t^3 }
  plateau <- function(x, c, a) 1 - smootherstep(abs(((x - c + 0.5) %% 1) - 0.5) / (a / 2))
  2 + 2 * plateau(z[, 1], 3 / 4, 3 / 8) - 2 * plateau(z[, 1], 1 / 4, 3 / 8)
}

rho0_2d_slope <- function(z)
  pmax(0, 30 - 90 * sn::dmsn(z[, 1:2, drop = FALSE], xi = c(0.3, 0.3), Omega = diag(0.5, 2),
                             alpha = c(-1, -1)))

rho0_2d_twobump <- function(z)
  6 * sn::dmsn(z[, 1:2, drop = FALSE], xi = c(0.3, 0.8), Omega = diag(0.03, 2), alpha = c(-1, -1)) +
  14 * sn::dmsn(z[, 1:2, drop = FALSE], xi = c(0.7, 0.2), Omega = diag(0.05, 2), alpha = c(3, -2))

rho0_2d_detcov <- function(z) pmax(0, 15 * z[, 2] * sn::dsn(z[, 1], xi = 0.8, omega = 0.3, alpha = -5))

rho0_2d_overparam <- function(z) 5 * sn::dsn(z[, 1], xi = 0.8, omega = 0.3, alpha = -5)

# ---------------------------------------------------------------- scenario catalogue
# type: "gp" = all covariates are random GP fields (main-study generator, get_replicates);
#       "detcov" = random field Z1 plus deterministic Z2(x) = x1 (get_replicates_detcov);
#       "overparam" = main-study 1d_skn dataset augmented with Z2(x) = x1.
ADDSIM_SCENARIOS <- list(
  "1d_negexp" = list(d = 1L, type = "gp", rho0 = rho0_1d_negexp, ls = 0.005,
    seed_base = 11000000L, section = "S4.4", equation = "1DTruth2",
    formula = "2 exp(3 (1 - z) - 1)"),
  "1d_minmax" = list(d = 1L, type = "gp", rho0 = rho0_1d_minmax, ls = 0.005,
    seed_base = 13000000L, section = "S4.4", equation = "1DTruth3",
    formula = "2 + 2 P(z; 3/4, 3/8) - 2 P(z; 1/4, 3/8), smootherstep 6t^5 - 15t^4 + 10t^3"),
  "2d_slope" = list(d = 2L, type = "gp", rho0 = rho0_2d_slope, ls = c(0.005, 0.05),
    seed_base = 21000000L, section = "S4.5", equation = "2DTruth3",
    formula = "max{0, 30 - 90 f_SN(z; (0.3, 0.3), 0.5 I, (-1, -1))}"),
  "2d_twobump" = list(d = 2L, type = "gp", rho0 = rho0_2d_twobump, ls = c(0.005, 0.05),
    seed_base = 23000000L, section = "S4.5", equation = "2DTruth2",
    formula = "6 f_SN(z; (0.3, 0.8), 0.03 I, (-1, -1)) + 14 f_SN(z; (0.7, 0.2), 0.05 I, (3, -2))"),
  "2d_detcov" = list(d = 2L, type = "detcov", rho0 = rho0_2d_detcov, ls = 0.005,
    seed_base = 25000000L, section = "S4.6", equation = "2DTruthCoord",
    formula = "15 z2 f_SN(z1; 0.8, 0.3, -5), Z2(x) = x1"),
  "2d_overparam" = list(d = 2L, type = "overparam", rho0 = rho0_2d_overparam, ls = 0.005,
    seed_base = 1000000L, section = "S4.6", equation = "1DTruth, over-parametrised",
    formula = "5 f_SN(z1; 0.8, 0.3, -5), data of the main-study 1d_skn datasets, Z2(x) = x1"))

ADDSIM_N <- c(50L, 250L, 500L, 1000L)

addsim_truth_spec <- function(scenario) {
  sc <- ADDSIM_SCENARIOS[[scenario]]
  list(name = scenario, section = sc$section, equation = sc$equation, formula = sc$formula)
}

addsim_data_seed <- function(scenario, n, rep)
  as.integer(ADDSIM_SCENARIOS[[scenario]]$seed_base + 1000L * as.integer(n) + as.integer(rep))

# The engine's scenario lookup (repstudy_parse_cells, documentation) sees the new scenarios.
experiment_scenarios <- function() ADDSIM_SCENARIOS

# ---------------------------------------------------------------- deterministic covariate
# Field with value x1 = x-coordinate of the pixel centre (rows = y, columns = x, as all fields).
x_coordinate_field <- function(m = 50) {
  xc <- (seq_len(m) - 0.5) / m
  M <- matrix(rep(xc, each = m), m, m)
  function() M
}

# Dataset cache shared by the custom loader: simulate once, save with checksum, reload after.
# Same file format and locking as get_replicates() in ../R/data_generation.R.
addsim_cached_raw <- function(f, spec, key, n, m, simulate) {
  data_dir <- dirname(f)
  if (!file.exists(f)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    lock <- paste0(f, ".lock"); started <- Sys.time()
    repeat {
      if (dir.create(lock, showWarnings = FALSE)) { on.exit(unlink(lock, recursive = TRUE), add = TRUE); break }
      if (file.exists(f)) break
      if (as.numeric(difftime(Sys.time(), started, units = "secs")) > 45)
        stop("Timed out waiting for dataset creation: ", lock)
      Sys.sleep(0.1)
    }
  }
  if (file.exists(f)) {
    obj <- readRDS(f)
    if (!identical(obj$specification, spec))
      stop("Cached dataset definition differs: ", f, ". Use a new data directory; existing observations were preserved.")
    if (!identical(obj$data_id, object_id(obj$payload))) stop("Dataset checksum failed: ", f)
    raw_data <- lapply(seq_len(n), function(i) list(
      code = matrix(readBin(obj$payload$codes[[i]], "integer", n = m * m * spec$d,
                            size = 2, signed = TRUE, endian = "little") + 32768L, ncol = spec$d),
      pix = obj$payload$pix[[i]], u = obj$payload$u[[i]]))
  } else {
    raw_data <- with_data_seed(spec$seed, simulate())
    payload <- list(
      codes = lapply(raw_data, function(r)
        writeBin(as.vector(r$code) - 32768L, raw(), size = 2, endian = "little")),
      pix = lapply(raw_data, `[[`, "pix"), u = lapply(raw_data, `[[`, "u"))
    obj <- list(specification = spec, key = key, data_id = object_id(payload), payload = payload,
                provenance = list(R = R.version.string, BLAS = extSoftVersion()[["BLAS"]]))
    tmp <- tempfile(tmpdir = data_dir, fileext = ".tmp")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(obj, tmp)
    if (!file.rename(tmp, f)) stop("Could not save dataset: ", f)
  }
  list(raw = raw_data, data_id = obj$data_id)
}

# Random field Z1 (length-scale ls) and deterministic Z2(x) = x1.
get_replicates_detcov <- function(key, n, rho0, ls, seed, m = 50, data_dir, truth_spec) {
  stopifnot(length(key) == 1L, grepl("^[A-Za-z0-9_-]+$", key), length(ls) == 1L)
  spec <- list(schema = DATA_SCHEMA_VERSION, n = as.integer(n), m = as.integer(m), d = 2L,
    covariates = c(sprintf("Z1: GP field, length-scale %g", ls), "Z2: x-coordinate of the pixel centre"),
    ls = as.numeric(ls), seed = as.integer(seed), window = c(0, 1, 0, 1), coding_max = Q16,
    rng = c("Mersenne-Twister", "Inversion", "Rejection"),
    truth = truth_definition(rho0, truth_spec),
    generator = lapply(list(make_covariate_sampler, x_coordinate_field, simulate_raw, pixel_centres),
                       function_definition))
  f <- file.path(data_dir, paste0(key, ".rds"))
  cached <- addsim_cached_raw(f, spec, key, n, m, function()
    simulate_raw(n, rho0, list(make_covariate_sampler(ls, m = m), x_coordinate_field(m)), m))
  reps <- assemble_replicates(cached$raw, m)
  attr(reps, "dataset") <- list(key = key, data_id = cached$data_id, specification = spec, file = f)
  reps
}

# Main-study 1d_skn dataset (identical observations: read from the main study's Data folder when
# present, otherwise simulated from the same seed into this folder's Data), augmented with the
# deterministic covariate Z2(x) = x1 at the quadrature nodes and at the events.
get_replicates_overparam <- function(n, rep, m = 50, data_dir, main_data_dir = NULL) {
  key <- sprintf("1d_skn_n%d_r%d", as.integer(n), as.integer(rep))
  seed <- repstudy_data_seed(1L, n, rep)
  stopifnot(identical(seed, addsim_data_seed("2d_overparam", n, rep)))
  src <- if (!is.null(main_data_dir) && file.exists(file.path(main_data_dir, paste0(key, ".rds"))))
    main_data_dir else data_dir
  reps1 <- get_replicates(key, n, rho0_1d_skn, 0.005, seed = seed, m = m, data_dir = src)
  ds <- attr(reps1, "dataset")
  pix <- readRDS(ds$file)$payload$pix                  # event pixels (file verified just above)
  xq <- round(pixel_centres(m)[, "x"] * Q16) / Q16      # coded like the random fields
  reps <- lapply(seq_along(reps1), function(i) {
    r <- reps1[[i]]
    stopifnot(nrow(r$Zq) == m^2, nrow(r$Zev) == length(pix[[i]]),
              identical(unname(r$Zev[, 1]), unname(r$Zpix[pix[[i]], 1])))
    r$Zq <- unname(cbind(r$Zq, xq))
    r$Zev <- unname(cbind(r$Zev, xq[pix[[i]]]))
    r$Zpix <- unname(cbind(r$Zpix, xq))
    r
  })
  ds$source <- if (identical(src, main_data_dir)) "main-study dataset (Simulation Study/Data)" else
    "simulated from the main-study seed into this folder's Data (main-study file not found)"
  ds$augmentation <- "Z2(x) = x-coordinate of the pixel centre, 16-bit coded"
  attr(reps, "dataset") <- ds
  reps
}

# ---------------------------------------------------------------- engine hook
# Replaces repstudy_load_data() of ../R/replicate_study.R (same return value).
repstudy_load_data <- function(cfg, scenario, n, rep) {
  sc <- ADDSIM_SCENARIOS[[scenario]]
  if (is.null(sc)) stop("Unknown scenario: ", scenario)
  seed <- addsim_data_seed(scenario, n, rep)
  key <- sprintf("%s_n%d_r%d", scenario, as.integer(n), as.integer(rep))
  m <- cfg$science$m
  reps <- switch(sc$type,
    gp = get_replicates(key, n, sc$rho0, sc$ls, seed = seed, m = m, data_dir = cfg$data_dir,
                        truth_spec = addsim_truth_spec(scenario)),
    detcov = get_replicates_detcov(key, n, sc$rho0, sc$ls, seed = seed, m = m,
                                   data_dir = cfg$data_dir, truth_spec = addsim_truth_spec(scenario)),
    overparam = get_replicates_overparam(n, rep, m = m, data_dir = cfg$data_dir,
                                         main_data_dir = cfg$main_data_dir),
    stop("Unknown scenario type: ", sc$type))
  stopifnot(all(vapply(reps, function(r) ncol(r$Zq) == sc$d && ncol(r$Zev) == sc$d, logical(1))))
  list(reps = reps, dataset = attr(reps, "dataset"), sc = sc, seed = seed)
}
