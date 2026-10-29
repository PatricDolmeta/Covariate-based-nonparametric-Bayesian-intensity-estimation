###############################################################################
# data_generation.R
# Covariate fields, ground truths and replicated point patterns.
#
# Conventions (all scripts in 'Revised Code'):
#   * observation window W = [0,1]^2, discretised into an m x m pixel grid;
#   * a covariate field is stored as an m x m matrix (rows = y, columns = x, as in
#     spatstat 'im'), pixel vector = as.vector(matrix) (column-major);
#   * a replicate is a list with
#       Zq    : (#quadrature nodes) x d matrix of covariate values at the quadrature nodes
#       qw    : quadrature weight per node (= pixel area)
#       Zev   : K x d matrix of covariate values at the K events
#       xy    : K x 2 matrix of event locations (needed only by spatstat-based competitors)
#       Zpix  : full-resolution pixel covariates (used by the exact sampler)
#       Zim   : list of d spatstat 'im' objects (only if keep_im = TRUE)
#
# Reproducibility (6 Oct 2026):
#   * covariate fields are drawn with a Cholesky factor (plus a 1e-6 diagonal jitter) instead of an
#     eigendecomposition: the eigenvectors of the near-singular SE covariance are not unique, so
#     different BLAS/LAPACK libraries (e.g. R's reference BLAS vs vecLib) gave different fields
#     from the same seed; the Cholesky factor is unique and practically BLAS-independent;
#   * fields are stored as 16-bit codes (resolution 1/65535) and the events are simulated from
#     the coded values, so a saved data set is exact;
#   * get_replicates() simulates a data set once, saves it under Data/<key>.rds and reloads
#     it afterwards, so every method and every rerun uses exactly the same data.
###############################################################################

# ---------------------------------------------------------------- covariates
# Square-exponential GP on the pixel grid, covariance exp(-|x-x'|^2 / (2 ls)) (same
# parametrisation as in the old code: ls = 0.005 for Z1, 0.05 for Z2). The SE kernel is
# separable, so the field is drawn as B E B^T with B B^T = 1D covariance matrix (+ jitter).
make_covariate_sampler <- function(ls, m = 50, jitter = 1e-6) {
  xc <- (seq_len(m) - 0.5) / m
  K1 <- exp(-outer(xc, xc, "-")^2 / (2 * ls))
  B <- t(chol(K1 + diag(jitter, m)))                 # deterministic factor up to roundoff
  function() {
    X <- B %*% matrix(rnorm(m * m), m, m) %*% t(B)
    # Var(X[y,x]) = (1 + jitter)^2. Population standardisation preserves
    # the Gaussian-CDF marginal construction; never standardise each realisation.
    pnorm(X / (1 + jitter))
  }
}

pixel_centres <- function(m) {
  xc <- (seq_len(m) - 0.5) / m
  # column-major order of an m x m matrix with rows = y, cols = x
  cbind(x = rep(xc, each = m), y = rep(xc, times = m))
}

# ---------------------------------------------------------------- ground truths
# z is always an N x d matrix
rho0_1d_skn <- function(z) 5 * sn::dsn(z[, 1], xi = 0.8, omega = 0.3, alpha = -5)   # ||rho0||_2 = 6.42

rho0_1d_sin <- function(z) 15 + 12 * sin(2 * pi * 5 * z[, 1])                         # sensitivity truth

# Anisotropic bump/hole truth. NOTE: the bumps are *unnormalised* Gaussians; this reproduces
# ||rho0||_2 = 10.8 reported in the paper's table (the paper's display uses the normal p.d.f.
# f_N, which would give ||rho0||_2 = 18.0 and needs correcting). The old code had amplitude 100.
rho0_2d_hole <- function(z) {
  g <- function(m1, m2) exp(-0.5 * (((z[, 1] - m1) / 0.08)^2 + ((z[, 2] - m2) / 0.5)^2))
  pmax(0, 10 - 10 * g(0.8, 0.3) + 10 * g(0.3, 0.8))
}

# ---------------------------------------------------------------- replicates
Q16 <- 65535L                                        # 16-bit coding of covariate values in [0, 1]
DATA_SCHEMA_VERSION <- 2L
DATA_DIR <- "Data"

# A single catalog for all experiment entry points. The observation window is [0,1]^2
# in both scenarios; d is the number of covariates, not the spatial dimension.
experiment_scenarios <- function() list(
  "1d_skn" = list(d = 1L, rho0 = rho0_1d_skn, ls = 0.005),
  "2d_hole" = list(d = 2L, rho0 = rho0_2d_hole, ls = c(0.005, 0.05)))

object_id <- function(x) {
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f))
  writeBin(serialize(x, NULL, version = 2), f)
  unname(tools::md5sum(f))
}

function_definition <- function(f) {
  # Exclude source locations, which differ between RStudio and Rscript.
  definition <- parse(text = deparse(f, width.cutoff = 500L,
    control = c("keepNA", "keepInteger", "niceNames")), keep.source = FALSE)[[1]]
  list(formals = definition[[2]], body = definition[[3]])
}

# Cached replicate-study experiments use self-contained, registered truths. Custom truths
# must provide an explicit list of all their parameters through truth_spec.
truth_definition <- function(rho0, truth_spec = NULL) {
  registered <- list("1d_skn" = rho0_1d_skn, "1d_sin" = rho0_1d_sin, "2d_hole" = rho0_2d_hole)
  if (is.null(truth_spec)) {
    matches <- vapply(registered, function(f)
      identical(function_definition(f), function_definition(rho0)), logical(1))
    if (!any(matches)) stop("Unregistered truth: supply truth_spec with its name and every parameter.")
    truth_spec <- list(name = names(registered)[which(matches)[1]])
  }
  list(parameters = truth_spec, implementation = function_definition(rho0),
       sn_version = if (identical(truth_spec$name, "1d_skn")) as.character(utils::packageVersion("sn")) else NULL)
}

# A saved dataset is the reproducibility authority. A seed alone is not a guarantee
# of identical floating-point calculations across R versions or BLAS backends.
dataset_specification <- function(n, rho0, ls, seed, m, truth_spec = NULL) {
  stopifnot(length(n) == 1L, n >= 1, n == as.integer(n), length(m) == 1L,
            m >= 2, m == as.integer(m), is.numeric(ls), length(ls) %in% 1:2,
            all(is.finite(ls)), all(ls > 0), length(seed) == 1L, is.finite(seed))
  list(schema = DATA_SCHEMA_VERSION, n = as.integer(n), m = as.integer(m),
       d = as.integer(length(ls)), ls = as.numeric(ls), seed = as.integer(seed),
       window = c(0, 1, 0, 1), coding_max = Q16,
       rng = c("Mersenne-Twister", "Inversion", "Rejection"),
       truth = truth_definition(rho0, truth_spec),
       generator = lapply(list(make_covariate_sampler, simulate_raw, pixel_centres), function_definition))
}

# Preserve the caller's RNG state: MCMC seeds must be set independently afterwards.
with_data_seed <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(code)
}

set_mcmc_seed <- function(seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
}

# Raw simulation of n replicates on the m x m pixel grid: 16-bit covariate codes (m^2 x d integer
# matrix), event pixel indices and the uniform within-pixel offsets of the events.
# 'ls' = vector of covariate length-scales (one per covariate), or a list of samplers from
# make_covariate_sampler() (as in the older scripts).
simulate_raw <- function(n, rho0, ls, m = 50) {
  samplers <- if (is.list(ls) && is.function(ls[[1]])) ls else lapply(ls, make_covariate_sampler, m = m)
  d <- length(samplers); area <- 1 / m^2
  lapply(seq_len(n), function(i) {
    code <- matrix(sapply(samplers, function(s) as.integer(round(as.vector(s()) * Q16))), ncol = d)
    Z <- code / Q16
    counts <- rpois(nrow(Z), rho0(Z) * area)
    pix <- rep(seq_len(nrow(Z)), counts)
    list(code = code, pix = pix, u = matrix(runif(2 * length(pix), -0.5, 0.5), ncol = 2))
  })
}

# Build the replicates used by the samplers from raw data.
# quad_m: resolution of the quadrature grid used by the sampler (must divide m); the events are
# always simulated from the full-resolution (m x m) piecewise-constant intensity, so for
# quad_m < m the likelihood integral is genuinely approximated (sensitivity analysis).
assemble_replicates <- function(raw, m = 50, quad_m = m, keep_im = FALSE) {
  pc <- pixel_centres(m)
  stopifnot(m %% quad_m == 0)
  step <- m / quad_m
  sub <- seq(ceiling(step / 2), m, by = step)        # pixel closest to each coarse-cell centre
  sub_idx <- as.vector(outer(sub, (sub - 1) * m, "+"))   # = as.vector(M[sub, sub]) indices
  xc <- (seq_len(m) - 0.5) / m
  lapply(raw, function(r) {
    Z <- r$code / Q16; d <- ncol(Z)
    Zq <- if (step == 1) Z else Z[sub_idx, , drop = FALSE]
    rep_i <- list(Zq = Zq, qw = 1 / quad_m^2, Zev = Z[r$pix, , drop = FALSE],
                  xy = pc[r$pix, , drop = FALSE] + r$u / m, Zpix = Z)
    if (keep_im)
      rep_i$Zim <- lapply(seq_len(d), function(j)
        spatstat.geom::im(matrix(Z[, j], m, m), xcol = xc, yrow = xc))
    rep_i
  })
}

# One-shot simulation without saving (as before): 'samplers' = list of samplers or length-scales.
simulate_replicates <- function(n, rho0, samplers, m = 50, quad_m = m, keep_im = FALSE)
  assemble_replicates(simulate_raw(n, rho0, samplers, m), m, quad_m, keep_im)

# Simulate once, save, reload afterwards. Files: <data_dir>/<key>.rds (fields as 2-byte codes).
get_replicates <- function(key, n, rho0, ls, seed, m = 50, quad_m = m, keep_im = FALSE,
                           data_dir = DATA_DIR, save = TRUE, truth_spec = NULL) {
  stopifnot(length(key) == 1L, grepl("^[A-Za-z0-9_-]+$", key), m %% quad_m == 0)
  spec <- dataset_specification(n, rho0, ls, seed, m, truth_spec)
  f <- file.path(data_dir, paste0(key, ".rds"))
  # Serialise creation of a shared key: sensitivity workers may request it together.
  if (save && !file.exists(f)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    lock <- paste0(f, ".lock")
    started <- Sys.time()
    repeat {
      if (dir.create(lock, showWarnings = FALSE)) {
        on.exit(unlink(lock, recursive = TRUE), add = TRUE)
        break
      }
      if (file.exists(f)) break
      if (as.numeric(difftime(Sys.time(), started, units = "secs")) > 45)
        stop("Timed out waiting for dataset creation: ", lock)
      Sys.sleep(0.1)
    }
  }
  if (save && file.exists(f)) {
    obj <- readRDS(f)
    if (!identical(obj$specification, spec))
      stop("Cached dataset definition differs or lacks version-2 metadata: ", f,
           ". Use a new key or data directory; existing observations were preserved.")
    if (!identical(obj$data_id, object_id(obj$payload))) stop("Dataset checksum failed: ", f)
    raw_data <- lapply(seq_len(n), function(i) list(
      code = matrix(readBin(obj$payload$codes[[i]], "integer", n = m * m * spec$d,
                            size = 2, signed = TRUE, endian = "little") + 32768L, ncol = spec$d),
      pix = obj$payload$pix[[i]], u = obj$payload$u[[i]]))
  } else {
    raw_data <- with_data_seed(seed, simulate_raw(n, rho0, ls, m))
    payload <- list(
      codes = lapply(raw_data, function(r)
        writeBin(as.vector(r$code) - 32768L, raw(), size = 2, endian = "little")),
      pix = lapply(raw_data, `[[`, "pix"), u = lapply(raw_data, `[[`, "u"))
    obj <- list(specification = spec, key = key, data_id = object_id(payload), payload = payload,
                provenance = list(R = R.version.string, BLAS = extSoftVersion()[["BLAS"]]))
    if (save) {
      tmp <- tempfile(tmpdir = data_dir, fileext = ".tmp")
      on.exit(unlink(tmp), add = TRUE)
      saveRDS(obj, tmp)
      if (!file.rename(tmp, f)) stop("Could not save dataset: ", f)
    }
  }
  reps <- assemble_replicates(raw_data, m, quad_m, keep_im)
  attr(reps, "dataset") <- list(key = key, data_id = obj$data_id,
                                specification = spec, file = if (save) f else NA_character_)
  reps
}

n_events <- function(reps) sum(sapply(reps, function(r) nrow(r$Zev)))

# spatstat objects for the competitors (requires keep_im = TRUE)
as_ppp_list <- function(reps)
  lapply(reps, function(r) spatstat.geom::ppp(r$xy[, 1], r$xy[, 2], window = spatstat.geom::owin()))

# evaluation grid on the covariate space (midpoint rule for L2 errors)
eval_grid <- function(d, k = if (d == 1) 500 else 60) {
  c <- (seq_len(k) - 0.5) / k
  if (d == 1) matrix(c, ncol = 1) else as.matrix(expand.grid(z1 = c, z2 = c))
}

rel_L2 <- function(est, truth) sqrt(mean((est - truth)^2)) / sqrt(mean(truth^2))
