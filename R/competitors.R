###############################################################################
# competitors.R -- frequentist baselines of the simulation study (also used for the kernel
# curves of the wildfire application)
#
# (1) Averaged ratio-form kernel estimator, paper Eqs. (AverageKernel)-(IndividualKernels): the
#     plain average over replicates of the per-replicate ratio estimates, Gaussian kernel, with
#     bandwidth selected by an oracle (minimal L2 error to rho0, simulations only) or by
#     leave-one-replicate-out Poisson likelihood cross-validation.
# (2) Log-polynomial Poisson intensity models of total degree 1-4, with all interactions in 2D,
#     fitted by exact Poisson maximum likelihood on the pixel quadrature.
#
# Consolidated 12 October 2026 from competitors.R (version 2.1), run_kernel_grid_addon.R and
# run_polynomial_addon.R (originals in archive/). Changes with respect to version 2.1:
#   * extended default bandwidth grids and finer bins, so that the selected bandwidths are
#     interior (in 2D the old minimum 0.02 was often selected at large n):
#       1D: 1000 bins, 25 log-spaced bandwidths in [0.003, 1]   (was 500 bins, 20 in [0.01, 1]);
#       2D: 160 x 160 bins, 12 x 12 log-spaced (h1, h2) in [0.008, 2]^2 (was 80 x 80, 8 x 8 in [0.02, 1]^2);
#   * vectorised separable 2D smoothing (conv2_all), checked on 11 October 2026 to reproduce the
#     version-2.1 kernel scores and curves to 1e-14;
#   * full log-polynomials of degree 1-4 (the additive log-cubic is dropped; in 1D full and
#     additive designs coincide);
#   * spatstat::mppm and inlabru/INLA variants removed (not used in the paper).
###############################################################################

KERNEL_GRID <- list(
  d1 = list(G = 1000L, h = exp(seq(log(0.003), log(1), length.out = 25))),
  d2 = list(G = 160L,  h = exp(seq(log(0.008), log(2), length.out = 12))))

# ---------------------------------------------------------------- kernel estimator
# Covariate values of the quadrature nodes (denominator: spatial occupation measure) and of the
# events (numerator) are binned on G cells per axis; cell id = (i1 - 1) G + i2 in 2D.
.bin_counts <- function(reps, d, G) {
  cell <- function(Z) { k <- pmin(floor(Z * G), G - 1); if (d == 1) k[, 1] + 1 else k[, 1] * G + k[, 2] + 1 }
  n <- length(reps); nc <- G^d
  Hd <- Matrix::sparseMatrix(i = unlist(lapply(seq_len(n), function(i) rep(i, nrow(reps[[i]]$Zq)))),
                             j = unlist(lapply(reps, function(r) cell(r$Zq))),
                             x = unlist(lapply(reps, function(r) rep(r$qw, nrow(r$Zq)))), dims = c(n, nc))
  Ho <- Matrix::sparseMatrix(i = unlist(lapply(seq_len(n), function(i) rep(i, nrow(reps[[i]]$Zev)))),
                             j = unlist(lapply(reps, function(r) cell(r$Zev))), x = 1, dims = c(n, nc))
  list(Hd = as.matrix(Hd), Ho = as.matrix(Ho), G = G, d = d,
       centres = if (d == 1) matrix((seq_len(G) - 0.5) / G, ncol = 1) else
         { c <- (seq_len(G) - 0.5) / G; cbind(rep(c, each = G), rep(c, times = G)) })
}

.gauss_matrix <- function(G, h) { c <- (seq_len(G) - 0.5) / G; exp(-0.5 * (outer(c, c, "-") / h)^2) }

# Separable 2D Gaussian smoothing of every row of H (n x G^2, cell id (i1-1)*G + i2) at once:
# row r -> K1 %*% M_r %*% t(K2), returned in the same cell ordering.
conv2_all <- function(H, K1, K2, G) {
  n <- nrow(H)
  Y <- H; dim(Y) <- c(n * G, G)               # rows (r, i2), columns i1
  Y <- Y %*% t(K1)                            # [r, i2, a]
  dim(Y) <- c(n, G, G); Y <- aperm(Y, c(1, 3, 2)); dim(Y) <- c(n * G, G)   # rows (r, a), cols i2
  Y <- Y %*% t(K2)                            # [r, a, b]
  dim(Y) <- c(n, G, G); Y <- aperm(Y, c(1, 3, 2)); dim(Y) <- c(n, G * G)   # column (a-1)*G + b
  Y
}

# n x G^d matrix of per-replicate ratio estimates for bandwidth h (scalar in 1D, (h1, h2) in 2D)
.kernel_per_replicate <- function(bc, h) {
  if (bc$d == 1) {
    K <- .gauss_matrix(bc$G, h)
    return(tcrossprod(bc$Ho, K) / pmax(tcrossprod(bc$Hd, K), 1e-300))
  }
  K1 <- .gauss_matrix(bc$G, h[1]); K2 <- .gauss_matrix(bc$G, h[2])
  conv2_all(bc$Ho, K1, K2, bc$G) / pmax(conv2_all(bc$Hd, K1, K2, bc$G), 1e-300)
}

# Averaged ratio-form estimator with oracle and leave-one-replicate-out CV bandwidth selection.
# rho0: true intensity (oracle selection). On real data pass any placeholder, e.g.
# function(z) rep(1, nrow(z)), and use only $cv. Returns both selections on the bin centres.
kernel_competitor <- function(reps, d, rho0, h_grid = NULL, G = NULL) {
  spec <- KERNEL_GRID[[paste0("d", d)]]
  if (is.null(G)) G <- spec$G
  if (is.null(h_grid))
    h_grid <- if (d == 1) as.list(spec$h) else asplit(as.matrix(expand.grid(spec$h, spec$h)), 1)
  bc <- .bin_counts(reps, d, G)
  truth <- rho0(bc$centres); n <- length(reps)
  best_o <- best_c <- NULL; scores <- vector("list", length(h_grid))
  for (k in seq_along(h_grid)) {
    h <- unname(unlist(h_grid[[k]]))
    R <- .kernel_per_replicate(bc, h)
    rbar <- colMeans(R)
    loo <- (n * matrix(rbar, n, ncol(R), byrow = TRUE) - R) / (n - 1)   # leave-one-replicate-out
    cvll <- rowSums(bc$Ho * log(pmax(loo, 1e-10))) - rowSums(bc$Hd * loo)
    rm(R, loo)
    l2 <- rel_L2(rbar, truth); cv <- mean(cvll)
    scores[[k]] <- c(h = h, relL2 = l2, cv = cv)
    if (is.null(best_o) || l2 < best_o$l2) best_o <- list(h = h, l2 = l2, cv = cv, rbar = rbar)
    if (is.null(best_c) || cv > best_c$cv) best_c <- list(h = h, l2 = l2, cv = cv, rbar = rbar)
  }
  hs <- if (d == 1) unlist(h_grid) else do.call(rbind, h_grid)
  at_edge <- function(h) if (d == 1) any(h %in% range(hs)) else
    any(h[1] %in% range(hs[, 1]), h[2] %in% range(hs[, 2]))
  best_o$at_edge <- at_edge(best_o$h); best_c$at_edge <- at_edge(best_c$h)
  list(oracle = best_o, cv = best_c, centres = bc$centres, truth = truth, G = G,
       scores = as.data.frame(do.call(rbind, scores)))
}

# Kernel estimates (oracle, CV) at arbitrary points z in [0,1]^d: linear (1D) or bilinear (2D)
# interpolation between bin centres, constant beyond the outermost centres.
kernel_on_grid <- function(kc, z) {
  z <- as.matrix(z); d <- ncol(z)
  if (d == 1) {
    x <- kc$centres[, 1]
    return(cbind(approx(x, kc$oracle$rbar, xout = z[, 1], rule = 2)$y,
                 approx(x, kc$cv$rbar, xout = z[, 1], rule = 2)$y))
  }
  g <- sort(unique(kc$centres[, 1])); r <- range(g)
  stopifnot(identical(unname(kc$centres), unname(cbind(rep(g, each = length(g)), rep(g, times = length(g))))))
  q <- (pmin(pmax(z, r[1]), r[2]) - r[1]) / diff(r)
  A <- interp_matrix(list(d = 2, type = "grid", g1 = g, g2 = g), q)   # from gp_sampler.R
  as.matrix(A %*% cbind(kc$oracle$rbar, kc$cv$rbar))
}

# ---------------------------------------------------------------- log-polynomial Poisson models
# Intercept plus all monomials of total degree 1..degree (interactions included in 2D).
poly_design <- function(Z, degree) {
  Z <- as.matrix(Z)
  if (ncol(Z) == 1)
    return(cbind(1, do.call(cbind, lapply(seq_len(degree), function(r) Z[, 1]^r))))
  ex <- expand.grid(a = 0:degree, b = 0:degree)
  ex <- ex[ex$a + ex$b >= 1 & ex$a + ex$b <= degree, ]
  ex <- ex[order(ex$a + ex$b, -ex$a), ]
  cbind(1, do.call(cbind, lapply(seq_len(nrow(ex)), function(k) Z[, 1]^ex$a[k] * Z[, 2]^ex$b[k])))
}

# Exact Poisson ML on the pixel quadrature: maximise sum_k eta(z_k) - sum_q w_q exp(eta(z_q)),
# eta = poly_design %*% b, by Newton-Raphson with step halving. Returns the fitted intensity at
# z_eval and convergence information.
poly_competitor <- function(reps, z_eval, degree, max_iter = 100L) {
  Xq <- poly_design(do.call(rbind, lapply(reps, `[[`, "Zq")), degree)
  wq <- unlist(lapply(reps, function(r) rep(r$qw, nrow(r$Zq))))
  Xe <- poly_design(do.call(rbind, lapply(reps, `[[`, "Zev")), degree)
  se <- colSums(Xe)
  obj <- function(b) sum(Xe %*% b) - sum(wq * exp(Xq %*% b))
  b <- c(log(nrow(Xe) / sum(wq)), rep(0, ncol(Xq) - 1))
  converged <- FALSE; iter <- 0L
  for (iter in seq_len(max_iter)) {
    mu <- as.vector(wq * exp(Xq %*% b))
    step <- solve(crossprod(Xq, Xq * mu), se - colSums(Xq * mu))
    t <- 1; f0 <- obj(b)
    while (obj(b + t * step) < f0 - 1e-12 && t > 1e-8) t <- t / 2
    b <- b + t * step
    if (max(abs(t * step)) < 1e-9) { converged <- TRUE; break }
  }
  list(curve = as.vector(exp(poly_design(z_eval, degree) %*% b)), coef = as.vector(b),
       iterations = as.integer(iter), converged = converged, loglik = obj(b))
}
