###############################################################################
# gp_sampler.R -- the GP sampler shared by the simulation study and the wildfire application
# (consolidated 12 October 2026; formerly Code for Data Analysis/.../R/gp_pcn_sampler_nd.R).
#
# Extension of the original simulation sampler (archive/simulation_R/gp_pcn_sampler.R) to d = 3:
#   * make_nodes(3, type = "grid", grid_m = m): m^3 tensor grid, node id
#     ((k1 - 1) m + (k2 - 1)) m + k3 (z1 slowest, z3 fastest), as in the 2D convention;
#   * trilinear interpolation (interp_matrix) and 3D covariate binning (bin_data);
#   * Kronecker factor L1 x L2 x L3 of (C1 + jI) x (C2 + jI) x (C3 + jI): products and solves
#     are applied mode-wise (kron_apply), log det from the factor diagonals (flogdiag);
#   * optional storage of thinned node draws (w, rho*) after burn-in (save_nodes_every), from
#     which any functional of rho can be computed afterwards (spatial maps, slices, bands).
# For d = 1 and the 2D tensor grid the code paths, random-number stream and output are those
# of the original sampler (checked on 11 October 2026 by archive/application/verify_nd_sampler.R).
# Default prior (12 October 2026): a_ell = 1, b_ell = 2.5, lambda = 0.1 (exponential gamma part
# of the tilted hyper-prior), rho* ~ Gamma(1, 0.1) truncated to [0, 100 + log n]. Callers in
# this repository always pass the prior explicitly. Documentation of the original file follows.
###############################################################################
# gp_pcn_sampler.R
# Metropolis-within-Gibbs sampler for the revised prior (quadrature likelihood).
#
#   rho(z)  = rho* sigmoid(w(z)),
#   rho*    ~ Gamma(a_rho, b_rho) truncated to [0, c_rho + log n],
#   w | ell ~ GP(0, k_ell),   k_ell(z, z') = exp(- sum_j ell_j^2 (z_j - z'_j)^2)   [paper, Eq. (ARD)]
#   ell     ~ Pi(ell) propto prod_j ell_j^(a_ell-1) e^(-b_ell ell_j) * exp(-lambda prod_j (1+ell_j))
#
# w is discretised on V nodes and linearly interpolated (1D: piecewise linear; 2D: either
# barycentric on a triangular mesh, or bilinear on an m x m tensor grid), as in the paper's
# Eq. (Discret). On a tensor grid the ARD covariance is a Kronecker product C1(ell1) x C2(ell2),
# so its Cholesky factor is L1 x L2 and every ell-update costs O(m^3) instead of O(m^6)
# (make_nodes(2, type = "grid")); this removes the 2D bottleneck. The likelihood integral uses the
# quadrature nodes of each replicate. For speed, covariate values are binned on a fine G^d grid
# (G = NULL disables binning, i.e. exact use of every quadrature node and event); the
# interpolation matrices are precomputed, so a likelihood evaluation costs O(#bins) and not
# O(n x #pixels).
#
# One iteration:
#   1. n_w_steps (default 10) pCN updates of w | ell, rho*  (step beta, Robbins-Monro tuned in
#      burn-in to ~30%). Several cheap w-moves per (expensive) ell-move greatly improve mixing of ell.
#   2. rho* | w                            (truncated gamma, exact)
#   3. ell_j | w, ell_-j, j = 1..d,        coordinatewise and *sequential* log-random-walk MH;
#      "noncentred" (default): w = L_ell xi with xi held fixed (acceptance uses the likelihood);
#      "centred": target N(w; 0, C_ell) Pi(ell), as in the paper's Eq. (FullCondEll).
#      Step sizes Robbins-Monro tuned in burn-in to ~44%; all adaptation frozen after burn-in at the
#      average over the second half of burn-in (removes chain-to-chain differences in the frozen step).
#
# Changes w.r.t. the old Adaptive_Anisotropic_GP.R: new hyperprior (no theta block), kernel in
# the paper's parametrisation, sequential (valid) ell updates, correctly scaled adaptation,
# truncated rho* prior, linear interpolation instead of nearest node, O(#bins) likelihood.
###############################################################################

suppressPackageStartupMessages(library(Matrix))

# ---------------------------------------------------------------- discretisation
make_nodes <- function(d, V = 200, max_area = 0.0014, type = c("mesh", "grid"), grid_m = 20) {
  type <- match.arg(type)
  if (d == 1) return(list(d = 1, type = "line", nodes = matrix(seq(0, 1, length.out = V), ncol = 1)))
  if (d == 3) {                            # 3D: tensor grid only; node id = ((k1-1) m + (k2-1)) m + k3
    if (type != "grid") stop("d = 3 requires type = \"grid\"")
    g <- seq(0, 1, length.out = grid_m)
    return(list(d = 3, type = "grid", g1 = g, g2 = g, g3 = g,
                nodes = cbind(rep(g, each = grid_m^2), rep(rep(g, each = grid_m), times = grid_m),
                              rep(g, times = grid_m^2))))
  }
  if (d != 2) stop("d must be 1, 2 or 3")
  if (type == "grid") {                    # tensor grid; node id = (k1 - 1) * m + k2 (z1 slow)
    g <- seq(0, 1, length.out = grid_m)
    return(list(d = 2, type = "grid", g1 = g, g2 = g,
                nodes = cbind(rep(g, each = grid_m), rep(g, times = grid_m))))
  }
  sq <- RTriangle::pslg(P = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
  mesh <- RTriangle::triangulate(sq, a = max_area, j = TRUE)
  list(d = 2, type = "mesh", nodes = mesh$P, tri = mesh$T)
}

# sparse interpolation matrix: rows = query points z (N x d), columns = nodes
interp_matrix <- function(nd, z) {
  N <- nrow(z)
  if (nd$d == 1) {
    V <- nrow(nd$nodes)
    t <- z[, 1] * (V - 1); k <- pmin(floor(t), V - 2); f <- t - k
    return(sparseMatrix(i = rep(seq_len(N), 2), j = c(k + 1, k + 2), x = c(1 - f, f),
                        dims = c(N, V)))
  }
  if (nd$d == 3) {                         # trilinear interpolation on the 3D tensor grid
    if (!identical(nd$type, "grid")) stop("d = 3 requires a tensor grid")
    m <- c(length(nd$g1), length(nd$g2), length(nd$g3))
    k <- f <- matrix(0, N, 3)
    for (a in 1:3) {
      ta <- z[, a] * (m[a] - 1); k[, a] <- pmin(floor(ta), m[a] - 2); f[, a] <- ta - k[, a]
    }
    j <- x <- vector("list", 8); h <- 0
    for (c1 in 0:1) for (c2 in 0:1) for (c3 in 0:1) {
      h <- h + 1
      j[[h]] <- ((k[, 1] + c1) * m[2] + (k[, 2] + c2)) * m[3] + (k[, 3] + c3) + 1
      x[[h]] <- (if (c1) f[, 1] else 1 - f[, 1]) * (if (c2) f[, 2] else 1 - f[, 2]) *
                (if (c3) f[, 3] else 1 - f[, 3])
    }
    return(sparseMatrix(i = rep(seq_len(N), 8), j = unlist(j), x = unlist(x), dims = c(N, prod(m))))
  }
  if (identical(nd$type, "grid")) {        # bilinear interpolation on the tensor grid
    m1 <- length(nd$g1); m2 <- length(nd$g2)
    t1 <- z[, 1] * (m1 - 1); k1 <- pmin(floor(t1), m1 - 2); f1 <- t1 - k1
    t2 <- z[, 2] * (m2 - 1); k2 <- pmin(floor(t2), m2 - 2); f2 <- t2 - k2
    j <- c(k1 * m2 + k2, k1 * m2 + k2 + 1, (k1 + 1) * m2 + k2, (k1 + 1) * m2 + k2 + 1) + 1
    x <- c((1 - f1) * (1 - f2), (1 - f1) * f2, f1 * (1 - f2), f1 * f2)
    return(sparseMatrix(i = rep(seq_len(N), 4), j = j, x = x, dims = c(N, m1 * m2)))
  }
  ts <- geometry::tsearch(nd$nodes[, 1], nd$nodes[, 2], nd$tri, z[, 1], z[, 2], bary = TRUE)
  if (anyNA(ts$idx)) stop("query points outside the mesh")
  sparseMatrix(i = rep(seq_len(N), 3), j = as.vector(nd$tri[ts$idx, ]), x = as.vector(ts$p),
               dims = c(N, nrow(nd$nodes)))
}

# ---------------------------------------------------------------- binning of the data
bin_data <- function(reps, d, G = if (d == 1) 2000 else if (d == 2) 200 else 100) {
  Zq <- do.call(rbind, lapply(reps, `[[`, "Zq"))
  qw <- unlist(lapply(reps, function(r) rep(r$qw, nrow(r$Zq))))
  Zev <- do.call(rbind, lapply(reps, `[[`, "Zev"))
  if (is.null(G)) return(list(z_dom = Zq, h_dom = qw, z_obs = Zev, h_obs = rep(1, nrow(Zev)),
                              n = length(reps), N = nrow(Zev)))
  cell <- function(Z) {
    k <- pmin(floor(Z * G), G - 1)
    if (d == 3) return(as.integer((k[, 1] * G + k[, 2]) * G + k[, 3] + 1))   # G^3 < 2^31
    if (d == 1) k[, 1] + 1 else k[, 1] * G + k[, 2] + 1
  }
  ctr <- (seq_len(G) - 0.5) / G
  centre <- function(id) {
    if (d == 1) return(matrix(ctr[id], ncol = 1))
    if (d == 3) return(cbind(ctr[(id - 1) %/% G^2 + 1], ctr[((id - 1) %/% G) %% G + 1],
                             ctr[(id - 1) %% G + 1]))
    cbind(ctr[(id - 1) %/% G + 1], ctr[(id - 1) %% G + 1])
  }
  hd <- rowsum(qw, cell(Zq)); ho <- table(cell(Zev))
  list(z_dom = centre(as.integer(rownames(hd))), h_dom = as.vector(hd),
       z_obs = centre(as.integer(names(ho))), h_obs = as.vector(ho),
       n = length(reps), N = nrow(Zev))
}

# ---------------------------------------------------------------- prior pieces
ard_chol <- function(D2, ell, jitter = 1e-6) {
  # Match the tensor factor convention for every discretisation:
  # product_j { exp(-ell_j^2 D_j) + jitter * 1[D_j = 0] }.
  # On a tensor grid this is exactly (C1 + jitter I) tensor (C2 + jitter I).
  # It differs from adding one scalar nugget after forming the full covariance.
  C <- Reduce(`*`, Map(function(D, l) exp(-l^2 * D) + jitter * (D == 0), D2, ell))
  t(chol(C))                                           # lower triangular
}
# Cholesky factor of the prior covariance on the nodes, as an object: dense (1D, mesh) or
# Kronecker (tensor grid; vectors in node order = row-major vec of the m1 x m2 node array).
make_factor <- function(nd, D2, ell, jitter = 1e-6) {
  if (identical(nd$type, "grid") && !identical(nd$factorisation, "dense")) {
    ch <- function(g, l) { C <- exp(-l^2 * outer(g, g, "-")^2); diag(C) <- diag(C) + jitter; t(chol(C)) }
    if (nd$d == 3)
      return(list(kron = TRUE, Ls = list(ch(nd$g1, ell[1]), ch(nd$g2, ell[2]), ch(nd$g3, ell[3])),
                  ms = c(length(nd$g1), length(nd$g2), length(nd$g3))))
    return(list(kron = TRUE, L1 = ch(nd$g1, ell[1]), L2 = ch(nd$g2, ell[2]),
                m1 = length(nd$g1), m2 = length(nd$g2)))
  }
  if (is.null(D2))
    D2 <- lapply(seq_len(nd$d), function(j) outer(nd$nodes[, j], nd$nodes[, j], "-")^2)
  list(kron = FALSE, L = ard_chol(D2, ell, jitter))
}
# Mode-wise application of a Kronecker factor list (d = 3). A node vector, ordered with z1
# slowest, is the column-major array X[k3, k2, k1]; factor a acts on array mode d + 1 - a.
mode_apply <- function(X, k, f) {
  dims <- dim(X); perm <- c(k, seq_along(dims)[-k])
  Y <- f(matrix(aperm(X, perm), dims[k]))
  aperm(array(Y, c(nrow(Y), dims[-k])), order(perm))
}
kron_apply <- function(F, x, f) {
  m <- F$ms; X <- array(x, rev(m))
  for (a in seq_along(m)) X <- mode_apply(X, length(m) + 1 - a, function(M) f(F$Ls[[a]], M))
  as.vector(X)
}
fmul <- function(F, x) {                    # L x
  if (!is.null(F$Ls)) return(kron_apply(F, x, function(L, M) L %*% M))
  if (!F$kron) return(as.vector(F$L %*% x))
  as.vector(t(F$L1 %*% matrix(x, F$m1, F$m2, byrow = TRUE) %*% t(F$L2)))
}
fsolve <- function(F, w) {                  # L^{-1} w
  if (!is.null(F$Ls)) return(kron_apply(F, w, function(L, M) forwardsolve(L, M)))
  if (!F$kron) return(forwardsolve(F$L, w))
  X <- forwardsolve(F$L1, matrix(w, F$m1, F$m2, byrow = TRUE))
  as.vector(forwardsolve(F$L2, t(X)))
}
flogdiag <- function(F) {                   # sum(log(diag(L))) = log det(C) / 2
  if (!is.null(F$Ls)) return(sum(vapply(seq_along(F$Ls), function(a)
    prod(F$ms[-a]) * sum(log(diag(F$Ls[[a]]))), numeric(1))))
  if (!F$kron) sum(log(diag(F$L))) else
    F$m2 * sum(log(diag(F$L1))) + F$m1 * sum(log(diag(F$L2)))
}

log_prior_ell <- function(ell, a, b, lambda) sum((a - 1) * log(ell) - b * ell) - lambda * prod(1 + ell)

rtgamma <- function(shape, rate, upper) {            # Gamma(shape, rate) truncated to [0, upper]
  Fu <- pgamma(upper, shape, rate)
  if (Fu < 1e-300) return(upper)
  min(qgamma(runif(1) * Fu, shape, rate), upper)
}

# draw from the hyperprior by Gibbs on its gamma conditionals (used for the cold start)
r_ell_prior <- function(d, a, b, lambda, sweeps = 50) {
  ell <- rgamma(d, a, b)
  for (s in seq_len(sweeps)) for (j in seq_len(d))
    ell[j] <- rgamma(1, a, b + lambda * prod(1 + ell[-j]))
  ell
}

# ---------------------------------------------------------------- sampler
gp_pcn_mcmc <- function(binned, nd, z_eval,
                        n_iter = 20000, burn = 5000, thin = 10,
                        a_ell = 1, b_ell = 2.5, lambda = 0.1,
                        a_rho = 1, b_rho = 0.1, c_rho = 100,  # c_rho: generous ceiling; see README
                        ell_update = c("noncentred", "centred"),
                        n_w_steps = 10,         # pCN steps for w per Gibbs sweep (cheap; see README)
                        beta_init = 0.1, target_w = 0.30, target_ell = 0.44,
                        freeze_average = TRUE,  # freeze step sizes at their mean over 2nd half of burn-in
                        jitter = 1e-6, verbose = TRUE,
                        save_nodes_every = NULL,  # store (w, rho*) every k sweeps after burn-in
                        checkpoint_path = NULL, checkpoint_every = 10000L,
                        checkpoint_guard = NULL) {
  ell_update <- match.arg(ell_update)
  stopifnot(n_iter > burn, burn >= 0, thin >= 1, n_iter >= thin,
            n_iter == as.integer(n_iter), burn == as.integer(burn),
            thin == as.integer(thin), checkpoint_every >= 1,
            checkpoint_every == as.integer(checkpoint_every),
            is.null(save_nodes_every) || (save_nodes_every >= 1 &&
              save_nodes_every == as.integer(save_nodes_every)))
  # n_iter is a cumulative target, deliberately excluded from the continuation contract.
  contract <- NULL; saved <- NULL
  if (!is.null(checkpoint_path)) {
    if (!dir.exists(dirname(checkpoint_path))) stop("Checkpoint directory does not exist.")
    contract <- list(version = 1L, binned = binned, nd = nd, z_eval = z_eval,
      burn = burn, thin = thin, a_ell = a_ell, b_ell = b_ell, lambda = lambda,
      a_rho = a_rho, b_rho = b_rho, c_rho = c_rho, ell_update = ell_update,
      n_w_steps = n_w_steps, beta_init = beta_init, target_w = target_w,
      target_ell = target_ell, freeze_average = freeze_average, jitter = jitter,
      save_nodes_every = save_nodes_every,
      R_version = R.version.string, BLAS = sessionInfo()$BLAS,
      packages = lapply(c("Matrix", "coda"), function(p) as.character(packageVersion(p))),
      implementation = lapply(c("gp_pcn_mcmc", "interp_matrix", "make_factor", "ard_chol",
        "fmul", "fsolve", "flogdiag", "log_prior_ell", "rtgamma", "r_ell_prior",
        "mode_apply", "kron_apply", "bin_data", "make_nodes"),
        function(name) list(formals = formals(get(name)), body = body(get(name)))),
      external_guard = checkpoint_guard)
    if (file.exists(checkpoint_path)) {
      saved <- readRDS(checkpoint_path)
      if (!identical(saved$contract, contract))
        stop("Checkpoint model, data, source or environment changed; existing checkpoint preserved.")
      if (saved$iteration > n_iter) stop("Target is shorter than the saved chain.")
      if (verbose) cat("Resuming after sweep", saved$iteration, "towards", n_iter, "\n")
    }
  }
  d <- nd$d
  A_dom <- interp_matrix(nd, binned$z_dom); h_dom <- binned$h_dom
  A_obs <- interp_matrix(nd, binned$z_obs); h_obs <- binned$h_obs
  A_eval <- interp_matrix(nd, z_eval)
  D2 <- if (identical(nd$type, "grid")) NULL else
    lapply(seq_len(d), function(j) outer(nd$nodes[, j], nd$nodes[, j], "-")^2)
  V <- nrow(nd$nodes)
  upper <- c_rho + log(binned$n)

  loglik <- function(w, rs)
    sum(h_obs * plogis(as.vector(A_obs %*% w), log.p = TRUE)) -
      rs * sum(h_dom * plogis(as.vector(A_dom %*% w)))

  n_save <- n_iter %/% thin
  saved_iterations <- seq(thin, n_iter, by = thin)
  tr <- list(rho_star = numeric(n_save), ell = matrix(NA, n_save, d), loglik = numeric(n_save))
  draws <- matrix(NA_real_, sum(saved_iterations > burn), nrow(z_eval))
  acc_w <- numeric(n_iter); acc_ell <- matrix(FALSE, n_iter, d)
  n_node <- if (is.null(save_nodes_every)) 0L else max(0L, (n_iter - burn) %/% save_nodes_every)
  node_w <- matrix(NA_real_, n_node, V); node_rs <- numeric(n_node)
  start <- 0L; elapsed_before <- 0
  if (is.null(saved)) {
    # The original cold start and random-number ordering are unchanged.
    ell <- r_ell_prior(d, a_ell, b_ell, lambda)
    L <- make_factor(nd, D2, ell, jitter)
    w <- fmul(L, rnorm(V))
    rs <- rtgamma(a_rho, b_rho, upper)
    ll <- loglik(w, rs)
    beta <- beta_init; s_ell <- rep(0.5, d)
    lb_sum <- 0; ls_sum <- rep(0, d); n_avg <- 0
    s <- 0L; q <- 0L; qn <- 0L
  } else {
    start <- saved$iteration; elapsed_before <- saved$elapsed
    for (name in c("ell", "L", "w", "rs", "ll", "beta", "s_ell", "lb_sum", "ls_sum", "n_avg", "s", "q", "qn"))
      assign(name, saved$state[[name]])
    stopifnot(s == start %/% thin, q == sum(saved_iterations[saved_iterations <= start] > burn),
      length(saved$acc_w) == start, nrow(saved$acc_ell) == start,
      nrow(saved$draws) == q, ncol(saved$draws) == nrow(z_eval),
      all(is.finite(c(w, ell, rs, ll, beta, s_ell))), nrow(saved$node_w) == qn)
    acc_w[seq_len(start)] <- saved$acc_w
    acc_ell[seq_len(start), ] <- saved$acc_ell
    tr$rho_star[seq_len(s)] <- saved$trace$rho_star
    tr$ell[seq_len(s), ] <- saved$trace$ell
    tr$loglik[seq_len(s)] <- saved$trace$loglik
    draws[seq_len(q), ] <- saved$draws
    node_w[seq_len(qn), ] <- saved$node_w; node_rs[seq_len(qn)] <- saved$node_rs
    do.call(RNGkind, as.list(saved$rng_kind))
    assign(".Random.seed", saved$rng_seed, envir = .GlobalEnv)
    rm(saved)
  }
  t0 <- Sys.time()
  save_checkpoint <- function(iteration) {
    if (is.null(checkpoint_path)) return(invisible(NULL))
    state <- list(ell = ell, L = L, w = w, rs = rs, ll = ll, beta = beta,
      s_ell = s_ell, lb_sum = lb_sum, ls_sum = ls_sum, n_avg = n_avg, s = s, q = q, qn = qn)
    record <- list(contract = contract, iteration = iteration, state = state,
      rng_kind = RNGkind(), rng_seed = get(".Random.seed", envir = .GlobalEnv),
      elapsed = elapsed_before + as.numeric(difftime(Sys.time(), t0, units = "secs")),
      trace = list(rho_star = tr$rho_star[seq_len(s)], ell = tr$ell[seq_len(s), , drop = FALSE],
                   loglik = tr$loglik[seq_len(s)]),
      draws = draws[seq_len(q), , drop = FALSE], acc_w = acc_w[seq_len(iteration)],
      acc_ell = acc_ell[seq_len(iteration), , drop = FALSE],
      node_w = node_w[seq_len(qn), , drop = FALSE], node_rs = node_rs[seq_len(qn)])
    pending <- paste0(checkpoint_path, ".pending_", Sys.getpid())
    on.exit(unlink(pending))
    saveRDS(record, pending, compress = FALSE)
    if (!file.rename(pending, checkpoint_path)) stop("Could not save checkpoint atomically.")
    if (verbose) cat("\nCheckpoint saved after sweep", iteration, "\n")
  }
  for (it in if (start < n_iter) seq.int(start + 1L, n_iter) else integer()) {
    # 1. pCN for w (n_w_steps steps per sweep)
    for (k in seq_len(n_w_steps)) {
      w_p <- sqrt(1 - beta^2) * w + beta * fmul(L, rnorm(V))
      ll_p <- loglik(w_p, rs)
      if (log(runif(1)) < ll_p - ll) { w <- w_p; ll <- ll_p; acc_w[it] <- acc_w[it] + 1 / n_w_steps }
    }
    # 2. rho* | w
    rs <- rtgamma(a_rho + binned$N, b_rho + sum(h_dom * plogis(as.vector(A_dom %*% w))), upper)
    ll <- loglik(w, rs)
    # 3. ell_j, sequential
    for (j in seq_len(d)) {
      ell_p <- ell; ell_p[j] <- ell[j] * exp(s_ell[j] * rnorm(1))
      L_p <- tryCatch(make_factor(nd, D2, ell_p, jitter), error = function(e) NULL)
      if (is.null(L_p)) next
      lr <- log_prior_ell(ell_p, a_ell, b_ell, lambda) - log_prior_ell(ell, a_ell, b_ell, lambda) +
        log(ell_p[j]) - log(ell[j])                       # log-scale random walk Jacobian
      xi <- fsolve(L, w)
      if (ell_update == "centred") {
        xi_p <- fsolve(L_p, w)
        lr <- lr - flogdiag(L_p) - 0.5 * sum(xi_p^2) + flogdiag(L) + 0.5 * sum(xi^2)
        if (log(runif(1)) < lr) { ell <- ell_p; L <- L_p; acc_ell[it, j] <- TRUE }
      } else {
        w_p <- fmul(L_p, xi); ll_p <- loglik(w_p, rs)
        if (log(runif(1)) < lr + ll_p - ll) {
          ell <- ell_p; L <- L_p; w <- w_p; ll <- ll_p; acc_ell[it, j] <- TRUE
        }
      }
    }
    # adaptation, burn-in only (diminishing Robbins-Monro on the log scale)
    if (it <= burn) {
      g <- it^(-0.6)
      beta <- min(max(exp(log(beta) + g * (acc_w[it] - target_w)), 1e-3), 0.99)
      s_ell <- pmin(pmax(exp(log(s_ell) + g * (acc_ell[it, ] - target_ell)), 1e-3), 5)
      if (freeze_average && it > burn %/% 2) {
        lb_sum <- lb_sum + log(beta); ls_sum <- ls_sum + log(s_ell); n_avg <- n_avg + 1
        if (it == burn) { beta <- exp(lb_sum / n_avg); s_ell <- exp(ls_sum / n_avg) }
      }
    }
    if (it %% thin == 0) {
      s <- s + 1
      tr$rho_star[s] <- rs; tr$ell[s, ] <- ell; tr$loglik[s] <- ll
      if (it > burn) { q <- q + 1; draws[q, ] <- rs * plogis(as.vector(A_eval %*% w)) }
    }
    if (!is.null(save_nodes_every) && it > burn && (it - burn) %% save_nodes_every == 0) {
      qn <- qn + 1L; node_w[qn, ] <- w; node_rs[qn] <- rs
    }
    if (!is.null(checkpoint_path) && (it %% checkpoint_every == 0 || it == n_iter))
      save_checkpoint(it)
    if (verbose && it %% 1000 == 0)
      cat(sprintf("\rit %d | rho* %.2f | ell %s | acc_w %.2f | beta %.3f",
                  it, rs, paste(sprintf("%.2f", ell), collapse = ","),
                  mean(acc_w[max(1, it - 999):it]), beta))
  }
  elapsed <- elapsed_before + as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) cat("\n")
  keep <- seq(thin, n_iter, by = thin) > burn
  list(post_mean = colMeans(draws),
       lower = apply(draws, 2, quantile, 0.025), upper = apply(draws, 2, quantile, 0.975),
       draws = draws, trace = tr, keep = keep,
       acc_w = mean(acc_w[seq.int(burn + 1L, n_iter)]), acc_ell = colMeans(acc_ell[seq.int(burn + 1L, n_iter), , drop = FALSE]),
       beta = beta, s_ell = s_ell, time = elapsed,
       node_w = node_w[seq_len(qn), , drop = FALSE], node_rs = node_rs[seq_len(qn)],
       node_iterations = burn + save_nodes_every * seq_len(qn),
       ess = c(rho_star = unname(coda::effectiveSize(tr$rho_star[keep])),
               setNames(coda::effectiveSize(tr$ell[keep, , drop = FALSE]), paste0("ell", seq_len(d))),
               loglik = unname(coda::effectiveSize(tr$loglik[keep]))),
       ess_rho_eval = coda::effectiveSize(draws[, unique(round(seq(1, ncol(draws), length.out = 10))),
                                                drop = FALSE]))
}
