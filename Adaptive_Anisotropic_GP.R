library(pracma)
library(cubature)
library(rmutil)
library(fields)
library(mvnfast)
library(purrr)
library(furrr)
library(pals)
library(tictoc)
library(future)
library(future.apply)

# kernels with distance entry, d is a list of (1,2 or 3) pairwise distance matrices (one along each dimension) 
# length_scale is a vector of the same length of the list
# squared exponential
exponential_cov <- function(d, length_scale, sigma_2 = 1, nu = 0) {
  d <- Reduce(`+`, Map(function(di, li) di / li, d, length_scale))
  sigma_2 * exp(-d / 2)
}
#general matern
matern_kernel <- function(d1, d2, d3, length_scale, sigma_2, nu){
  d <-sqrt(Reduce(`+`, Map(function(di, li) di / li, d, length_scale)))
  d[d == 0] <- 1e-10
  b <- besselK(sqrt(2 * nu)*d, nu)
  c <- (sqrt(2 * nu)*d)^nu
  a <- 2^(1-nu) / gamma(nu)
  return(sigma_2 * a * c * b)
}
#matern 32
matern32_kernel <- function(d1, d2, d3, length_scale, sigma_2, nu){
  d <- sqrt(Reduce(`+`, Map(function(di, li) di / li, d, length_scale)))
  b <- 1+sqrt(3)*d
  a <- exp(-sqrt(3)*d)
  return(sigma_2 * a * b)
}
#matern 52
matern52_kernel <- function(d1, d2, d3, length_scale, sigma_2, nu){
  d <- sqrt(Reduce(`+`, Map(function(di, li) di / li, d, length_scale)))
  b <- 1 + sqrt(5)*d + 5/3*d^2
  a <- exp(-sqrt(5)*d)
  return(sigma_2 * a * b)
}

# (log) likelihood ratio by pointwise 
lik_ratio <- function(index_dom, index_loc, gp_new, gp_old, link, area_per_pixels, lambda_star){
  
  if(identical(link, "exponential")){
    log_rho_n <- gp_new[1,unlist(index_loc)]
    log_rho_o <- gp_old[1,unlist(index_loc)]
    I_new <- lapply(index_dom, function(i) {
      sum(exp(gp_new[1, i])) * area_per_pixels})
    I_old <- lapply(index_dom, function(i) {
      sum(exp(gp_old[1, i])) * area_per_pixels})
  }else{
    log_rho_n <- log(lambda_star * 1/(1 + exp(-gp_new[1,unlist(index_loc)])))
    log_rho_o <- log(lambda_star * 1/(1 + exp(-gp_old[1,unlist(index_loc)])))
    #### Integral over the observation window of the intensity function 
    I_new <- lapply(index_dom, function(i) {
      sum(lambda_star * 1/(1 + exp(-gp_new[1,i]))) * area_per_pixels})
    I_old <- lapply(index_dom, function(i) {
      sum(lambda_star * 1/(1 + exp(-gp_old[1,i]))) * area_per_pixels})
  }
  
  return(min(1, exp(sum(log_rho_n - log_rho_o) - sum(unlist(I_new)) + sum(unlist(I_old)))))
}

# anisotropic propsal for length scale
ind_anisotropic_proposal <- function(d_list, gp_cur, length_scale, index,
                                     theta, sigma_2, nu, Chol, 
                                     mu, kernel, jitter, alpha1, alpha2, adapt_cov){
  
  length_scale_prop <- length_scale
  length_scale_prop[index] <- exp(log(length_scale[index]) + rnorm(1, 0, adapt_cov[index])) 
  
  sigma_prop2 <- kernel(d_list, length_scale_prop, sigma_2, nu) + jitter 
  Chol2 <- Matrix::chol(sigma_prop2)
  
  MH_acc <- ((2*alpha1/theta  - 1) * log(length_scale_prop[index])) -  
    alpha2 * (length_scale_prop[index] ^(2/theta)) -
    ((2 * alpha1/theta - 1) * log(length_scale[index] )) +  
    alpha2*(length_scale[index] ^(2/theta)) + 
    dmvn(gp_cur[1,], mu, sigma = Chol2, isChol = TRUE, log = TRUE) -
    dmvn(gp_cur[1,], mu, sigma = Chol, isChol = TRUE, log = TRUE) + 
    dlnorm(length_scale[index], log(length_scale_prop[index]), adapt_cov[index], log = TRUE) -
    dlnorm(length_scale_prop[index], log(length_scale[index]), adapt_cov[index], log = TRUE) 
  
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    length_scale_j <- length_scale_prop[index]
    Chol <- Chol2
  }else{
    length_scale_j <- length_scale[index]
  }
  return(list(length_scale_j, Chol))
}

# length scale exponents proposal
ind_theta_proposal <- function(theta, length_scale, exp_param, alpha1, alpha2){
  # browser()
  theta_prop <- rbeta(1, exp_param[1], exp_param[2])
  MH_acc <- (dbeta(theta, exp_param[1], exp_param[2], log = TRUE)) -         # proposal contribution 
    (dbeta(theta_prop, exp_param[1], exp_param[2], log = TRUE)) +
    (alpha1/theta_prop*2 * log(length_scale)) -                                # gamma likelihood kernel
    alpha2*(length_scale^(2/theta_prop)) -
    (alpha1/theta*2 * log(length_scale)) +  
    alpha2*(length_scale^(2/theta)) +                                           
    ((exp_param[1]-2)*log(theta_prop) + (exp_param[2]-1)*log(1 - theta_prop)) -   # prior beta kernel
    ((exp_param[1]-2)*log(theta) + (exp_param[2]-1)*log(1 - theta))
  
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    theta <- theta_prop 
  }
  return(theta)
}


MCMC <- function(run_name, niter, 
                 loc, cov_lists, 
                 discretization,  window, kernel, 
                 beta, alpha1, alpha2, sigma_2, nu, 
                 exp_param, shape, rate, link,
                 index_dom = NA, index_loc = NA){
  
  P <- length(cov_lists)
  K <- dim(discretization)[1]
  jitter <- 1e-10 * diag(1,K) 
  
  tic()
  # pairwise distances between points on the discretization grid
  d_list <- lapply(seq_len(ncol(discretization)), function(j) {
    outer(discretization[, j], discretization[, j], FUN = function(a, b) (a - b)^2)
  })
  # covariate values at the locations of observations
  plan(multisession, workers = availableCores())
  
  Z_list <- future_lapply(cov_lists, function(covariate_process, loc_arg){
    library(spatstat)
    Map(function(cov, lo){cov[lo]},covariate_process, loc_arg)
  }, loc_arg = loc)

  area <- area.owin(window)
  pixels <- length(cov_lists[[1]][[1]]$v) 
  area_per_pixels <- area/pixels
  browser()
  # indices of the discretization closest to the covariate values over the whole observation window
  index_dom <- future_pmap(cov_lists, function(...) {
    covariates <- list(...)
    mat <- do.call(cbind, lapply(covariates, function(cov) na.omit(as.vector(cov$v))))
    k <- ncol(mat)

    apply(mat, 1, function(w) {
      which.min(
        rowSums((discretization -
                   matrix(w, nrow = nrow(discretization), ncol = k, byrow = TRUE))^2)
      )
    })
  }, .options = furrr::furrr_options(packages = c("spatstat.geom"),
                                     seed = TRUE,
                                     globals = list(discretization = discretization))
  )
  
  plan(sequential)
  assign("Index_over_dom", index_dom, envir = .GlobalEnv)
  # indices of the discretization closest to the covariate values realized at event locations
  index_loc <- pmap(Z_list, function(...) {
    Zs <- list(...)

    mat <- do.call(cbind, Zs)
    k <- ncol(mat)

    apply(mat, 1, function(w) {
      which.min(
        rowSums((discretization -
                   matrix(w, nrow = nrow(discretization), ncol = k, byrow = TRUE))^2)
      )
    })
  })
  assign("Index_over_loc", index_loc, envir = .GlobalEnv)
  toc()
  
  # mean function 
  mu = rep(0, K)
  
  # output files 
  folder = paste("~/Research/IPP/", run_name, sep = "")
  if(!dir.exists(folder)){
    dir.create(folder)
  }
  setwd(folder)
  
  file.create(description = "Post_gp.csv")
  file.create(description = "Post_ls.csv")
  file.create(description = "Post_l_star.csv")
  file.create(description = "Post_theta.csv")										   
  
  theta <- rbeta(P, exp_param[1], exp_param[2])
  length_scale <- as.vector((rgamma(P, shape = alpha1, rate = alpha2))^(theta/P))
  
  while(any(length_scale>1 | length_scale<0.05)){
            theta <- rbeta(P, exp_param[1], exp_param[2])
            length_scale <- as.vector((rgamma(P, shape = alpha1, rate = alpha2))^(theta/P))
  }
  mean_ls <- length_scale
  cov_ls <- rep(0,P)
  
  Cov <- kernel(d_list, length_scale, sigma_2, nu) + jitter
  Chol <- Matrix::chol(Cov)
  gp_cur <- rnorm(n = K, mean = 0, sd = 1) %*% Chol + mu
  lambda_star <- rgamma(1, shape, rate)

  gp <- file(description = "Post_gp.csv", open = "a")
  ls <- file(description = "Post_ls.csv", open = "a")
  lstar <- file(description = "Post_l_star.csv", open = "a")
  dir_theta <- file(description = "Post_theta.csv", open = "a")

  N <- sum(sapply(loc, npoints))
  
  time <- Sys.time()
  
  for(i in 1:niter){
    # pCN proposal
    # a new Covariance matrix is needed only if the length scale changes, so postponed to latter update
    gp_prior <- rnorm(n = K, mean = 0, sd = 1) %*% Chol + mu
    # proposal (log scale)
    gp_prop <- sqrt(1 - beta^2) * gp_cur + beta * gp_prior
    # likelihood ratio
    lr <- lik_ratio(index_dom, index_loc, 
                    gp_new = gp_prop, gp_old = gp_cur, 
                    link, area_per_pixels, lambda_star)
    if(runif(1,0,1) < lr){
      gp_cur <- gp_prop 
    }
    if(identical(link, "exponential")){
      write.table(t(exp(gp_cur[1,])), file = gp, sep=',', append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }else{
      write.table(t(lambda_star * 1/(1 + exp(-gp_cur[1,]))), file = gp, sep=',', append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
      I <- lapply(index_dom, function(i) {			  
        sum(1/(1 + exp(-gp_cur[1, i]))) * area / pixels})
      lambda_star <- rgamma(1, shape = shape + N, scale = 1/(rate + sum(unlist(I))))
      write.table(lambda_star, file = lstar, sep=',', append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }									 
    
    # length_scale parameter update with independent prior
    if(i > niter/5){
      adapt_cov <- 2.4^2 * (cov_ls + rep(0.01,P)) 
    }else{
      adapt_cov <- rep(0.2,P)
    }
    
    # on parallel architecture, uncomment future_lapply and comment for loop
    {
    parallel_res <- future_lapply(1:P, function(j) {
      res <- ind_anisotropic_proposal(
        d_list, gp_cur, length_scale, j,
        theta[j], sigma_2, nu, Chol,
        mu, kernel, jitter, alpha1, alpha2, adapt_cov
      )

      length_scale_j <- res[[1]]
      Chol_j <- res[[2]]

      theta_j <- ind_theta_proposal(
        theta[j], length_scale_j, exp_param, alpha1, alpha2
      )
      list(l_j = length_scale_j, t_j = theta_j, C_j = Chol_j)
    }, future.seed=TRUE)

    theta <- vapply(parallel_res, function(x) x$t_j, numeric(1))
    length_scale_new <- vapply(parallel_res, function(x) x$l_j, numeric(1))

    check <- (length_scale_new != length_scale)
    if(sum(check) == 1){
      Chol <- parallel_res[[which(check)]]$C_j
    }else if(sum(check) == 0){

    }else{
      Chol <- Matrix::chol(kernel(d_list, length_scale_new, sigma_2, nu) + jitter) }
    }
    
    
    # {
    # length_scale_new <- rep(NA,P)
    # Chol_list <- list()
    # for(j in 1:P){
    #     res <- ind_anisotropic_proposal(
    #       d_list, gp_cur, length_scale, j,
    #       theta[j], sigma_2, nu, Chol,
    #       mu, kernel, jitter, alpha1, alpha2, adapt_cov
    #     )
    # 
    #     length_scale_j <- res[[1]]
    #     Chol_j <- res[[2]]
    # 
    #     theta[j] <- ind_theta_proposal(
    #       theta[j], length_scale_j, exp_param, alpha1, alpha2
    #     )
    #     length_scale_new[j] <- length_scale_j
    #     Chol_list[[j]] <- Chol_j
    #   }
    #   
    # # check which length scales updated and if Covaraince matrix has to be recomputed
    # check <- (length_scale_new != length_scale)
    # if(sum(check) == 1){
    #   Chol <- Chol_list[[which(check)]]
    # }else if(sum(check) == 0){
    #   
    # }else{
    #   Chol <- Matrix::chol(kernel(d_list, length_scale_new, sigma_2, nu) + jitter) 
    # }
    # }
    
    length_scale <- length_scale_new
    
    mean_old <- mean_ls
    mean_ls <- 1/(i+1) * (i * mean_old + length_scale)
    cov_ls <- 1/i * ((i-1) * cov_ls + length_scale^2 + 
                       i * mean_old^2 - (i+1) * mean_ls^2) 
    
    write.table(t(theta), file = dir_theta, sep=',', append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    write.table(t(length_scale), file = ls, sep=',', append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
    cat(paste("\rIteration", i))
  }
  print(difftime(Sys.time(), time, units = "secs"))
  close(lstar)
  close(ls)
  close(gp)
  close(dir_theta)
}
