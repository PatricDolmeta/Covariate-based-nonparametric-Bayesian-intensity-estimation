# Naive preconditioned Crank–Nicolson algorithm 
library(pracma)
library(cubature)
library(rmutil)
library(fields)
library(mvnfast)
# kernels with distance entry
# squared exponential
exponential_cov <- function(d1, d2, d3, length_scale = c(1,1,1), sigma_2 = 1, nu = 0) {
  d <- (d1/length_scale[1] + d2/length_scale[2] + d3/length_scale[3])
  sigma_2 * exp(-d / 2)
}
#general matern
matern_kernel <- function(d1, d2, d3, length_scale, sigma_2, nu){
  d <-sqrt(d1/length_scale[1] + d2/length_scale[2] + d3/length_scale[3])
  d[d == 0] <- 1e-10
  b <- besselK(sqrt(2 * nu)*d, nu)
  c <- (sqrt(2 * nu)*d)^nu
  a <- 2^(1-nu) / gamma(nu)
  return(sigma_2 * a * c * b)
}

#matern 32
matern32_kernel <- function(d1, d2, d3, length_scale, sigma_2, nu){
  d <- sqrt(d1/length_scale[1] + d2/length_scale[2] + d3/length_scale[3])
  b <- 1+sqrt(3)*d
  a <- exp(-sqrt(3)*d)
  return(sigma_2 * a * b)
}

#matern 52
matern52_kernel <- function(d1, d2, d3, length_scale, sigma_2, nu){
  d <- sqrt(d1/length_scale[1] + d2/length_scale[2] + d3/length_scale[3])
  b <- 1 + sqrt(5)*d + 5/3*d^2
  a <- exp(-sqrt(5)*d)
  return(sigma_2 * a * b)
}

# likelihood ratio for IPP
lik_ratio <- function(vertices, index_dom, index_loc, K, gp_new, gp_old, link, area, pixels, lambda_star){
  # browser()
  if(identical(link, "exponential")){
    log_rho_n <- gp_new[1,unlist(index_loc)]
    log_rho_o <- gp_old[1,unlist(index_loc)]
    I_new <- lapply(index_dom, function(i) {
       sum(exp(gp_new[1, i])) * area / pixels})
    I_old <- lapply(index_dom, function(i) {
       sum(exp(gp_old[1, i])) * area / pixels})
  }else{
  log_rho_n <- log(lambda_star * 1/(1 + exp(-gp_new[1,unlist(index_loc)])))
  log_rho_o <- log(lambda_star * 1/(1 + exp(-gp_old[1,unlist(index_loc)])))
  #### Integral over the observation window of the intensity function 
  I_new <- lapply(index_dom, function(i) {
       sum(lambda_star * 1/(1 + exp(-gp_new[1,i]))) * area / pixels})
  I_old <- lapply(index_dom, function(i) {
		sum(lambda_star * 1/(1 + exp(-gp_old[1,i]))) * area / pixels})
  }
  return(min(1, exp(sum(log_rho_n - log_rho_o) - sum(unlist(I_new)) + sum(unlist(I_old)))))
}

MCMC <- function(run_name, niter, beta, K, loc, window,
				 covariate_process1, covariate_process2, covariate_process3, 
				 vertex, kernel, alpha1, alpha2, sigma_2, nu, 
				 exp_param, shape, rate, link, LS,
				 index_dom = NA, index_loc = NA){
  
  # browser()
  
  jitter <- 1e-10 * diag(1,K) 
  
  # pairwise distances between points on the grid
  d1 <- outer(vertex[,1], vertex[,1], FUN = function(a, b) (a - b)^2)
  d2 <- outer(vertex[,2], vertex[,2], FUN = function(a, b) (a - b)^2)
  d3 <- outer(vertex[,3], vertex[,3], FUN = function(a, b) (a - b)^2)
  # covariate values at the locations of observations
  # browser()
  Z1 <- Map(function(cov, lo){cov[lo]}, covariate_process1, loc)
  Z2 <- Map(function(cov, lo){cov[lo]}, covariate_process2, loc)
  Z3 <- Map(function(cov, lo){cov[lo]}, covariate_process3, loc)
  
  # extremes of integration 
  # xl <- covariate_process1$xrange[1]
  # yl <- covariate_process1$yrange[1]
  # xu <- covariate_process1$xrange[2]
  # yu <- covariate_process1$yrange[2]
  # area <- (yu-yl)*(xu-xl)
  area <- area.owin(window)
  pixels <- length(covariate_process2[[1]]$v) #- sum(is.na(covariate_process2[[1]]$v))
  
  
  # covariate values on the obsgrid
  #x <- lapply(covariate_process1, function(c1) sapply(c1$v, function(w) d_k$x[which.min(abs(d_k$x - w))]))
  #y <- lapply(covariate_process2, function(c2) sapply(c2$v, function(w) d_k$y[which.min(abs(d_k$y - w))]))
  #x_indexes <- lapply(x, function(xi) {match(xi, discr$x)})
  #y_indexes <- lapply(y, function(yi) {match(yi, discr$y)})
  #indexes <- Map(cbind, y_indexes, x_indexes)
  
   # Precomputed
    time <- Sys.time()
    index_dom <- Map(function(cov1, cov2, cov3){unlist(apply(
      cbind(na.omit(as.vector(cov1$v)),na.omit(as.vector(cov2$v)),na.omit(as.vector(cov3$v))),
                                    1, function(w){which.min(rowSums((vertex - matrix(w, nrow = dim(vertex)[1], ncol = 3, byrow = TRUE))^2))}))},
  					              covariate_process1, covariate_process2, covariate_process3)
    print(difftime(Sys.time(), time, "secs"))
    assign("Index3d_over_dom_prec", index_dom, envir = .GlobalEnv)
    index_loc <- Map(function(z1, z2, z3){apply(cbind(z1,z2,z3), 1,
                       function(w){which.min(rowSums((vertex - matrix(w, nrow = dim(vertex)[1], ncol = 3, byrow = TRUE))^2))})},
  					           Z1, Z2, Z3)
    assign("Index3d_over_loc_prec", index_loc, envir = .GlobalEnv)

  # browser()
  
  # discretization grid for the GP prior on rho
  # mean function 
  mu = rep(0, K)
  
  # output files 
  folder = paste("~/Research/IPP/","3d_mesh", run_name, sep = "")
  if(!dir.exists(folder)){
    dir.create(folder)
  }
  setwd(folder)
  
  if(file.exists("Post_gp.csv")) {
    browser()
    # Read existing data
    lambda_post = read.table("Post_l_star.csv", sep = ",")
	  expon = read.table("Post_theta.csv", sep = ",")
    											   
    Ls = read.table("Post_ls.csv", sep = ",")
    gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
    gs <- apply(gs, 2, as.numeric)
    # Extract the last saved value
    lambda_star <- as.numeric(tail(lambda_post, 1)) 
	  theta <- as.numeric(tail(expon, 1)) 
    length_scale <- as.numeric(tail(Ls, 1)) 
    
    mean_ls <- colMeans(Ls) 
    cov_ls <- c(var(Ls[,1]), var(Ls[,2]))
    gp_cur <- t(matrix(tail(gs,1)))
    # mean function 
    mu = rep(0, K)
    # covariance matrix of the GP generating the covariate process, discretizd over the grid
    Cov <- kernel(d1, d2, d3, length_scale, sigma_2, nu) + jitter
    Chol <- Matrix::chol(Cov)
  } else {
    file.create(description = "Post_gp.csv")
    file.create(description = "Post_ls.csv")
    file.create(description = "Post_l_star.csv")
    file.create(description = "Post_theta.csv")										   
    # mean function 
    mu = rep(0, K)
	if(identical(LS, "isotropic")){
	  length_scale <- rgamma(1, shape = alpha1, rate = alpha2)
	  length_scale <- rep(length_scale, 3)
	}else{
	  # theta <- rdirichlet(1, exp_param)
	  theta <- rbeta(3, exp_param[1], exp_param[2])
	  length_scale <- as.vector((rgamma(3, shape = alpha1, rate = alpha2))^(theta/3))
	  length_scale
	  
	  mean_ls <- length_scale
	  cov_ls <- rep(0,3)
	  
	  while(length_scale[1]>1 | length_scale[2]>1 | length_scale[3]>1 | 
	        length_scale[1]<0.05 | length_scale[2]<0.05 | length_scale[3]<0.05){
	    theta <- rbeta(3, exp_param[1], exp_param[2])
	    length_scale <- as.vector((rgamma(3, shape = alpha1, rate = alpha2))^(theta/3))
	    length_scale
	    
	    mean_ls <- length_scale
	    cov_ls <- rep(0,3)
	  }
	  
	}	
    # covariance matrix of the GP generating the covariate process, discretizd over the grid
    # browser()
    Cov <- kernel(d1, d2, d3, length_scale, sigma_2, nu) + jitter
    Chol <- Matrix::chol(Cov)
    # log-GP prior on rho
    gp_cur <- rnorm(n = K, mean = 0, sd = 1) %*% Chol + mu

    # length_scale <- rexp(1, rate = 1)
    lambda_star <- rgamma(1, shape, rate)
  }
  gp <- file(description = "Post_gp.csv", open = "a")
  ls <- file(description = "Post_ls.csv", open = "a")
  lstar <- file(description = "Post_l_star.csv", open = "a")
  dir_theta <- file(description = "Post_theta.csv", open = "a")
  
  time <- Sys.time()
  Acc <- 1	 

  N <- sum(sapply(loc, npoints))
  
  # browser()
  
  for(i in 1:niter){
    # pCN proposal
    # a new Covariance matrix is needed only if the length scale changes, so postponed to latter update
    # sample form the prior 
    # browser()
    gp_prior <- rnorm(n = K, mean = 0, sd = 1) %*% Chol + mu
    # proposal (log scale)
    gp_prop <- sqrt(1 - beta^2) * gp_cur + beta * gp_prior
    lr <- lik_ratio(vertex, index_dom, index_loc,
					K, gp_new = gp_prop, gp_old = gp_cur, link,
					area, pixels, lambda_star)
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
    
    # length_scale parameter update
    # browser()
	if(identical(LS, "isotropic")){
      res <- ls_isotropic_proposal(Acc, d1, d2, d3, gp_cur, length_scale, 
                                   sigma_2, nu, K, Chol, Cov, mu, kernel, 
                                   jitter, alpha1, alpha2)
      length_scale <- res[[1]]
      Cov <- res[[2]]
      Chol <- res[[3]]
    }else{
      # { # with dirichelet hierarchy
      # res <- ls_anisotropic_proposal(Acc, d1, d2, gp_cur, length_scale, theta,
      #                                sigma_2, nu, K, Chol, Cov, mu, kernel,
      #                                jitter, alpha1, alpha2)
      # length_scale <- res[[1]]
      # Cov <- res[[2]]
      # Chol <- res[[3]]
      # # browser()
      # theta <- theta_proposal(theta, length_scale, exp_param, alpha1, alpha2)
      # }
      { # with independent priors
        if(i > niter/5){
          adapt_cov <- 2.4^2 * (cov_ls + rep(0.01,3)) 
        }else{
          adapt_cov <- rep(0.2,3)
        }
        
        res <- ind_anisotropic_proposal(Acc, d1, d2, d3, gp_cur, length_scale,
                                        theta, sigma_2, nu, K, Chol, Cov,
                                        mu, kernel, jitter, alpha1, alpha2, adapt_cov)
        length_scale <- res[[1]]
        Cov <- res[[2]]
        Chol <- res[[3]]
        # Acc <- res[[4]]
        
        # sequential update of length scale proposal 
        # mean 
        # browser()
        mean_old <- mean_ls
        mean_ls <- 1/(i+1) * (i * mean_old + length_scale)
        # covariance 
        cov_ls <- 1/i * ((i-1) * cov_ls + length_scale^2 + 
                              i * mean_old^2 - (i+1) * mean_ls^2) 
        
        theta <- ind_theta_proposal(theta, length_scale, exp_param, alpha1, alpha2)
		 
      }
    
	write.table(t(theta), file = dir_theta, sep=',', append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }																				
    write.table(t(length_scale), file = ls, sep=',', append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
																			
    cat(paste("\rIteration", i, Acc))
  }
  print(difftime(Sys.time(), time, units = "secs"))
  close(lstar)
  close(ls)
  close(gp)
  close(dir_theta)
}


########################### ISOTROPIC ###############################
ls_isotropic_proposal <- function(Acc, d1, d2, d3, gp_cur, length_scale, 
                                  sigma_2, nu, K, Chol, Cov, mu, kernel, 
                                  jitter, alpha1, alpha2){
  # browser()
  length_scale_prop <- rep(exp(log(length_scale[1]) + rnorm(1, 0, 0.2)),3)
  sigma_prop <- kernel(d1, d2, d3, length_scale_prop, sigma_2, nu) + jitter 
  # kernel(distance_matrix, length_scale_prop, sigma_2, nu) + jitter
  sigma_old <- Cov
  MH_acc <- dgamma(length_scale_prop[1], shape = alpha1, rate = alpha2, log = TRUE) -
    dgamma(length_scale[1], shape = alpha1, rate = alpha2, log = TRUE) + 
    dmvn(gp_cur[1,], mu, sigma = sigma_prop, isChol = FALSE, log = TRUE) -
    dmvn(gp_cur[1,], mu, sigma = Chol, isChol = TRUE, log = TRUE) + 
    dlnorm(length_scale[1], log(length_scale_prop[1]), 0.2, log = TRUE) -
    dlnorm(length_scale_prop[1], log(length_scale[1]), 0.2, log = TRUE) 
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    length_scale <- length_scale_prop
    # new covariance matrix 
    Cov <- kernel(d1, d2, d3, length_scale, sigma_2, nu) + jitter 
    # kernel(distance_matrix, length_scale, sigma_2, nu) + jitter
    Chol <- Matrix::chol(Cov)
    # Acc <- Acc + 1
  }
  return(list(length_scale, Cov, Chol, Acc))
}

############################ ANISOTROPIC ###############################
ls_anisotropic_proposal <- function(Acc, d1, d2, d3, gp_cur, length_scale, 
                                    theta, sigma_2, nu, K, Chol, Cov, 
                                    mu, kernel, jitter, alpha1, alpha2){
  # browser()
  length_scale_prop <- rep(NA,3)
  length_scale_prop[1] <- exp(log(length_scale[1]) + rnorm(1, 0, 0.2)) 
  length_scale_prop[2] <- exp(log(length_scale[2]) + rnorm(1, 0, 0.2)) 
  length_scale_prop[3] <- exp(log(length_scale[3]) + rnorm(1, 0, 0.2)) 

  sigma_prop2 <- kernel(d1, d2, d3, length_scale_prop, sigma_2, nu) + jitter 
  # kernel(distance_matrix, length_scale_prop, sigma_2, nu) + jitter
  sigma_old <- Cov
  MH_acc <- sum((alpha1/theta - 1) * log(length_scale_prop)) -  # posterior kernel
    alpha2 * sum(length_scale_prop^(1/theta)) -
    sum((alpha1/theta-1) * log(length_scale)) +  
    alpha2*sum(length_scale^(1/theta)) + 
    dmvn(gp_cur[1,], mu, sigma = sigma_prop2, isChol = FALSE, log = TRUE) -
    dmvn(gp_cur[1,], mu, sigma = Chol, isChol = TRUE, log = TRUE) + 
    dlnorm(length_scale[1], log(length_scale_prop[1]), 0.2, log = TRUE) -
    dlnorm(length_scale_prop[1], log(length_scale[1]), 0.2, log = TRUE) +
    dlnorm(length_scale[2], log(length_scale_prop[2]), 0.2, log = TRUE) -
    dlnorm(length_scale_prop[2], log(length_scale[2]), 0.2, log = TRUE) +
	dlnorm(length_scale[3], log(length_scale_prop[3]), 0.2, log = TRUE) -
    dlnorm(length_scale_prop[3], log(length_scale[3]), 0.2, log = TRUE)
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    length_scale <- length_scale_prop 
    # new covariance matrix 
    Cov <- kernel(d1, d2, length_scale, sigma_2, nu) + jitter
    # kernel(distance_matrix, length_scale, sigma_2, nu) + jitter
    Chol <- Matrix::chol(Cov)
    # Acc <- Acc+1
  }
  return(list(length_scale, Cov, Chol, Acc))
}

############################ THETA DIRICHLET ###############################
theta_proposal <- function(theta, length_scale, exp_param, alpha1, alpha2){
  # browser()
  theta_prop <- rdirichlet(1, exp_param)
  MH_acc <- ddirichlet(theta, exp_param, log = TRUE) -  # proposal contribution 
    ddirichlet(theta_prop, exp_param, log = TRUE) + 
    sum(alpha1/theta_prop * log(length_scale)) -  # posterior kernel
    alpha2*sum(length_scale^(1/theta_prop))+
    sum((exp_param - 2) * log(theta_prop)) -
    sum(alpha1/theta * log(length_scale)) +  
    alpha2*sum(length_scale^(1/theta))-
    sum((exp_param - 2) * log(theta))
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    theta <- theta_prop 
  }
  return(theta)
}

############################ IND ANISOTROPIC ###############################
ind_anisotropic_proposal <- function(Acc, d1, d2, d3, gp_cur, length_scale,
                                    theta, sigma_2, nu, K, Chol, Cov,
                                    mu, kernel, jitter, alpha1, alpha2, adapt_cov){
  
  length_scale_prop <- rep(NA,3)
  length_scale_prop[1] <- exp(log(length_scale[1]) + rnorm(1, 0, adapt_cov[1])) 
  length_scale_prop[2] <- exp(log(length_scale[2]) + rnorm(1, 0, adapt_cov[2])) 
  length_scale_prop[3] <- exp(log(length_scale[3]) + rnorm(1, 0, adapt_cov[3])) 
  
  sigma_prop2 <- kernel(d1, d2, d3, length_scale_prop, sigma_2, nu) + jitter 
  sigma_old <- Cov
			 
  
  MH_acc <- sum(theta/3 * (alpha1 - 1) * log(length_scale_prop)) -  
    alpha2 * sum(length_scale_prop^(theta/3)) -
    sum(theta/3 * (alpha1 - 1) * log(length_scale)) +  
    alpha2*sum(length_scale^(theta/3)) + 
    dmvn(gp_cur[1,], mu, sigma = sigma_prop2, isChol = FALSE, log = TRUE) -
    dmvn(gp_cur[1,], mu, sigma = Chol, isChol = TRUE, log = TRUE) + 
    dlnorm(length_scale[1], log(length_scale_prop[1]), adapt_cov[1], log = TRUE) -
    dlnorm(length_scale_prop[1], log(length_scale[1]), adapt_cov[1], log = TRUE) +
    dlnorm(length_scale[2], log(length_scale_prop[2]), adapt_cov[2], log = TRUE) -
    dlnorm(length_scale_prop[2], log(length_scale[2]), adapt_cov[2], log = TRUE) +
	dlnorm(length_scale[3], log(length_scale_prop[3]), adapt_cov[3], log = TRUE) -
    dlnorm(length_scale_prop[3], log(length_scale[3]), adapt_cov[3], log = TRUE)
  
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    length_scale <- length_scale_prop
    # new covariance matrix
    Cov <- kernel(d1, d2, d3, length_scale, sigma_2, nu) + jitter
																 
    Chol <- Matrix::chol(Cov)
    # Acc <- Acc+1
  }
  return(list(length_scale, Cov, Chol, Acc))
}
############################# IND EXPONENTS ############################

ind_theta_proposal <- function(theta, length_scale, exp_param, alpha1, alpha2){
  # browser()
  theta_prop <- rbeta(3, exp_param[1], exp_param[2])
  MH_acc <- sum(dbeta(theta, exp_param[1], exp_param[2], log = TRUE)) -         # proposal contribution 
    sum(dbeta(theta_prop, exp_param[1], exp_param[2], log = TRUE)) +
    sum(alpha1*theta_prop/3 * log(length_scale)) -                                # gamma likelihood kernel
    alpha2*sum(length_scale^(theta_prop/3)) -
    sum(alpha1*theta/3 * log(length_scale)) +  
    alpha2*sum(length_scale^(theta/3)) +                                        # change of measure
    sum(log(theta_prop) - log(theta)) +                                             
    sum((exp_param[1]-1)*log(theta_prop) + (exp_param[2]-1)*log(1 - theta_prop)) -   # prior beta kernel
    sum((exp_param[1]-1)*log(theta) + (exp_param[2]-1)*log(1 - theta))
  
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    theta <- theta_prop 
  }
  return(theta)
}