# Naive preconditioned Crank–Nicolson algorithm 
library(pracma)
library(cubature)
library(rmutil)
library(mvnfast)
# square exponential kernel with scalar entries 
squared_exponential_kernel <- function(x, y, length_scale, sigma_2, nu){
  return(sigma_2 * exp(-(x-y)^2/(2 * length_scale) ))
}
matern_kernel <- function(x, y, length_scale, sigma_2, nu){
  d <-sqrt((x-y)^2)  
  d[d == 0] <- 1e-10				
  b <- besselK(sqrt(2 * nu)*d/length_scale, nu)
  c <- (sqrt(2 * nu)*d/length_scale)^nu
  a <- 2^(1-nu) / gamma(nu)
  return(sigma_2 * a * c * b)
}
matern32_kernel <- function(x, y, length_scale, sigma_2, nu){
  d <-sqrt((x-y)^2) 
  b <- 1+sqrt(3)*d/length_scale
  a <- exp(-sqrt(3)*d/length_scale)
  return(sigma_2 * a * b)
}
matern52_kernel <- function(x, y, length_scale, sigma_2, nu){
  d <-sqrt((x-y)^2) 
  b <- 1+sqrt(5)*d/length_scale+5*d^2/(3*length_scale^2)
  a <- exp(-sqrt(5)*d/length_scale)
  return(sigma_2 * a * b)
}

# functions returning the value of the covariate process at given location (x,y)
# by bilinear interpolation
covariate_function2 <- function(x, y, grid, covar){
  index <- which.min(sqrt((grid$x - x)^2 + ((grid$y) - y)^2))
  return(interp.im(covar, x = grid$x[index], y = grid$y[index]))
}
# by nearest neighbour in the grid
covariate_function <- function(x, y, covariates, grid){
  index <- which.min(sqrt((grid$x - x)^2 + ((grid$y) - y)^2))
  X_index <- which.min(abs(covariates$xcol - grid$x[index]))
  Y_index <- which.min(abs(covariates$yrow - grid$y[index]))
  return(covariates[Y_index, X_index])
}
# for given discretization of the intensity function rho over the grid d_k
# return the value of the intenity function at given coordinates (x,y) by linear interpolation
f_intensity <- function(x, y, rho, d, xl, yl, xu, yu, grid, covar){
  return(approx(d, rho, covariate_function(xl + x * (xu-xl), yl + y * (yu-yl), covariates = covar, grid = grid))$y)
}
# likelihood ratio for IPP
lik_ratio <- function(d_k, K, rho_new, rho_old, z_loc, grid, covar, lambda_star, area, index_loc, index_dom){
  #value of the proposed and old intensity functions at the locations of observations
  # rho_n <- approx(d_k, rho_new, xout = unlist(z_loc))$y
  # rho_o <- approx(d_k, rho_old, xout = unlist(z_loc))$y
  # browser()
  rho_n <- log(lambda_star * 1/(1 + exp(-rho_new[1,unlist(index_loc)])))
  rho_o <- log(lambda_star * 1/(1 + exp(-rho_old[1,unlist(index_loc)])))
  # rho_n <- log(lambda_star * 1/(1 + exp(-approx(d_k, rho_new, xout = unlist(z_loc))$y)))
  # rho_o <- log(lambda_star * 1/(1 + exp(-approx(d_k, rho_old, xout = unlist(z_loc))$y)))
  # #integral over the observation window of the rho_new function
  # I_new <- lapply(covar, function(c){sum(exp(rho_new)[index_dom])*area/length(grid$x)})
  I_new <- lapply(index_dom, function(i){sum(lambda_star * 1/(1 + exp(-rho_new)[i]), na.rm = TRUE)*area/length(grid$x)})
  # precise counterpart 
  # lapply(covar, function(CC){(yu-yl) * (xu-xl) * integrate2d(f_intensity, rho = exp(rho_new), d = d_k,
  #                                        xl = xl, yl = yl, xu = xu, yu = yu, grid = grid, covar = CC, error = 1e-05)$value})
  # I_old <- lapply(covar, function(c){sum(exp(rho_old)[index_dom])*area/length(grid$x)})
  I_old <- lapply(index_dom, function(i){sum(lambda_star * 1/(1 + exp(-rho_old)[i]), na.rm = TRUE)*area/length(grid$x)})
  # lapply(covar, function(CC){(yu-yl) * (xu-xl) * integrate2d(f_intensity, rho = exp(rho_old), d = d_k, 
  #                                                                     xl = xl, yl = yl, xu = xu, yu = yu, grid = grid, covar = CC, error = 1e-05)$value})
  return(min(1, exp(sum(rho_n - rho_o) - sum(unlist(I_new)) + sum(unlist(I_old)))))
}

adapt_multi_MCMC <- function(folder_name, niter, beta, K, loc, grid, covariate_process, d_k, 
                       kernel, alpha1, alpha2, sigma_2, nu, 
						           exp_param, shape, rate, link){
	# output files 
  # browser()

  folder = paste("Adaptive_MCMC_post_multi/", folder_name, sep = "")
  dir.create(paste("~/Research/IPP/", folder, sep = ""))
  setwd(folder)
  file.create(description = "Post_gp.csv")
  file.create(description = "Post_ls.csv")
  file.create(description = "Post_l_star.csv")
  file.create(description = "Post_theta.csv")
  
  length_scale <- 0
  
  while(length_scale>1 | length_scale<0.1){
    theta <- rbeta(1, exp_param[1], exp_param[2])
    length_scale <- (rgamma(1, shape = alpha1, rate = alpha2))^(theta)
    length_scale
    
    mean_ls <- length_scale
    cov_ls <- 0
  }
  
  gp	<- file(description = "Post_gp.csv", open = "a")
  ls <- file(description = "Post_ls.csv", open = "a")
  lstar <- file(description = "Post_l_star.csv", open = "a")
  dir_theta <- file(description = "Post_theta.csv", open = "a")																   
  time <- Sys.time()
  Acc <- 1
  # covariate values at the locations of observations
  # Z <- Map(function(cov, lo){
  #   mapply(covariate_function, lo$x, lo$y, cov, grid)
  # }, covariate_process, loc)
  Z <-Map(function(cov, lo){cov[lo]}, covariate_process, loc)
  
  # extremes of integration 
  xl <- covariate_process[[1]]$xrange[1]
  yl <- covariate_process[[1]]$yrange[1]
  xu <- covariate_process[[1]]$xrange[2]
  yu <- covariate_process[[1]]$yrange[2]
  # parameters for spatial integration 
  area <- (yu-yl)*(xu-xl)
  range <- max(d_k) - min(d_k)
  
  # browser()
  index_dom <- lapply(covariate_process, function(c){pmin(round((c$v- min(d_k))*K/range), K)})				   
  index_loc <- lapply(Z, function(c){pmin(round((c- min(d_k))*K/range), K)})		
  
  # browser()
  
  # discretization grid for the GP prior on rho
  mean <- rep(0, K)
  # covariance matrix (with jitter for non singularity)
  jitter <- 1e-10 * diag(1,K) 
  cov <- outer(d_k, d_k, kernel, length_scale, sigma_2, nu) + jitter
  cov <- (cov + t(cov))/2
  Chol <- Matrix::chol(cov)
  # log-GP prior on rho
  log_rho_cur <- (rnorm(K, 0, 1) %*% Chol + mean)
  # length_scale <- rexp(1, rate = 1)
  N <- sum(sapply(loc, npoints))
  lambda_star <- rgamma(1, shape, rate)
  for(i in 1:niter){
    # pCN proposal
    # sample form the prior 
    log_rho_prior <- (rnorm(K, 0, 1) %*% Chol + mean)
    # proposal (log scale)
    log_rho_prop <- sqrt(1 - beta^2) * log_rho_cur + beta * log_rho_prior
    lr <- lik_ratio(d_k, K, rho_new = log_rho_prop, rho_old = log_rho_cur, 
                    Z, grid, covariate_process, lambda_star, area, index_loc, index_dom)
    if(runif(1,0,1) < lr){
      log_rho_cur <- log_rho_prop 
      Acc <- Acc+1
	}
	# write.table(exp(log_rho_cur), file = gp, sep=',', append = TRUE, quote = FALSE,
    #             col.names = FALSE, row.names = FALSE)
    write.table(lambda_star * 1/(1 + exp(-log_rho_cur)), file = gp, sep=',', append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
    # max lambda parameter update
    # rho_eval <- log(lambda_star * plogis(approx(d_k, log_rho_cur, xout = z)$y))
    I_j <- lapply(index_dom, function(i){sum(1/(1 + exp(-log_rho_cur[i])), na.rm = TRUE)*area/length(grid$x)})
    # log_lik <- sum(rho_eval) - lambda_star * I
    # browser()
    lambda_star <- rgamma(1, shape = shape + N, scale = 1/(rate + sum(unlist(I_j))))
    write.table(lambda_star, file = lstar, sep=',', append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    
    # length_scale parameter update
    if(i > niter){
      adapt_cov <- 2.4^2 * (cov_ls + 0.01) 
    }else{
      adapt_cov <- 0.2
																			 
																	 
																	
										   
										
																				   
							   
    }
        
    res <- ind_anisotropic_proposal(Acc, d_k, log_rho_cur, length_scale,
                                    theta, sigma_2, nu, K, Chol, cov,
                                    mean, kernel, jitter, alpha1, alpha2, adapt_cov)
    length_scale <- res[[1]]
    cov <- res[[2]]
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

      
  # browser()
    write.table(t(theta), file = dir_theta, sep=',', append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
	
    write.table(length_scale, file = ls, sep=',', append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    # cat(paste("\rLikelihood is", log_lik, "at iteration", i, Acc))
    cat(paste("\riteration", i, "accept rate is", Acc/i, "max intensity is", lambda_star))
  }
  print(difftime(Sys.time(), time, units = "secs"))
  close(lstar)
  close(ls)
  close(gp)
  close(dir_theta)		   
}
############################ IND ANISOTROPIC ###############################
ind_anisotropic_proposal <- function(Acc, d_k, gp_cur, length_scale,
                                    theta, sigma_2, nu, K, Chol, Cov,
                                    mu, kernel, jitter, alpha1, alpha2, adapt_cov){
  
  length_scale_prop <- exp(log(length_scale) + rnorm(1, 0, adapt_cov)) 
  
  sigma_prop2 <- outer(d_k, d_k, kernel, length_scale_prop, sigma_2, nu) + jitter 
  sigma_old <- Cov
  
  MH_acc <- (alpha1/theta - 1) * log(length_scale_prop) -  
    alpha2 * (length_scale_prop^(1/theta)) -
    (alpha1/theta- 1) * log(length_scale) +  
    alpha2*(length_scale^(1/theta)) + 
    dmvn(gp_cur[1,], mu, sigma = sigma_prop2, isChol = FALSE, log = TRUE) -
    dmvn(gp_cur[1,], mu, sigma = Chol, isChol = TRUE, log = TRUE) + 
    dlnorm(length_scale, log(length_scale_prop), adapt_cov, log = TRUE) -
    dlnorm(length_scale_prop, log(length_scale), adapt_cov, log = TRUE)
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    length_scale <- length_scale_prop
    # new covariance matrix
    Cov <- outer(d_k, d_k, kernel, length_scale, sigma_2, nu) + jitter
    Chol <- Matrix::chol(Cov)
    # Acc <- Acc+1
  }
  return(list(length_scale, Cov, Chol, Acc))
}
############################# IND EXPONENTS ############################

ind_theta_proposal <- function(theta, length_scale, exp_param, alpha1, alpha2){
  # browser()
  theta_prop <- rbeta(1, exp_param[1], exp_param[2])
  MH_acc <- (dbeta(theta, exp_param[1], exp_param[2], log = TRUE)) -         # proposal contribution 
    (dbeta(theta_prop, exp_param[1], exp_param[2], log = TRUE)) +
    (alpha1/theta_prop * log(length_scale)) -                                # gamma likelihood kernel
    alpha2*(length_scale^(1/theta_prop)) -
    (alpha1/theta * log(length_scale)) +  
    alpha2*(length_scale^(1/theta)) +                                        # change of measure                                       
    ((exp_param[1]-2)*log(theta_prop) + (exp_param[2]-1)*log(1 - theta_prop)) -   # prior beta kernel
    ((exp_param[1]-2)*log(theta) + (exp_param[2]-1)*log(1 - theta))
  
  if(runif(1,0,1) < min(1, exp(MH_acc))){
    theta <- theta_prop 
  }
  return(theta)
}			   				   