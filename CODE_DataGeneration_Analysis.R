########################################### LOAD LIBRARIES #######################################

library(spatstat)
library(MASS)
library(ggplot2)
library(plotly)
library(zeallot)
library(fourPNO)
library(fMultivar)
library(akima)
library(sn)
library(Matrix)
library(profvis)
library(kernlab)
library(mvtnorm)
library(xtable)
library(pracma)
library(patchwork)
library(viridis)
library(scico)
library(RColorBrewer)
library(pals)

### SET WORKING DIRECTORY TO SOURCE FILE LOCATION

#################### DEFINITION OF OBSERVATION WINDOW AND COVARIATE PROCESS GENERATION ##############

# define observation window as owin object  
window <- owin(xrange = c(0,1), yrange = c(0,1))

# define grid for covariate process evaluation on center of each pixels 
grid_points <- gridcentres(window, nx = 50, ny = 50)

# pairwise distance between points in the grid
distance_matrix <- pairdist(grid_points, squared = TRUE)

# covariance function
exponential_cov <- function(d, length_scale = 1, sigma_2 = 1) {
  sigma_2 * exp(-d / (2 * length_scale))
}

# covariance matrix of the GP generating the covariate process, discretized over the grid
Cov <- exponential_cov(distance_matrix, length_scale = 0.005, sigma_2 = 1) + 1e-10 * diag(1,length(grid_points$x))
Chol <- Matrix::chol(Cov)
plot(im(Cov[dim(Cov)[1]:1,]))

# mean function 
mu = rep(0, length(grid_points$x))

# N replicates of covariate process Z_1 
M <- 1e+01
cov1_list <- list()
dim = length(grid_points$x)
for(i in 1:M){
  covariate_process <- as.vector(rnorm(n = dim, mean = 0, sd = 1) %*% Chol)
  covariate_process <- pnorm((covariate_process - mean(covariate_process))/sd(covariate_process))
  int_i <- im(matrix(covariate_process, sqrt(dim), sqrt(dim))[sqrt(dim):1,],
                            xrange = c(0,1), yrange = c(0,1))
  cov1_list[[i]] <- int_i
}

#################### UNIVARIATE INTENSITY FUNCTIONS #############################

# Choose among next 3 functions 
{# negative exponential
  rho <- function(z){
    2 * exp(3*(1-z)-1)
  }
}
{# skew normal bump 
  rho <- function(z){
    5 * dsn(z, xi = 0.8, omega = 0.3, alpha = -5)
  }
}
# negative/positive deviation from plateau
{ 
  smootherstep <- function(x) {
    x <- pmin(pmax(x, 0), 1)
    6*x^5 - 15*x^4 + 10*x^3
  }
  
  # periodic distance on [0,1] from center c
  periodic_dist <- function(z, c) {
    d <- abs(((z - c + 0.5) %% 1) - 0.5)
    d
  }
  
  # plateau builder: returns value in [0,1]; ~=1 near center, ~=0 outside
  plateau <- function(z, center = 0.5, width = 0.1) {
    # width = full plateau width (value ~1 for |d| <= width/2)
    d <- periodic_dist(z, center)
    t <- d / (width/2)       # normalized distance:
    1 - smootherstep(t)      # = 1 when d=0, ->0 as d >= width/2
  }
  
  # final function: baseline 5, min plateau at 1/4 => 0, max plateau at 3/4 => 15
  rho <- function(z, width = 3/8) {
    # ensure z in [0,1] (vectorized)
    z <- z %% 1
    baseline <- 2
    # bump amplitudes:
    amp_max <- 4 - baseline   # +10
    amp_min <- baseline - 0    # 5 (we will subtract this plateau to reach 0)
    baseline + amp_max * plateau(z, center = 3/4, width = width) -
      amp_min * plateau(z, center = 1/4, width = width)
  }
}

# plot of true intensity 
plot(seq(0, 1, length.out = 200), 
     rho(seq(0, 1, length.out = 200)), 
     type ="l", main = "True intensity function rho", ylab = "intensity", xlab = "covariate values")


#################### POINT PATTERN GENERATION ##################################
loc_list <- list()
for(i in 1:M){
  int_i <- cov1_list[[i]]
  loc_list[[i]] <- rpoispp(lambda = eval.im(rho(int_i)), win= window, nsim = 1)
}

#################### Figure 1: covariate process ################################
plot_list_cov <-list()
for(i in 1:3){
covariate_process <- cov1_list[[i]]
loc <- loc_list[[i]]

covariate_process_ppobj <- ppp(x = grid_points$x, y = grid_points$y,
                               marks = as.vector(t(covariate_process$v)), window = window)
covariate_process_plot <- Smooth(covariate_process_ppobj, kernel = "gaussian", dimyx = c(200,200), method = "pwl")


img_df <- data.frame(value = as.vector(covariate_process_plot$v))
img_df$x <- gridcentres(window, nx = 200, ny = 200)$y-0.5
img_df$y <- gridcentres(window, nx = 200, ny = 200)$x-0.5

points <- data.frame(x = loc$x-0.5,
                     y = loc$y-0.5)

plot_list_cov[[i]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  
  scale_fill_gradientn(colors =colorRampPalette(coolwarm(11))(100), name = expression(Z[1](x))) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = expression(x[1]), y = expression(x[2])) +
  geom_point(data = points, aes(x = x, y = y), shape = 19, color="black", size = 1, alpha = 1) 
}
combined_plot <- wrap_plots(plot_list_cov, nrow = 1)
combined_plot

#################### FREQUENTIS INFERENCE for INTENSITY FUNCTION ###########################
# Average kernel Eq 14 
{
  kernel_function <- function(loc_i, cov_i, N) {
    res <- tryCatch(
      {
        rhohat(loc_i, cov_i, n = 200)
      },
      error = function(e) {
        data.frame(rho = rep(NA,200),
                   hi = rep(NA,200),
                   lo = rep(NA,200))  
      }
    )
    return(cbind(res$rho, res$hi, res$lo))
  }
  
  # L2 errors and relative counterpart across 100 experiments 
  # (need M = 100*max(N) covariate process and events replicates)
  L2_exp = matrix(NA,4,100)
  L2_rel_exp = matrix(NA,4,100)
  i = 1
  discr <- seq(0,1, length.out = 200)
  for(N in c(50,250,500,1000)){
    for(exp in 1:100){
      multi_kernel <- Map(function(loc_i, cov_i){kernel_function(loc_i, cov_i)}, 
                          loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], 
                          cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)])
      arr <- array(unlist(multi_kernel), dim = c(200, 3, length(multi_kernel)))
      avg_matrix <- apply(arr, c(1, 2), mean, na.rm = TRUE)
      L2_exp[i,exp]<-(sqrt(sum((avg_matrix[,1] - rho(discr))^2)/length(discr)))
      L2_rel_exp[i,exp]<-(sqrt(sum((avg_matrix[,1] - rho(discr))^2)/(sum(rho(discr)^2))))
    }
  }

  # Avarege L2 errros (and sd) as in Table 1 
  funcs <- list(mean = mean, sd = sd)
  result <- lapply(funcs, function(f) apply(L2_exp, 1, f, na.rm = TRUE))
  print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))
  
  result <- lapply(funcs, function(f) apply(L2_rel_exp, 1, f, na.rm = TRUE))
  print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))
}

#################### BAYESIAN ANALYSIS ################################

L2_bay_exp = matrix(NA,4,100)
L2_bay_rel_exp = matrix(NA,4,100)

for(N in c(50,250,500,1000)){
  setwd("~/")
  source(file = "pCN_MCMC_GP_mult.R")
  folder = paste("Adaptive_MCMC_post_multi/", N, sep = "")
  dir.create(paste("~/", folder, sep = ""))
  setwd(folder)
  discr <- seq(0,1, length.out = 200)
    n_iter <- 20000
  for(exp in 1:100){
    setwd("~/Research/IPP/")
    source(file = "pCN_MCMC_GP_mult.R")
    adapt_multi_MCMC(paste0(N, "/", exp, sep = ""), 
                     n_iter, beta = 0.08, K = length(discr), 
                     loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], grid_points, 
                     cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], d_k = discr,
                     squared_exponential_kernel, alpha1 = 1, alpha2 = 1, 
                     sigma_2 = 1, nu = 3/2, exp_param = c(2,2),
                     shape = 1, rate = 2, "lambda")
    
    # upload posterior draws for rho
    gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
    gs <- apply(gs, 2, as.numeric)
    post_mean <- colMeans(gs[seq(n_iter/4,n_iter),],na.rm = TRUE)
    
    L2_bay_exp[i,exp]<- sqrt(sum(post_mean - rho(discr))^2)/length(discr)
    L2_bay_rel_exp[i,exp]<- sqrt(sum((post_mean - rho(discr))^2)/(sum(rho(discr)^2)))
  }
}   

# Avarege L2 errros (and sd) as in Table 1 
funcs <- list(mean = mean, sd = sd)
result <- lapply(funcs, function(f) apply(L2_bay_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))

result <- lapply(funcs, function(f) apply(L2_bay_rel_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))

#################### Figure 2: functional estimates ################################ 
plot_list <-list()
for(i in 1:3){
  N <- c(250,500,1000)[i]
  exp <- 1
  rug <- unlist(Map(function(c,l){return(c[l])} , cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], 
                    loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)]))
  multi_kernel <- Map(function(loc_i, cov_i){kernel_function(loc_i, cov_i)}, 
                      loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], 
                      cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)])
  arr <- array(unlist(multi_kernel), dim = c(200, 3, length(multi_kernel)))
  avg_matrix <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  
  spline_fit_est <- smooth.spline(discr, avg_matrix[,1], spar = 0.8)
  
  folder = paste("~/Adaptive_MCMC_post_multi/neg_exp/", M, "/", exp, sep = "")
  setwd(folder)
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs <- apply(gs, 2, as.numeric)
  post_mean <- colMeans(gs[seq(5/10*dim(gs)[1],dim(gs)[1]),],na.rm = TRUE)
      
  upper = apply(gs[seq(3/10*dim(gs)[1],dim(gs)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
  lower = apply(gs[seq(3/10*dim(gs)[1],dim(gs)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)
      
  plot_list[[i]] <- 
    ggplot(data.frame(x=discr, y=rho(discr)), aes(x = x, y = y)) +
    geom_line(linetype = "solid", color = "black", size = 1) +
    theme_minimal(base_family = "sans") +  # Clean white background
    theme(
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      text = element_text(color = "black")
    ) + labs(x = expression(z), y = expression(rho(z))) + 
    geom_line(data = data.frame(x=discr, y=predict(spline_fit_est, x = discr)$y), mapping = aes(x = x, y = y), 
              color = "firebrick", size = 1) + ylim(0,16) + 
    geom_rug(data = data.frame(x = rug, y = 0), sides = "b") + ylim(0,10) + 
    geom_rug(data = data.frame(x = rug, y = 0), sides = "b") + 
    geom_line(data = data.frame(x=discr, y=post_mean), mapping = aes(x = x, y = y),
              color = "dodgerblue3", size = 1) +
    geom_ribbon(data = data.frame(x = discr, y2 = upper, y1 = lower, y=post_mean), mapping = aes(x = x, ymin = y1, y = post_mean, ymax = y2),
                fill = "dodgerblue3", alpha = 0.2)
}
combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot

#################### CONVERGENCE DIAGNOSTICS #######################################
col = brewer.pal(12,"Paired")
trace_plot1 <-  ggplot() + theme_minimal(base_family = "sans") +  
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) 
trace_plot2 <- trace_plot1
trace_plot3 <- trace_plot1
trace_plot4 <- trace_plot1
trace_plot5 <- trace_plot1
trace_plot6 <- trace_plot1
trace_plot7 <- trace_plot1

for(exp in c(2,4,6,8,10,12)){
  folder = paste("~/Adaptive_MCMC_post_multi/1000", "/", exp, sep = "")
  setwd(folder)
  
  lambda_post = read.table("Post_l_star.csv", sep = ",")
  theta = read.table("Post_theta.csv", sep = ",")
  Ls = read.table("Post_ls.csv", sep = ",")

  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs <- apply(gs, 2, as.numeric)
  
  Z <-Map(function(cov, lo){cov[lo]},
          cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
          loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  index_dom <- lapply(cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
                      function(c){pmin(round(c$v*200), 200)})
  index_loc <- lapply(Z, function(c){pmin(round(c*200), 200)})
  log_lik <- apply(as.data.frame(gs[seq(1,n_iter,1),]), 1, function(rho){
    rho_obs <- log(rho[unlist(index_loc)])
    I <- lapply(index_dom, function(i){sum(rho[i]-1, na.rm = TRUE)/2500})
    sum(rho_obs) - sum(unlist(I))
  })
  true_lik <- sum(log(rho(discr)[unlist(index_loc)])) -
    sum(unlist(lapply(index_dom, function(i){sum(rho(discr)[i]-1, na.rm = TRUE)/2500})))
  
  # Figure 15
  trace_plot1 <- trace_plot1 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(lambda_post)[1:n_iter]), mapping = aes(x = x, y = y, color = col[exp]),
  linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho^"*"))
  trace_plot2 <- trace_plot2 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(Ls)[1:n_iter]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = "ℓ")
  trace_plot3 <- trace_plot3 + geom_line(data = data.frame(x=seq(2:n_iter), y=t(log_lik)[2:n_iter]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3) + labs(x = "Number of iterations", y = "log-likelihood") +
                               geom_hline(yintercept = true_lik, color = col[exp], linetype = "dashed", size = 0.5)

  trace_plot4 <- trace_plot4 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(theta)[1:n_iter]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(theta))
  # Figure 16
  trace_plot5 <- trace_plot5 + geom_line(data = data.frame(x=seq(1:n_iter), y=gs[1:n_iter,130]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho [(0.65)])) +
    geom_hline(yintercept = rho(discr[130]), color = "black", linetype = "dashed", size = 0.5)
  
  trace_plot6 <- trace_plot6 + geom_line(data = data.frame(x=seq(1:n_iter), y=gs[1:n_iter,31]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho [(0.15)])) +
    geom_hline(yintercept = rho(discr[31]), color = "black", linetype = "dashed", size = 0.5)
  
  trace_plot7 <- trace_plot7 + geom_line(data = data.frame(x=seq(1:n_iter), y=gs[1:n_iter,190]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho[(0.95)])) +
    geom_hline(yintercept = rho(discr[190]), color = "black", linetype = "dashed", size = 0.5)
}

combined_plot <- wrap_plots(trace_plot5, trace_plot6, trace_plot7, nrow = 1)
combined_plot

#################### DEFINITION OF 2-DIMENSIONAL COVARIATE PROCESS ###########################
# Smaller length scale 
Cov2 <- exponential_cov(distance_matrix, length_scale = 0.05, sigma_2 = 1) + 1e-10 * diag(1,length(grid_points$x))
Chol2 <- Matrix::chol(Cov2)

cov2_list <- list()
for(i in 1:M){
  covariate_process <- as.vector(rnorm(n = dim, mean = 0, sd = 1) %*% Chol2)
  covariate_process <- pnorm((covariate_process - mean(covariate_process))/sd(covariate_process))
  int_i <- im(matrix(covariate_process, sqrt(dim), sqrt(dim))[sqrt(dim):1,],
              xrange = c(0,1), yrange = c(0,1))
  cov2_list[[i]] <- int_i
}

#################### BIVARIATE INTENSITY FUNCTIONS #############################

# Choose among next 3 functions 

{ # isotropic double bump function
  rho1 <- function(z1, z2){
    14 * dmsn(cbind(z1, z2), xi=c(0.7, 0.2), Omega = diag(0.05, nrow = 2), alpha = c(3,-2), tau=1, dp=NULL, log=FALSE)
    }
  rho2 <- function(z1, z2){
    6 * dmsn(cbind(z1, z2), xi=c(0.3, 0.8), Omega = diag(0.03, nrow = 2), alpha = c(-1,-1), tau=0, dp=NULL, log=FALSE)
  }
  rho2d <- function(z1, z2){
    pmax(0, rho1(z1,z2) + rho2(z1,z2))
  }
}

{ # monotonic slope function 
  rho1 <- function(z1, z2){
    90 * dmsn(cbind(z1, z2), xi=c(0.3, 0.3), Omega = diag(0.5, nrow = 2), alpha = c(-1, -1), tau=1, dp=NULL, log=FALSE)
  }
  rho2d <- function(z1, z2){
    pmax(0, 30 - rho1(z1,z2))
  }
}

{ # anisotropic minimum and maximum 
    rho2d <- function(x, y) {
    baseline <- 10
    
    # bump parameters (positive Gaussian)
    A_b <- baseline
    x_b <- 0.3; y_b <- 0.8
    sx_b <- 0.08; sy_b <- 0.50
    
    # hole parameters (negative Gaussian, deep enough to reach 0)
    A_h <- baseline   # ensures depth reaches zero
    x_h <- 0.8; y_h <- 0.3
    sx_h <- 0.08; sy_h <- 0.50   # parallel anisotropy (same orientation as bump)
    
    # helper: anisotropic Gaussian
    gauss2d <- function(x, y, A, x0, y0, sx, sy) {
      A * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2)))
    }
    
    bump <- gauss2d(x, y, A_b, x_b, y_b, sx_b, sy_b)
    hole <- gauss2d(x, y, A_h, x_h, y_h, sx_h, sy_h)
    
    f_raw <- baseline + bump - hole
    
    # clip at 0 (non-smooth floor)
    pmax(f_raw, 0)
  }
}

# plot of true intensity 
exp_cov <- expand.grid(x = seq(0,1, length.out = 800), 
            y = seq(0,1, length.out = 800))
rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                      sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
               xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
plot(rho_imag, axes=TRUE, main = "Rho function", xlab = "First covariate", 
     ylab = "Second covariate")


#################### POINT PATTERN GENERATION #############################

loc2d_list <- list()
for(i in 1:M){
  loc2d_list[[i]] <- rpoispp(lambda = im(matrix(rho2d(as.vector(cov1_list[[i]]$v), 
                                                      as.vector(cov2_list[[i]]$v)),
                                                50,50),
                                         xrange = c(0,1), yrange = c(0,1)),
                                         win= window, nsim = 1)
}

#################### TRIANGULAR MESH CONSTRUCTION #################################

# convex hull as whole domain 
points <- matrix(c(0,0,0,1,1,1,1,0), ncol = 2, byrow = T)
conv_hull <- chull(points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))

# triangular tesseletion of the convex hull with maximum volume a = 0.0014

tri_tess <- RTriangle::pslg(P = cbind(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y))
mesh <- RTriangle::triangulate(tri_tess, a = 0.0014, j = TRUE)  #Y = TRUE, D = TRUE # smaller 'a' = finer mesh
triangle_centroids <- t(apply(mesh$T, 1, function(tri_idx) {
  colMeans(mesh$P[tri_idx, ])
}))

triangle_area <- function(a, b, c) {
  0.5 * abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]))
}
triangle_area <- apply(mesh$T, 1, function(tri) {
  a <- mesh$P[tri[1], ]
  b <- mesh$P[tri[2], ]
  c <- mesh$P[tri[3], ]
  triangle_area(a, b, c)
})

plot(mesh, asp = 1, main = "Uniform Triangular Tessellation of Convex Hull")
vertices <- mesh$P
extremes <- mesh$P[mesh$PB==1,]
colnames(extremes) <- c("x", "y")
colnames(vertices) <- c("x", "y")

# Convert to data frame for plotting
triangles_df <- as.data.frame(do.call(rbind, lapply(1:nrow(mesh$T), function(i) {
  tri <- mesh$T[i, ]  # Indices of the triangle's 3 vertices
  coords <- mesh$P[tri, ]  # Get x, y coords of those vertices
  
  # Add triangle ID and ensure polygons are closed
  data.frame(
    triangle_id = i,
    x = c(coords[, 1], coords[1, 1]),  # close the triangle
    y = c(coords[, 2], coords[1, 2])
  )
})))

#################### FREQUENTIS INFERENCE for INTENSITY FUNCTION ###########################

# average kernel 
kern2_est <- rho2hat(loc2d_list[[1]], cov1_list[[1]], cov2_list[[1]], method = "ratio", dimyx = c(800,800))

# L2 errors and relative counterpart across 100 experiments 
# (need M = 100*max(N) covariate process and events replicates)
L2d_exp = matrix(NA,4,100)
L2d_rel_exp = matrix(NA,4,100)
i = 1
for(N in c(10,50,250,1000)){
  for(exp in 1:100){
      kern_years <- lapply(c(((exp-1)*1000 + 1):((exp-1)*1000 + N)), function(i){rho2hat(loc2d_list_slope[[i]], 
                              cov1_list[[i]], cov2_list[[i]], method = "ratio", dimyx = c(200,200))$v})
      kern2_est$v <- apply(simplify2array(kern_years), c(1, 2), mean, na.rm = TRUE)
      
      L2d_exp[i,exp]<-(sqrt(sum((kern2_est$v - rho_imag)^2)/length(exp_cov$x)))
      L2d_rel_exp[i,exp]<-(sqrt(sum((kern2_est$v - rho_imag)^2)/sum(rho_imag^2)))
  }
  i <- i+1
}

# Avarege L2 errros (and sd) as in Table 2 
funcs <- list(mean = mean, sd = sd)
result_rel <- lapply(funcs, function(f) apply(L2d_rel_exp, 1, f, na.rm = TRUE))
print(paste(round(result_rel$mean,2), "(", round(result_rel$sd,2), ")&"))
result <- lapply(funcs, function(f) apply(L2d_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,2), "(", round(result$sd,2), ")&"))

#################### BAYESIAN ANALYSIS ################################

L2_2dbay_exp = matrix(NA,5,10)
L2_2dbay_rel_exp = matrix(NA,5,10)

bet = c(0.3, 0.15, 0.1, 0.08)

for(N in c(10,50,250,1000)){
  setwd("~/")
  source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
  folder = paste("2dAdaptive_MCMC_post_mesh_rep/", N, sep = "")
  dir.create(paste("~/Research/IPP/", folder, sep = ""))
  setwd(folder)
  for(exp in 1:100){
    setwd("~/Research/IPP/")
    source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")

    MCMC(paste("2dAdaptive_MCMC_post_mesh_rep/",N, "/", exp,sep = ""), n_iter, beta = bet[i],
       K = dim(vertices)[1], loc2d_list_doublebump[((exp-1)*1000 + 1):((exp-1)*1000 + N)], 
       grid_points, window,
       cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], 
       cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)],
       as.data.frame(vertices), triangle_centroids,
       exponential_cov, alpha1=1, alpha2=1, sigma_2 = 1, nu = NA,
       exp_param = c(2,2), shape = 1, rate = 2, 
       link = "lambda", LS = "anisotropic")
    
    gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
    gs_m <- apply(gs_m, 2, as.numeric)
    pp <- ppp(x = vertices[,1], y = vertices[,2],
              marks = colMeans(gs_m[seq(5/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
    post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
    
    L2_2dbay_exp[i,exp]<- sqrt(sum((post_mean_m - rho_imag)^2)/length(exp_cov$x))
    L2_2dbay_rel_exp[i,exp]<-  sqrt(sum((post_mean_m - rho_imag)^2))/sqrt(sum((rho_imag)^2))
  }
}

# Avarege L2 errros (and sd) as in Table 2
funcs <- list(mean = mean, sd = sd)
result_rel <- lapply(funcs, function(f) apply(L2_2dbay_rel_exp, 1, f, na.rm = TRUE))
print(paste(round(result_rel$mean,2), "(", round(result_rel$sd,2), ")&"))
result <- lapply(funcs, function(f) apply(L2_2dbay_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,2), "(", round(result$sd,2), ")&"))
  

#################### Figure 3: functional estimates ################################ 
plot_list2 <- list()
i = 1
for(M in c(50, 250, 1000)){
  setwd("~/")
  folder = paste("2dAdaptive_MCMC_post_mesh_repl/", M, sep = "")
  setwd(folder)
  exp = 1
  folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_repl/", M, "/", exp, sep = "")
  setwd(folder)
  
  gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs_m <- apply(gs_m, 2, as.numeric)
  pp <- ppp(x = vertices[,1], y = vertices[,2],
            marks = colMeans(gs_m[seq(5/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
  post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
  
  exp_cov <- expand.grid(x = seq(1.084202e-19, 1, length.out = 800), 
                         y = seq(1.084202e-19, 1, length.out = 800))
  
  img_df <- data.frame(value = as.vector(post_mean_m$v))
  img_df$x <- exp_cov$y
  img_df$y <- exp_cov$x
  
  plot_list2[[i]] <-
    ggplot() +
    geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  
    scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"), 
                         name = expression(rho(z))) + # , limits = c(0, 55)
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = expression(z[1]), y = expression(z[2])) + theme(legend.position = "none")
  
  i<-i+1
}

# add ground truth 
img_df <- data.frame(value = as.vector(rho_imag$v))
img_df$x <- gridcentres(window, nx = 200, ny = 200)$y
img_df$y <- gridcentres(window, nx = 200, ny = 200)$x

plot_list2[[4]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) + #, limits = c(0, 55)
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = expression(z[1]), y = expression(z[2]))


combined_plot <- wrap_plots(plot_list2, nrow = 1)
combined_plot

#################### CONVERGENCE DIAGNOSTICS #######################################
trace_plot1 <-  ggplot() + theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) 
trace_plot2 <- trace_plot1
trace_plot3 <- trace_plot1
trace_plot4 <- trace_plot1
trace_plot5 <- trace_plot1
trace_plot6 <- trace_plot1

trace_plot7 <- trace_plot1
trace_plot8 <- trace_plot1
trace_plot9 <- trace_plot1

# find index for point-wise log-likelihood evaluation
which.min(rowSums((sweep(vertices, 2, c(0.8,0.3), "-"))^2))
which.min(rowSums((sweep(vertices, 2, c(0.3,0.8), "-"))^2))
which.min(rowSums((sweep(vertices, 2, c(0.5,1), "-"))^2))

RHO2 <- rho2d(vertices[,1], vertices[,2])

for(exp in c(1:5)){
  N = 1000
  setwd("~/Research/IPP/2dAdaptive_MCMC_post_mesh_repl/1000")
  folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_repl/1000", "/", exp, sep = "")
  setwd(folder)
  
  lambda_post = read.table("Post_l_star.csv", sep = ",")
  theta = read.table("Post_theta.csv", sep = ",")
  Ls = read.table("Post_ls.csv", sep = ",")
  
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs <- apply(gs, 2, as.numeric)
   
  Z1 <-Map(function(cov, lo){cov[lo]},
          cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)],
          loc2d_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)])
  Z2 <-Map(function(cov, lo){cov[lo]},
           cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)],
           loc2d_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)])

  list_index_dom <- Map(function(cov1, cov2){unlist(apply(cbind(na.omit(as.vector(cov1$v)),na.omit(as.vector(cov2$v))),
                                                          1, function(w){which.min(rowSums((vertices - matrix(w, nrow = dim(vertices)[1], ncol = 2, byrow = TRUE))^2))}))},
                        cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)], cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + N)])
  list_index_loc <- Map(function(z1, z2){apply(cbind(z1,z2), 1,
                        function(w){which.min(rowSums((vertices - matrix(w, nrow = dim(vertices)[1], ncol = 2, byrow = TRUE))^2))})},
                        Z1, Z2)
  
  log_lik <- rbind(log_lik,apply(as.data.frame(gs), 1, function(r){
    rho_obs <- log(r[unlist(list_index_loc)])
    I <- lapply(list_index_dom, function(i){sum(r[i]-1, na.rm = TRUE)/2500})
    (sum(rho_obs) - sum(unlist(I)))/1000
  }))
  
  true_lik[exp] <- (sum(log(RHO2[unlist(list_index_loc)])) -
    sum(unlist(lapply(list_index_dom, function(i){sum(RHO2[i]-1, na.rm = TRUE)/2500}))))/1000
  
  # max intensity
  trace_plot1 <- trace_plot1 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(lambda_post)[1:n_iter]), mapping = aes(x = x, y = y),
                      linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho^"*"))
  # length scales
  trace_plot2 <- trace_plot2 + geom_line(data = data.frame(x=seq(1:n_iter), y=Ls[1:n_iter,1]), mapping = aes(x = x, y = y),
                      linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(ℓ[1]))
  trace_plot4 <- trace_plot4 + geom_line(data = data.frame(x=seq(1:n_iter), y=Ls[1:n_iter,2]), mapping = aes(x = x, y = y),
  linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(ℓ[2]))
  # log-likelihood
  trace_plot3 <- trace_plot3 + geom_line(data = data.frame(x=seq(2:n_iter), y=log_lik[exp+6,2:n_iter]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 1) + labs(x = "Number of iterations", y = "log-likelihood") +
                               geom_hline(yintercept = true_lik[exp], color = col[exp], linetype = "dashed", size = 1)
  # length scale exponents
  trace_plot5 <- trace_plot5 + geom_line(data = data.frame(x=seq(1:n_iter), y=theta[1:n_iter,1]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(theta[1]))
  trace_plot6 <- trace_plot6 + geom_line(data = data.frame(x=seq(1:n_iter), y=theta[1:n_iter,2]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(theta[2]))
  # point-wise evaluations of the log-lik
  trace_plot7 <- trace_plot7 + geom_line(data = data.frame(x=seq(1:n_iter), y=gs[1:n_iter,546]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho[(0.8 * "," ~ 0.3)])) +
    geom_hline(yintercept = rho2d(vertices[546,1], vertices[546,2]), color = "black", linetype = "dashed", size = 0.5)
  
  trace_plot8 <- trace_plot8 + geom_line(data = data.frame(x=seq(1:n_iter), y=gs[1:n_iter,354]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho[(0.3 * "," ~ 0.8)])) +
    geom_hline(yintercept = rho2d(vertices[354,1], vertices[354,2]), color = "black", linetype = "dashed", size = 0.5)
  
  trace_plot9 <- trace_plot9 + geom_line(data = data.frame(x=seq(2:n_iter), y=gs[2:n_iter,135]), mapping = aes(x = x, y = y),
                                         linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho[(0.5 * "," ~ 1)])) +
    geom_hline(yintercept = rho2d(vertices[135,1], vertices[135,2]), color = "black", linetype = "dashed", size = 0.5)
  
}

trace_plot3  <- trace_plot3 + labs(x = "Number of iterations", y = expression(log-likelihood ~ (x ~ 10^3)))

combined_plot <- wrap_plots(trace_plot2, trace_plot4, nrow = 1)
combined_plot


#################### POSTERIOR PROJECTIONS fig 13 ############################# 
marg_list <- list()
i <- 1
for(N in c(50, 250, 1000)){
  
  setwd("~/")
  folder = paste("2dAdaptive_MCMC_post_mesh_repl/", N, sep = "")
  setwd(folder)
  
  exp = 1
  
  folder = paste("~/2dAdaptive_MCMC_post_mesh_repl/", N, "/", exp, sep = "")
  setwd(folder)
  gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs_m <- apply(gs_m, 2, as.numeric)
  pp <- ppp(x = vertices[,1], y = vertices[,2],
            marks = colMeans(gs_m[seq(5/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
  post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
  
  exp_cov <- expand.grid(x = seq(1.084202e-19, 1, length.out = 800), 
                         y = seq(1.084202e-19, 1, length.out = 800))
  
  rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                        sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
                 xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
  
  img_df <- data.frame(value = as.vector(post_mean_m$v))
  img_df$x <- exp_cov$y
  img_df$y <- exp_cov$x
  
  # change commented "y" for different diagonal 
  x = 1:800
  y = x 
  # y = 801 - x
  
  pp_upper <- ppp(x = vertices[,1], y = vertices[,2],
            marks = apply(gs_m[seq(3/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.9999, na.rm = TRUE), window = owin(poly_win))
  
  upper = Smooth(pp_upper, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
  
  pp_lower <- ppp(x = vertices[,1], y = vertices[,2],
                  marks = apply(gs_m[seq(3/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.01, na.rm = TRUE), window = owin(poly_win))
  
  lower = Smooth(pp_lower, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
  
  plot_M <- ggplot(data.frame(x=seq(0,1,length.out = 800), y=post_mean_m$v[cbind(x,y)]), aes(x = x, y = y)) +
    geom_line(linetype = "solid", color = "dodgerblue3", size = 1) +
    theme_minimal(base_family = "sans") +  # Clean white background
    theme(
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      text = element_text(color = "black")
    ) + labs(x = expression(z[1]), y = expression(rho(z[1], z[1]))) + 
    geom_line(data = data.frame(x=seq(0,1,length.out = 800), y=rho_imag$v[cbind(x,y)]), mapping = aes(x = x, y = y), 
              color = "black", linetype = "solid", size = 1) + 
    geom_ribbon(data = data.frame(x = seq(0,1,length.out = 800), y2 = upper$v[cbind(x,y)], 
                                                                 y1 = lower$v[cbind(x,y)], 
                                                                 y=post_mean_m$v[cbind(x,y)]), 
                mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "dodgerblue3", alpha = 0.2)
  marg_list[[i]] <- plot_M
  i <- i+1
}

combined_plot <- wrap_plots(marg_list, nrow = 1)
combined_plot
