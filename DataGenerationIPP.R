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
# define observation window as owin object  
window <- owin(xrange = c(0,1), yrange = c(0,1))
# define grid for covariate process evaluation
# centers of the pixels 
grid_points <- gridcentres(window, nx = 50, ny = 50)
# grid_points$x <- jitter(grid_points$x, amount = 0.02)
# grid_points$y <- jitter(grid_points$y, amount = 0.02)
# extremes of the pixels
# grid_points <- expand.grid(x = seq(0,1,length.out = 3), y = seq(0,1,length.out = 3))
# ppp_grid <- as.ppp(grid_points, W = owin(c(0, 1), c(0, 1)))
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

# with matern kernel instead
{
Cov <- matern32_kernel(distance_matrix, length_scale = 0.5, sigma_2 = 1, nu = 1.5) 
plot(im(Cov[dim(Cov)[1]:1,]))
Chol <- Matrix::chol(Cov)
# alternative with kernel function form kernlab
length_scale <- 0.005
rbf_kernel <- rbfdot(sigma = 1 / (2 * length_scale))
Cov2 <- kernelMatrix(rbf_kernel, matrix(c(grid_points$x, grid_points$y), nrow = length(grid_points$x)))
plot(im(Cov2[dim(Cov2)[1]:1,]))
# matern counterpart 
Cov2 <- Matern(sqrt(distance_matrix), alpha = sqrt(2 * 1.5)/0.5, smoothness = 1.5, phi = 1) +
        1e-10 * diag(1,length(grid_points$x)) #alpha = sqrt(2 * 1.5)/0.2
plot(im(Cov2[dim(Cov2)[1]:1,]))
# mean function 
}
mu = rep(0, length(grid_points$x))
#realization of the covariate process
{
covariate_process1 <- im(matrix(rnorm(n = length(grid_points$x), mean = 0, sd = 1) %*% Chol + mu, 
                               sqrt(length(grid_points$x)), sqrt(length(grid_points$x)))[sqrt(length(grid_points$x)):1,],
                      xrange = c(0,1), yrange = c(0,1))
# standardization and [0,1] mapping of the covariate process 
covariate_process1 <- ((covariate_process1 - mean(covariate_process1))/sd(covariate_process1))
covariate_process1 <- eval.im(pnorm(covariate_process1))

covariate_process_ppobj <- ppp(x = grid_points$x, y = grid_points$y,
          marks = as.vector(t(covariate_process1$v)), window = window)
covariate_process_plot <- Smooth(covariate_process_ppobj, kernel = "gaussian", dimyx = c(200,200), method = "pwl")

plot(covariate_process_plot, main = "Covariate process on unit square")
contour(covariate_process1, nlevels = 10)
plot(levelset(covariate_process1, 0, ">"))

plot_ly(x = ~seq(0,1,length.out = 100), y = ~seq(0,1,length.out = 100), z = ~covariate_process1$v, 
        type = 'surface', colorscale = 'Viridis', showscale = TRUE) %>%
  layout(scene = list(
    xaxis = list(title = 'X'),
    yaxis = list(title = 'Y'),
    zaxis = list(title = 'GP Value')
  ))
# single observation, inflated rho 

# intensity function from covariate space to [0, infty]
rho <- function(z){
  # exponential, tuning for reasonable intensity
  2 * exp(3*(1-z)-1)
}
plot(seq(min(covariate_process1)-1, max(covariate_process1)+1, length.out = 200), 
     rho(seq(min(covariate_process1)-1, max(covariate_process1)+1, length.out = 200)), 
     type ="l", main = "True intensity function rho", ylab = "intensity", xlab = "covariate values")
plot(rho(covariate_process1))
}
#realization of point pattern
{
rho <- function(z){
  # bump function
  5 * dsn(z, xi = 0.8, omega = 0.3, alpha = -5)
}
  
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

plot(seq(0, 1, length.out = 200), 
     rho(seq(0, 1, length.out = 200)), 
     type ="l", main = "True intensity function rho", ylab = "intensity", xlab = "covariate values")

intensity_process <- covariate_process1
intensity_process$v <- apply(covariate_process1$v, c(1, 2), function(x){rho((x))})
intensity_process_ppobj <- ppp(x = grid_points$x, y = grid_points$y,
                               marks = as.vector(t(intensity_process$v)), window = window)
intensity_process_plot <- Smooth(intensity_process_ppobj, kernel = "gaussian", dimyx = c(200,200), method = "pwl")

plot(intensity_process_plot, main = "Intensity on unit square")
# functions returning intensity for given coordinates (x,y)
# by returning intensity of the closest grid point
intensity <- function(x, y, r, grid, covariate){
  index <- which.min(sqrt((grid$x - x)^2 + ((grid$y) - y)^2))
  X_index <- which.min(abs(covariate$xcol - grid$x[index]))
  Y_index <- which.min(abs(covariate$yrow - grid$y[index]))
  return(r(covariate[Y_index, X_index]))
}
# by bilinear interpolation over the grid  
intensity2 <- function(x, y, r, grid, covariate){
  index <- which.min(sqrt((grid$x - x)^2 + ((grid$y) - y)^2))
  return(r(interp.im(covariate, x = grid$x[index], y = grid$y[index])))
}
# integral of the intensity function over the observation domain 
integrate2d(intensity2, rho, grid_points, covariate_process1, error = 1.0e-5)$value
integrate2d(intensity, rho, grid_points, covariate_process1, error = 1.0e-5)$value

# Sample a realization of the IPP 
loc <- rpoispp(lambda = Vectorize(intensity, c("x", "y")), win= window, nsim = 1, r = rho, grid = grid_points, covariate = covariate_process1)
loc <- rpoispp(lambda = eval.im(rho(covariate_process1)), win= window, nsim = 1, r = rho, grid = grid_points, covariate = covariate_process1)
plot(loc, add = TRUE, pch = 16, col = "white")
}

################################################# GGplot covariate process ###################################
plot_list_cov <-list()
for(i in 1:3){
covariate_process <- cov1_list[[i]]
loc <- loc2d_list_doublebump[[i]]

covariate_process_ppobj <- ppp(x = grid_points$x, y = grid_points$y,
                               marks = as.vector(t(covariate_process$v)), window = window)
covariate_process_plot <- Smooth(covariate_process_ppobj, kernel = "gaussian", dimyx = c(200,200), method = "pwl")


img_df <- data.frame(value = as.vector(covariate_process_plot$v))
img_df$x <- gridcentres(window, nx = 200, ny = 200)$y-0.5
img_df$y <- gridcentres(window, nx = 200, ny = 200)$x-0.5

points <- data.frame(x = loc$x-0.5,
                     y = loc$y-0.5)

plot_list_cov[[i]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
  scale_fill_gradientn(colors =colorRampPalette(coolwarm(11))(100), name = expression(Z[1](x))) +
  # scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  name = "Posterior \n width", 
  #                      limits = c(0, 6.1), 
  #                      breaks = c(0, 1, max(img_df$value), 3, 4, 5 ,6),  # Add max explicitly
  #                      labels = c("0","1",paste0("max = ", round(max(img_df$value),1)),"3","4", "5", "6")) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = expression(x[1]), y = expression(x[2])) #+
  # geom_point(data = points, aes(x = x, y = y), shape = 19, color="black", size = 1, alpha = 1) #+
# geom_polygon(data = triangles_df, aes(x = x, y = y, group = triangle_id), 
# fill = NA, color = "aquamarine2", linewidth = 0.1)
}
combined_plot <- wrap_plots(plot_list_cov, nrow = 1)
combined_plot

############################### posterior inference single obs #############################################
{
source(file = "pCN_MCMC_GP.R")
# discretization of the covariate space
discr <- seq(0,1, length.out = 300)
# posterior sampling with pCN MCMC
n_iter <- 20000
post_sample <- MCMC(n_iter, beta = 0.1, K = 300, loc, grid_points, covariate_process1, d_k = discr, 
                    squared_exponential_kernel, length_scale = 1, sigma_2 = 1, nu = 3/2)

# posterior mean and credible bands against truth and kernel estimators 
plot(rhohat(loc, covariate_process1), xlab = "Covariate values (with observations Rug)",
     ylab = "intensity", main = "Comarison of freq kernel estimator and BNP posterior inference",
     legend = FALSE, ylim = c(0,800))
lines(discr, 
      (colMeans(post_sample[seq(n_iter/2,n_iter),])), col = "blue", lwd = 2)
lines(discr, 
      rho(discr), 
      lwd = 2, col = "red")
legend("topleft", legend = c("kernel est", "truth", "post mean"), col = c("black", "red", "blue"), lty = 1, lwd= 1)
# rug(mapply(covariate_function, loc$x, loc$y))
# for(i in seq(1000,2000,2)){
#   lines(discr,
#         post_sample[i,], lwd = 0.1, col= rgb(0, 0, 0, alpha = 0.1))
# }
upper = apply(post_sample[seq(n_iter/2,n_iter),], 2, quantile, probs = 0.975)
lower = apply(post_sample[seq(n_iter/2,n_iter),], 2, quantile, probs = 0.025)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)
}
# adaptive sampler with sigmoid intensity function 
{
setwd("~/Research/IPP/")
source("pCN_adaptive_MCMC_GP.R")
discr <- seq(0,1, length.out = 300)
n_iter = 30000
adapt_MCMC(n_iter, beta = 0.08, K = length(discr), loc, grid_points, covariate_process2, d_k = discr, 
           matern_kernel, length_scale = 0.3, sigma_2 = 1, nu = 7/2, shape = 7, rate = 0.01)

lambda_post = read.table("Post_l_star.csv", sep = ",")
plot(t(lambda_post)[1:n_iter], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
Ls[Ls>10] <- 0
plot(t(Ls)[1:n_iter], type = "l")

gs2 = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs2 <- apply(gs2, 2, as.numeric)
post_mean2 <- colMeans(gs2[seq(4/5 * dim(gs2)[1], dim(gs2)[1], 1),],na.rm = TRUE)

plot(rhohat(loc, covariate_process1), xlab = "Covariate values (with observations Rug)",
     ylab = "intensity", main = "Comarison of freq kernel estimator and BNP posterior inference",
     legend = FALSE, ylim = c(0,800))
lines(discr, post_mean1, col = "blue", lwd = 2)
lines(discr, 
      rho(discr), 
      lwd = 2, col = "red")
legend("topleft", legend = c("kernel est", "truth", "post mean"), col = c("black", "red", "blue"), lty = 1, lwd= 1)# multiple observations 

upper = apply(gs1[seq(3/5 * n_iter,n_iter),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower = apply(gs1[seq(3/5 * n_iter,n_iter),], 2, quantile, probs = 0.025, na.rm = TRUE)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)


}
###################################### multiple observations ################################################
{
N <- 1e+04
cov_list <- list()
loc_list <- list()
dim = length(grid_points$x)
for(i in 1:N){
  # covariate_process <- as.vector(rnorm(n = dim, mean = 0, sd = 1) %*% Chol)
  # covariate_process <- pnorm((covariate_process - mean(covariate_process))/sd(covariate_process))
  # int_i <- im(matrix(covariate_process, sqrt(dim), sqrt(dim))[sqrt(dim):1,],
  #                           xrange = c(0,1), yrange = c(0,1))
  # cov_list[[i]] <- int_i
  int_i <- cov_list[[i]]
  loc_list[[i]] <- rpoispp(lambda = eval.im(rho(int_i)), win= window, nsim = 1)
  print(i)
}
}
# naive kernel counterpart
kernel_function <- function(loc_i, cov_i, M) {
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
L2_exp = matrix(NA,6,100)
L2_rel_exp = matrix(NA,6,100)
i = 1
discr <- seq(0,1, length.out = 200)
for(M in c(500)){
  for(exp in 1:10){
    multi_kernel <- Map(function(loc_i, cov_i){kernel_function(loc_i, cov_i)}, 
                        loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
                        cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
    arr <- array(unlist(multi_kernel), dim = c(200, 3, length(multi_kernel)))
    avg_matrix <- apply(arr, c(1, 2), mean, na.rm = TRUE)
    L2_exp[i,exp]<-(sqrt(sum((avg_matrix[,1] - rho(discr))^2)/length(discr)))
    L2_rel_exp[i,exp]<-(sqrt(sum((avg_matrix[,1] - rho(discr))^2)/(sum(rho(discr)^2))))
    print(paste(M, exp))
  }
  plot(discr, avg_matrix[,1], type = "l")
  lines(discr, 
        avg_matrix[,2], 
        lwd = 2, col = "red")
  lines(discr, 
        avg_matrix[,3], 
        lwd = 2, col = "red")
  lines(discr, 
        rho(discr), 
        lwd = 2, col = "blue")
  i <- i+1
}

funcs <- list(mean = mean, sd = sd)
result <- lapply(funcs, function(f) apply(L2_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))

result <- lapply(funcs, function(f) apply(L2_rel_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))

for(M in c(1000)){
  setwd("~/Research/IPP/")
  source(file = "pCN_MCMC_GP_mult.R")
  folder = paste("Adaptive_MCMC_post_multi/", M, sep = "")
  dir.create(paste("~/Research/IPP/", folder, sep = ""))
  setwd(folder)
  for(exp in 1){
    setwd("~/Research/IPP/")
    source(file = "pCN_MCMC_GP_mult.R")
    discr <- seq(0,1, length.out = 200)
    n_iter <- 20000
    adapt_multi_MCMC(paste0(M, "/", exp, sep = ""), n_iter, beta = 0.08, K = length(discr), 
                     loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], grid_points, 
                   cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], d_k = discr,
                   squared_exponential_kernel, alpha1 = 1, alpha2 = 1, 
                   sigma_2 = 1, nu = 3/2, exp_param = c(2,2),
                   shape = 1, rate = 2, "lambda")
  }
}   


lambda_post = read.table("Post_l_star.csv", sep = ",")
plot(t(lambda_post)[1:n_iter], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
plot(t(Ls)[5000:n_iter], type = "l")

L2_bay_exp = matrix(NA,5,100)
L2_bay_rel_exp = matrix(NA,5,100)
i = 1
discr <- seq(0,1, length.out = 200)

for(M in c(10)){
  setwd("~/Research/IPP/")
  folder = paste("Adaptive_MCMC_post_multi/", M, sep = "")
  setwd(folder)
  for(exp in 1:2){
    folder = paste("Adaptive_MCMC_post_multi/", M, "/", exp, sep = "")
    setwd(folder)
    gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
    gs <- apply(gs, 2, as.numeric)
    post_mean <- colMeans(gs[seq(n_iter/4,n_iter),],na.rm = TRUE)
    
    L2_bay_exp[i,exp]<- sqrt(sum(post_mean - rho(discr))^2)/length(discr)
    L2_bay_rel_exp[i,exp]<- sqrt(sum((post_mean - rho(discr))^2)/(sum(rho(discr)^2)))
  }
  plot(discr, post_mean, col = "blue", lwd = 2, type = "l", xlab = "Covariate values (with observations Rug)",
       ylab = "intensity", main = "Comarison of freq kernel estimator and BNP posterior inference",
       legend = FALSE)#, ylim = c(0,max(rho(discr))+2))
  lines(discr,
        rho(discr),
        lwd = 2, col = "red")
  legend("topleft", legend = c("kernel est", "truth", "post mean"), col = c("black", "red", "blue"), lty = 1, lwd= 1)# multiple observations
  
  upper = apply(gs[seq(8/10*n_iter,n_iter),], 2, quantile, probs = 0.95, na.rm = TRUE)
  lower = apply(gs[seq(8/10*n_iter,n_iter),], 2, quantile, probs = 0.05, na.rm = TRUE)
  polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)
  rug <- unlist(Map(function(c,l){return(c[l])} , cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
                    loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)]))
  rug(rug)
  i<-i+1
}

funcs <- list(mean = mean, sd = sd)
result <- lapply(funcs, function(f) apply(L2_bay_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))

result <- lapply(funcs, function(f) apply(L2_bay_rel_exp, 1, f, na.rm = TRUE))
print(paste(round(result$mean,3), "(", round(result$sd,3), ")&"))

plot_list <-list()
for(i in 1:3){
  M <- c(250,500,1000)[i]
  exp <- 1
  rug <- unlist(Map(function(c,l){return(c[l])} , cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
                    loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)]))
  multi_kernel <- Map(function(loc_i, cov_i){kernel_function(loc_i, cov_i)}, 
                      loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
                      cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  arr <- array(unlist(multi_kernel), dim = c(200, 3, length(multi_kernel)))
  avg_matrix <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  
  spline_fit_est <- smooth.spline(discr, avg_matrix[,1], spar = 0.8)
  spline_fit_upp <- smooth.spline(discr, avg_matrix[,2], spar = 0.8)
  spline_fit_low <- smooth.spline(discr, avg_matrix[,3], spar = 0.8)
  
 
  # ggplots 
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
    geom_rug(data = data.frame(x = rug, y = 0), sides = "b") #+ 
    geom_ribbon(data = data.frame(x = discr, 
                                  y2 = predict(spline_fit_upp, x = discr)$y, 
                                  y1 = pmax(0,predict(spline_fit_low, x = discr)$y), 
                                  y = predict(spline_fit_est, x = discr)$y), mapping = aes(x = x, ymin = y1, y = y, ymax = y2),
                            fill = "white", alpha = 0.2) + ylim(0,10) + 
    geom_rug(data = data.frame(x = rug, y = 0), sides = "b")
}
combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot

# compute likelihood of posterior draws 
{
exp = 1
M = 250
Z <-Map(function(cov, lo){cov[lo]}, 
        cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
        loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
index_dom <- lapply(cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
                    function(c){pmin(round(c$v*200), 200)})				   
index_loc <- lapply(Z, function(c){pmin(round(c*200), 200)})		

log_lik <- apply(as.data.frame(gs), 1, function(rho){
  rho_obs <- log(rho[unlist(index_loc)])
  I <- lapply(index_dom, function(i){sum(rho[i]-1, na.rm = TRUE)/2500})
  sum(rho_obs) - sum(unlist(I))
})
true_lik <- sum(log(rho(discr)[unlist(index_loc)])) - 
  sum(unlist(lapply(index_dom, function(i){sum(rho(discr)[i]-1, na.rm = TRUE)/2500})))
}
plot(2:20000,log_lik[2:20000], type = "l")
abline(h = true_lik, col  ="red")
############################## 2-DIMENSIONAL COVARIATE PROCESS ####################################
# covariate process
{
  # covariance matrix of the GP generating the covariate process, discretized over the grid
  Cov2 <- exponential_cov(distance_matrix, length_scale = 0.05, sigma_2 = 1) + 1e-10 * diag(1,length(grid_points$x))
  Chol2 <- Matrix::chol(Cov2)
  plot(im(Cov2[dim(Cov2)[1]:1,]))
  
  covariate_process2 <- im(matrix(rnorm(n = length(grid_points$x), mean = 0, sd = 1) %*% Chol2 + mu, 
                               sqrt(length(grid_points$x)), sqrt(length(grid_points$x)))[sqrt(length(grid_points$x)):1,],
                        xrange = c(0,1), yrange = c(0,1))

  covariate_process2 <- ((covariate_process2 - mean(covariate_process2))/sd(covariate_process2))
  # to unit interval
  covariate_process2 <- eval.im(pnorm(covariate_process2))
  
  covariate_process_ppobj2 <- ppp(x = grid_points$x, y = grid_points$y,
                                 marks = as.vector(t(covariate_process2$v)), window = window)
  covariate_process_plot2 <- Smooth(covariate_process_ppobj2, kernel = "gaussian", dimyx = c(200,200), method = "pwl")
  
  plot(covariate_process_plot, main = "First covariate process on unit square")
  plot(covariate_process_plot2, main = "Second covariate process on unit square")
}
#point pattern realization by intensity function from covariate space to [0, infty]
{
rho1 <- function(z1, z2){
  # bump function
  # 500 * dmsn(cbind(z1, z2), xi=c(0.4, 0.4), Omega = diag(0.03, nrow = 2), alpha = c(0,0), tau=1, dp=NULL, log=FALSE)
  14 * dmsn(cbind(z1, z2), xi=c(0.7, 0.2), Omega = diag(0.05, nrow = 2), alpha = c(3,-2), tau=1, dp=NULL, log=FALSE)
  # 90 * dmsn(cbind(z1, z2), xi=c(0.3, 0.3), Omega = diag(0.5, nrow = 2), alpha = c(-1, -1), tau=1, dp=NULL, log=FALSE)
}
rho2 <- function(z1, z2){
  # bump function
  # 400 * dmsn(cbind(z1, z2), xi=c(0.8, 0.9), Omega = diag(0.01, nrow = 2), alpha = c(0,0), tau=0, dp=NULL, log=FALSE)
  6 * dmsn(cbind(z1, z2), xi=c(0.3, 0.8), Omega = diag(0.03, nrow = 2), alpha = c(-1,-1), tau=0, dp=NULL, log=FALSE)
  # 10 * dmsn(cbind(z1, z2), xi=c(0.8, 0.8), Omega = diag(0.01, nrow = 2), alpha = c(0,0), tau=0, dp=NULL, log=FALSE)
}
rho2d <- function(z1, z2){
  # bump function
  pmax(0, rho1(z1,z2) + rho2(z1,z2))
}

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


exp_cov <- expand.grid(x = seq(0,1, length.out = 200), 
            y = seq(0,1, length.out = 200))

rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                      sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
               xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))

plot(rho_imag, axes=TRUE, main = "Rho function", xlab = "First covariate", 
     ylab = "Second covariate")


############################################# ggplot ###############################################
img_df <- data.frame(value = as.vector(rho_imag$v))
img_df$x <- gridcentres(window, nx = 200, ny = 200)$y
img_df$y <- gridcentres(window, nx = 200, ny = 200)$x

ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z)), limits = c(0, 17)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = expression(Z[1](x)), y = expression(Z[2](x)))



# functions returning intensity for given coordinates (x,y)
# by returning intensity of the closest grid point
intensity2d <- function(x, y, r, grid, covar1, covar2){
  index <- which.min(sqrt((grid$x - x)^2 + ((grid$y) - y)^2))
  X_index <- which.min(abs(covar1$xcol - grid$x[index]))
  Y_index <- which.min(abs(covar1$yrow - grid$y[index]))
  return(r(covar1[Y_index, X_index], covar2[Y_index, X_index]))
}

Vectorize(intensity2d, c("x", "y"))

loc <- rpoispp(lambda = im(matrix(rho2d(as.vector(covariate_process1$v), as.vector(covariate_process2$v)),50,50)),
               win= window, nsim = 1, 
               grid = grid_points, covar1 = covariate_process1, covar2 = covariate_process2, r = rho2d)
loc_complete <- do.call(superimpose, loc)
loc_complete <- unmark(loc_complete)
points(covariate_process1[loc], covariate_process2[loc],pch = "+", col = "white", cex = 0.5)
plot(rhohat(loc, covariate_process1), col = "black", main = "Marginal intensity first covariate", ylim = c(0,20), 
     ylab = "Marginal intensity", xlab = "First covariate value", lwd = 2, legend = FALSE)
lines(seq(0,1,length.out = 200), 
      colMeans(matrix(rho2d(rep(seq(0,1,length.out = 200), 200), rep(seq(0,1,length.out = 200), each = 200)), 
                      200 , byrow = TRUE)), col = "red", lwd = 2)
lines(seq(0,1,length.out = 200), 
      colMeans(rho2hat(loc_complete, covariate_process1, covariate_process2, method = "ratio", dimyx = 200)$v),
      col = "purple", lwd = 2)
points(covariate_process1[loc], covariate_process2[loc], pch = "+", col = "white", cex = 0.5)
plot(im(matrix(mapply(rho2d, covariate_process1$v, covariate_process2$v), nrow = 50), 
        xrange = c(0,1), yrange = c(0,1)), main = "Intensity on unit square")
plot(loc, add = TRUE, pch = "+", col = "white")

plot(covariate_process1, main = "First covaraite process on spatial domain with observatios")
plot(loc, add = TRUE, pch = "+", col = "white")
plot(covariate_process2, main = "First covaraite process on spatial domain with observatios")
plot(loc, add = TRUE, pch = "+", col = "white")


}
# marginals intensities
{
plot(rhohat(loc, covariate_process1), col = "black", main = "Marginal intensity first covariate", ylim = c(100,1300), 
     ylab = "Marginal intensity", xlab = "First covariate value", lwd = 2, legend = FALSE)
lines(seq(0,1,length.out = 200), 
      colMeans(matrix(rho2d(rep(seq(0,1,length.out = 200), 200), rep(seq(0,1,length.out = 200), each = 200)), 
                                                  200 , byrow = TRUE)), col = "red", lwd = 2)
lines(seq(0,1,length.out = 200), 
      colMeans(rho2hat(loc, covariate_process1, covariate_process2, method = "ratio", dimyx = 200)$v),
      col = "purple", lwd = 2)
lines(discr, post_mean1, col = "blue", lwd = 2, type = "l", xlab = "Covariate values (with observations Rug)",
     ylab = "intensity", main = "Comarison of freq kernel estimator and BNP posterior inference",
     legend = FALSE)#, ylim = c(0,max(rho(discr))+2))
# lines(seq(0,1,length.out = 50), 
#       colMeans(t(apply(gs[seq(3/5*dim(gs)[1],dim(gs)[1]),], 1, 
#                        function(x){colMeans(matrix(x, nrow = 50, byrow = TRUE), na.rm = TRUE)}))),
#       col = "forestgreen", lwd = 2)
lines(post_mean_m$yrow, 
      colMeans(post_mean_m$v, na.rm = TRUE),col = "forestgreen", lwd = 2)

upper = apply(gs1[seq(3/5 *dim(gs1)[1],dim(gs1)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower = apply(gs1[seq(3/5 *dim(gs1)[1],dim(gs1)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)

# marginal 2d on grid
# upper1 = apply(t(apply(gs[seq(3/5*dim(gs)[1],dim(gs)[1]),], 1, function(x){colMeans(matrix(x, nrow = 50, byrow = TRUE), na.rm = TRUE)})),
#                  2, quantile, probs = 0.975, na.rm = TRUE)
# lower1 = apply(t(apply(gs[seq(3/5*dim(gs)[1],dim(gs)[1]),], 1, function(x){colMeans(matrix(x, nrow = 50, byrow = TRUE), na.rm = TRUE)})),
#                2, quantile, probs = 0.025, na.rm = TRUE)
# polygon(c(discr2d$x, rev(discr2d$x)), c(upper1, rev(lower1)), col = rgb(0.13, 0.55, 0.13, alpha = 0.2), border = NA)

list_run <- apply(gs_m[seq(30000, 50000, 5),], 1, function(row){
  t(colMeans(Smooth(ppp(x = vertices[,1], y = vertices[,2],
                        marks = row, window = owin(poly_win)), 
                    kernel = "gaussian", dimyx = c(100,100), method = "pwl")$v))
})
marg_quant <- apply(t(list_run), 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
polygon(c(seq(0,1,length.out = 100), rev(seq(0,1,length.out = 100))), 
        c(marg_quant[2,], rev(marg_quant[1,])), col = rgb(0.13, 0.55, 0.13, alpha = 0.2), border = NA)

legend("topr", legend = c("univar kernel est","marginal kernel est", "truth", "univariate post mean", "marginal posterior"), col = c("black", "purple", "red", "blue", "forestgreen"), lty = 1, lwd= 1)# multiple observations 


plot(rhohat(loc, covariate_process2), col = "black", lwd = 2, 
      type = "l", main = "Marginal intensity second covariate", ylim = c(50,1300), 
      ylab = "Marginal intensity", xlab = "Second covariate value", legend = FALSE)
lines(seq(0,1,length.out = 100), 
      rowMeans(matrix(rho2d(rep(seq(0,1,length.out = 100), 100), rep(seq(0,1,length.out = 100), each = 100)),100 , byrow = TRUE)),, col = "red", lwd = 2)

lines(seq(0,1,length.out = 100), rowMeans(rho2hat(loc, covariate_process1, covariate_process2, method = "ratio", dimyx = 100)$v),
      col = "purple", lwd = 2)
lines(discr, post_mean2, col = "blue", lwd = 2, type = "l")#, ylim = c(0,max(rho(discr))+2))
# lines(seq(0,1,length.out = 50), 
#       colMeans(t(apply(gs[seq(3/5*dim(gs)[1],dim(gs)[1]),], 1, function(x){rowMeans(matrix(x, nrow = 50, byrow = TRUE))}))), col = "forestgreen", lwd = 2)
lines(post_mean_m$xcol, 
      rowMeans(post_mean_m$v, na.rm = TRUE),col = "forestgreen", lwd = 2)

legend("topr", legend = c("univar kernel est","marginal kernel est", "truth", "univariate post mean", "marginal posterior"), col = c("black", "purple", "red", "blue", "forestgreen"), lty = 1, lwd= 1)# multiple observations 


upper = apply(gs2[seq(8/10 * dim(gs2)[1],dim(gs2)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower = apply(gs2[seq(8/10 * dim(gs2)[1],dim(gs2)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0,0,1, alpha = 0.1), border = NA)

list_run2 <- apply(gs_m[seq(30000, 50000, 5),], 1, function(row){
  t(rowMeans(Smooth(ppp(x = vertices[,1], y = vertices[,2],
                        marks = row, window = owin(poly_win)), 
                    kernel = "gaussian", dimyx = c(100,100), method = "pwl")$v))
})
marg_quant2 <- apply(t(list_run2), 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
polygon(c(seq(0,1,length.out = 100), rev(seq(0,1,length.out = 100))), 
        c(marg_quant2[2,], rev(marg_quant2[1,])), col = rgb(0.13, 0.55, 0.13, alpha = 0.2), border = NA)


}
# posterior inference
{
source(file = "pCN_MCMC_2dGP.R")

discr <- data.frame(x = seq(0,1,length.out = 30), y = seq(0,1,length.out = 30))

n_iter = 10000
post_sample <- MCMC(n_iter, beta = 0.1, K = dim(discr)[1]^2, loc, grid_points, 
                   covariate_process1, covariate_process2, 
                   discr, matern32_kernel, length_scale = c(0.35, 0.35), sigma_2 = 1, nu = 3/2, lambda_star = 1000)
post_mean <- matrix(colMeans(post_sample[seq(n_iter/2,n_iter),]), nrow = 30, byrow = FALSE)
plot(im(post_mean, xrange = c(0,1), yrange = c(0,1)), main = "Posterior intensity mean on unit square")

rho_hat_2 <- rho2hat(loc, covariate_process1, covariate_process2, method = "ratio", dimyx = 100)

intensity2d_est <- function(x, y, rho, d, grid, Z1, Z2){
  Z <- covariate_function2(x, y, grid, Z1, Z2)
  return(interp2(d$x, d$y, rho, Z[1], Z[2], method = "constant"))
}
plot(as.im(Vectorize(intensity2d_est, c("x", "y")), rho = rho_hat_2$v, 
           d = discr, grid = grid_points, Z1 = covariate_process1, Z2 = covariate_process2, 
           W = window),
     main = "Posterior mean of the intensity function with observations")
plot(loc, add = TRUE, pch = "+", cex = .7, col = "white")

plot(rho_hat_2, zlim = range(c(0,1500)))
points(covariate_process1[loc], covariate_process2[loc], pch="+", cex = 0.7, col="white")
}
# 2-d adaptive 
{
  setwd("~/Research/IPP/")
  source(file = "pCN_adaptive_MCMC_2dGP_index.R")
  
  discr2d <- expand.grid(x = seq(0,1,length.out = 10), y = seq(0,1,length.out = 10))
    # data.frame(x = seq(0,1,length.out = 10), y = seq(0,1,length.out = 10))
  
  # posterior inference
  n_iter = 50000
  
  plot(seq(0,10000,0.01), dgamma(seq(0,10000,0.01), shape = 7, rate = 0.01), type = "l")
  
  setwd("~/Research/IPP/")
  MCMC("test_index_with_mesh1", n_iter, beta = 0.15, K = dim(discr2d)[1],
        loc, grid_points, 
        covariate_process1, covariate_process2, 
        vertex = discr2d, centroids = discr2d, 
        matern52_kernel, length_scale = c(0.35,0.35), sigma_2 = 1, nu = 3/2,
        shape = 5, rate = 0.001)
  
  lambda_post = read.table("Post_l_star.csv", sep = ",", nrows = n_iter)
  lambda_post[lambda_post>10000] <- 0 
  plot(t(lambda_post)[1:n_iter], type = "l")
  
  Ls = read.table("Post_ls.csv", sep = ",", nrows = n_iter)
  Ls[Ls>10] <- 0 
  plot(t(Ls)[1:n_iter], type = "l")
  
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE, nrows = n_iter)
  gs <- apply(gs, 2, as.numeric)
  post_mean <- matrix(colMeans(gs[seq(3/5 * n_iter,n_iter),],na.rm = TRUE), nrow = dim(discr2d)[1], byrow = TRUE)
  plot(im(post_mean, xrange = c(0,1), yrange = c(0,1)), main = "Posterior intensity mean on unit square", zlim = range(c(0,1500)))
  points(covariate_process1[loc], covariate_process2[loc], pch = "+", cex = .7, col = "white")
  
  post_lower <- matrix(apply(gs[seq(n_iter/2,n_iter),], 2, quantile, probs = 0.05), nrow = dim(discr2d)[1], byrow = TRUE)
  post_upper <- matrix(apply(gs[seq(n_iter/2,n_iter),], 2, quantile, probs = 0.95), nrow = dim(discr2d)[1], byrow = TRUE)
  
  post_uncert <- im(post_upper - post_lower, xrange = c(min(discr2d$x),max(discr2d$x)), yrange = c(min(discr2d$y),max(discr2d$y)))
  plot(post_uncert, 
       main = "Posterior bandwidth on unit square")
  
  
  rho_hat_2 <- rho2hat(loc, covariate_process1, covariate_process2, method = "ratio", dimyx = 200)
  plot(rho_hat_2, zlim = range(c(0,1500)))
  z <- mapply(covariate_function2, loc$x, loc$y, MoreArgs = list(grid_points, covariate_process1, covariate_process2))
  points(z[1,], z[2,], pch="+", cex = 0.7, col="white")
  
  intensity2d_est <- function(x, y, rho, d, grid, Z1, Z2){
    Z <- covariate_function2(x, y, grid, Z1, Z2)
    return(interp2(d$x, d$y, rho, Z[1], Z[2], method = "constant"))
  }
  lambda_est <- as.im(Vectorize(intensity2d_est, c("x", "y")), rho = post_mean, 
                   d = discr2d, #data.frame(x = rho_hat_2$yrow, y = rho_hat_2$xcol)
                   grid = grid_points, Z1 = covariate_process1, Z2 = covariate_process2, 
                   W = window)
  lambda_true <- as.im(Vectorize(intensity2d_est, c("x", "y")), rho = rho_imag$v, 
                       d = data.frame(x = seq(0,1, length.out = 50), 
                                      y = seq(0,1, length.out = 50)), 
                       grid = grid_points, Z1 = covariate_process1, Z2 = covariate_process2, 
                       W = window)
  
  plot(lambda_est, main = "Posterior mean of the intensity function with observations", zlim = range(c(0,1500)))
  plot(lambda_true, main = "Posterior mean of the intensity function with observations", zlim = range(c(0,1500)))
  plot(loc, add = TRUE, pch = "+", cex = .7, col = "white")
}
# multiple observations 2d adaptive
{
# all different covariate process 
{
N <- 2e+04
load("~/Research/IPP/Adaptive_MCMC_post_multi/res_100exp.Rdata")
cov1_list <- cov_list
cov2_list <- list()
loc2d_list_bumphole <- list()
dim = length(grid_points$x)
for(i in 1:N){
  # covariate_process2 <- as.vector(rnorm(n = dim, mean = 0, sd = 1) %*% Chol2)
  # covariate_process2 <- pnorm((covariate_process2 - mean(covariate_process2))/sd(covariate_process2))
  # int_i <- im(matrix(covariate_process2, sqrt(dim), sqrt(dim))[sqrt(dim):1,],
  #             xrange = c(0,1), yrange = c(0,1))
  # cov2_list[[i]] <- int_i
  
  loc2d_list_bumphole[[i]] <- rpoispp(lambda = im(matrix(rho2d(as.vector(cov1_list[[i]]$v), 
                                                                 as.vector(cov2_list[[i]]$v)),50,50),
                                                    xrange = c(0,1), yrange = c(0,1)),
                                        win= window, nsim = 1)
  print(i)
}
}
# multiple observations, on same covariate process, equivalent to infill asymptotics, but amenable to replicated observations 
{ 
  cov_list1_alleq <- rep(list(covariate_process1), 1500)
  cov_list2_alleq <- rep(list(covariate_process2), 1500)
  all_loc_inf2 <- rpoispp(lambda = Vectorize(intensity2d, c("x", "y")), 
                         win= window, nsim = 10, 
                         grid = grid_points, covar1 = covariate_process1, covar2 = covariate_process2, r = rho2d)
  loc_complete <- do.call(superimpose, all_loc_inf2)
  loc_complete <- unmark(loc_complete)
  plot(rho_imag, axes=TRUE, main = "Rho function", xlab = "First covariate", ylab = "Second covariate")
  points(covariate_process1[loc_complete], covariate_process2[loc_complete], pch = "+", col = "white")
  
  rho_hat_2 <- rho2hat(loc_complete, covariate_process1, covariate_process2, method = "ratio", dimyx = 50)
  plot(rho_hat_2/10, zlim = c(0,150), main = "Kernel estimator devided by N")
  
  setwd("~/Research/IPP/")
  source(file = "pCN_adaptive_MCMC_2dGP.R")
  discr <- data.frame(x = seq(0,1,length.out = 50), y = seq(0,1,length.out = 50))
  n_iter = 10000
  MCMC("_rep_same_cov", n_iter, beta = 0.08, K = dim(discr)[1]^2, loc_complete, grid_points, 
       covariate_process1, covariate_process2, 
       discr, matern52_kernel, length_scale = c(0.35,0.35), sigma_2 = 1, nu = 3/2,
       shape = 5, rate = 0.001)
  
  lambda_post = read.table("Post_l_star.csv", sep = ",", nrows = n_iter)
  lambda_post[lambda_post>10000] <- 0 
  plot(t(lambda_post)[1000:n_iter], type = "l")
  
  Ls = read.table("Post_ls.csv", sep = ",", nrows = n_iter)
  Ls[Ls>10] <- 0 
  plot(t(Ls)[1:n_iter], type = "l")
  
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE, nrows = n_iter)
  gs <- apply(gs, 2, as.numeric)
  post_mean <- matrix(colMeans(gs[seq(3/5 * n_iter,9183),],na.rm = TRUE), nrow = dim(discr)[1], byrow = FALSE)
  plot(im(post_mean, xrange = c(0,1), yrange = c(0,1)), main = "Posterior intensity mean/N on unit square")
  points(covariate_process1[loc_complete], covariate_process2[loc_complete], pch = "+", cex = .4, col = "white")
}

i = 0
############################################## posterior inference ##########################################
for(M in c(2, 12)){
i <- i+1
setwd("~/Research/IPP/")
source("Adaptive_pCN_MCMC_2dGP_mult.R")
discr <- data.frame(x = seq(0,1,length.out = 50), y = seq(0,1,length.out = 50))

n_iter = 10000

plot(seq(0,100,0.01), dgamma(seq(0,100,0.01), 3, 0.001), type = "l")
# profvis(
# Adapt2d_MCMC(paste("50x50_inflat_replicated_cov",M,sep = "_"), n_iter, beta = 0.1, K = dim(discr)[1]^2, loc2d_list_inflated[1:M], grid_points,
#              cov1_list[1:M], cov2_list[1:M],
#      discr, matern52_kernel, length_scale = 0.25, sigma_2 = 1, nu = 5/2,
#      shape = 3, rate = 0.001)
# )

# see multiple observations on same covariate process as single observations on multiple (all equal) covariate processes  
# profvis(
Adapt2d_MCMC(paste("50x50_inflate_same_cov_",M,sep = "_"), n_iter, beta = 0.1, K = dim(discr)[1]^2, all_loc_inf[1:M], grid_points,
             cov_list1_alleq[1:M], cov_list2_alleq[1:M],
             discr, matern52_kernel, length_scale = 0.3, sigma_2 = 1, nu = 5/2,
             shape = 6, rate = 0.001)
# )


lambda_post = read.table("Post_l_star.csv", sep = ",")
plot(t(lambda_post)[1:n_iter], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
plot(t(Ls)[1:n_iter], type = "l")

gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs <- apply(gs, 2, as.numeric)
post_mean <- matrix(colMeans(gs[(7/10 * n_iter):dim(gs)[1],],na.rm = TRUE), nrow = dim(discr2d)[1], byrow = TRUE)
plot(im(post_mean, xrange = c(0,1), yrange = c(0,1)), main = "Posterior intensity mean on unit square", zlim = range(c(0,1500)))
Map(function(c1,c2,l) points(c1[l], c2[l], ex = 0.7, col = "white", pch = "+"),
    cov1_list[1:M], cov2_list[1:M], loc2d_list_inflated[1:M])

exp_cov <- expand.grid(x = seq(0,1, length.out = 200), 
                       y = seq(0,1, length.out = 200))

rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                      sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
               xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))

# d2_l2_dist_replicatedcov <- data.frame(N = c(2,10,50,250), L2 = rep(NA,4))
d2_l2_dist_replicatedcov$L2[i] <- sum((post_mean/100 - rho_imag)^2)/length(exp_cov$x)
}
print(xtable(d2_l2_dist_replicatedcov), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = I)

############################################# ggplot ###############################################
img_df <- data.frame(value = as.vector(rho_imag$v))
img_df$x <- gridcentres(window, nx = 800, ny = 800)$y
img_df$y <- gridcentres(window, nx = 800, ny = 800)$x

default_colormap <- function(x) hcl(h = 240 - 240 * x, c = 100, l = 65)
colors <- default_colormap(seq(0, 1, length.out = 100))

 plot_list2[[6]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  name = expression(rho(z)), limits = c(0, 40)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = expression(Z[1](x)), y = expression(Z[2](x)))


}

########################################### average kernel ######################################################
kern2_est <- rho2hat(loc2d_list_slope[[1]], cov1_list[[1]], cov2_list[[1]], method = "ratio", dimyx = c(200,200))

L2d_exp = matrix(NA,5,100)
L2d_rel_exp = matrix(NA,5,100)
i = 1
for(M in c(2,10,50,250,1000)){
  for(exp in 1:10){
      kern_years <- lapply(c(((exp-1)*1000 + 1):((exp-1)*1000 + M)), function(i){rho2hat(loc2d_list_slope[[i]], 
                              cov1_list[[i]], cov2_list[[i]], method = "ratio", dimyx = c(200,200))$v})
      kern2_est$v <- apply(simplify2array(kern_years), c(1, 2), mean, na.rm = TRUE)
      L2_exp[i,exp]<-(sqrt(sum((kern2_est$v - rho_imag)^2)/length(exp_cov$x)))
      L2_rel_exp[i,exp]<-(sqrt(sum((kern2_est$v - rho_imag)^2)/sum(rho_imag^2)))
      print(c(M,exp))
  }
  plot(kern2_est)
  Map(function(c1,c2,l) points(c1[l], c2[l], ex = 0.7, col = "white", pch = "+"),
      cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
      cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
      loc2d_list_slope[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  i <- i+1
}

############################# triangular mesh construction and inference ######################################
{
  # convex hull of existing covariates 
  # x <- as.vector(covariate_process1$v) 
  # y <- as.vector(covariate_process2$v)
  # points <- cbind(x,y)
  # points <- cbind((x-min(elev))/(max(elev)- min(elev)),(y-min(slope))/(max(slope)- min(slope)))
  # conv_hull <- chull(points)
  # conv_hull <- rev(c(conv_hull, conv_hull[1]))
  # hull_points <- cbind(x, y)[conv_hull,]
  # hull_points <- cbind((y-min(slope))/(max(slope)-min(slope)), (x-min(elev))/(max(elev)-min(elev)))[conv_hull,]
  # poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
  
  # convex hull as whole domain 
  points <- matrix(c(0,0,0,1,1,1,1,0), ncol = 2, byrow = T)
  conv_hull <- chull(points)
  conv_hull <- rev(c(conv_hull, conv_hull[1]))
  poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
  
  plot(poly_win, main = "Covariate domain", xlab = "Slope", ylab = "Elevetion")
  lines(points, col = "forestgreen", lwd = 2)
  # points(hull_points, col = "black", cex = 1, pch = 16)
  # points(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "orange", cex = 1, pch = 16)
  
  # triangular tesseletion of the convex hull, independent of points 
  
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
  points(vertices, col = "black", pch = 16, cex = .4)
  vertices <- mesh$P
  extremes <- mesh$P[mesh$PB==1,]
  colnames(extremes) <- c("x", "y")
  colnames(vertices) <- c("x", "y")
  
  inside.owin(x = extremes[,1], y = extremes[,2], w = poly_win)
  
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
}
{
  setwd("~/Research/IPP/")
  source(file = "pCN_adaptive_MCMC_3mesh_GP.R")
  
  n_iter = 20000
  
  MCMC("_mesh3_ind_aniso", n_iter, beta = 0.2, K = dim(vertices)[1], 
       loc, grid_points, window,
       covariate_process1, covariate_process2, 
       as.data.frame(vertices), triangle_centroids,
       exponential_cov, alpha1=1, alpha2=2, sigma_2 = 1, nu = NA,
       exp_param = c(2,2), shape = 1, rate = 0.1, 
       link = "lambda", LS = "anisotropic")
  
  lambda_post = read.table("Post_l_star.csv", sep = ",")
  plot(t(lambda_post)[1:dim(lambda_post)[1]], type = "l")
  
  Ls = read.table("Post_ls.csv", sep = ",", fill = TRUE)
  plot(t(Ls)[1,1:dim(Ls)[1]], type = "l")
  
  theta = read.table("Post_theta.csv", sep = ",", fill = TRUE)
  plot(t(theta)[2,1:dim(theta)[1]], type = "l")
  
  gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs_m <- apply(gs_m, 2, as.numeric)
  pp <- ppp(x = vertices[,1], y = vertices[,2],
            marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
  post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(dim(mesh$P)[1], dim(mesh$P)[1]), method = "pwl")
  plot(post_mean_m, zlim = range(c(0,1000)), main = "Posterior for aggregated data")
  
  lines(hull_points, col = "black", lwd = 2)
  points(vertices, col = "black", cex = 1, pch = 16)
  points(covariate_process1[loc], covariate_process2[loc], pch = "+", col = "white", cex = .8)
  
  apply(mesh$T, 1, function(tri_idx) {
    tri_coords <- mesh$P[tri_idx, ]
    polygon(tri_coords, border = "black", lwd = 1)
  })
  
  pp_uncert <- ppp(x = vertices[,1], y = vertices[,2],
                   marks = apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.95) - 
                     apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.05), window = owin(poly_win))
  
  post_uncert <- Smooth(pp_uncert, kernel = "gaussian", dimyx = c(dim(mesh$P)[1], dim(mesh$P)[1]), method = "pwl")
  
  plot(post_uncert, 
       main = "Posterior bandwidth on unit square")
  
  intensity2d_est <- function(rho, Z1, Z2){
    return(matrix(interp2(rho$yrow, rho$xcol, rho$v, Z1, Z2, method = "constant"),byrow = FALSE, nrow = covariate_process1$dim[1]))
  }
  plot(as.im(intensity2d_est(rho = post_mean_m, 
                             Z1 = as.vector(covariate_process1$v), Z2 = as.vector(covariate_process2$v)), 
             W = window),
       main = "Posterior mean of the intensity function with observations")
  plot(loc, add = TRUE, pch = "+", cex = .7, col = "white")
  
}
# multiple observations 2d with mesh
{
  for(M in c(2,10,50)){
    setwd("~/Research/IPP/")
    source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
    folder = paste("2dAdaptive_MCMC_post_mesh_replicated_obs/repexp/", M, sep = "")
    dir.create(paste("~/Research/IPP/", folder, sep = ""))
    setwd(folder)
    for(exp in 1:50){
      setwd("~/Research/IPP/")
      source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
      i <- i+1
      n_iter  = 20000

      MCMC(paste("2dAdaptive_MCMC_post_mesh_replicated_obs/repexp/",M, "/", exp,sep = ""), n_iter, beta = 1.5/M,
         K = dim(vertices)[1], loc2d_list_doublebump[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
         grid_points, window,
         cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
         cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
         as.data.frame(vertices), triangle_centroids,
         exponential_cov, alpha1=1, alpha2=2, sigma_2 = 1, nu = NA,
         exp_param = c(2,2), shape = 1, rate = 1, 
         link = "lambda", LS = "anisotropic")
    }
  }
  
  plot_list2 <- list()
  i = 0
  for(M in c(2, 12, 60, 300, 900)){
    i = i+1
    folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_replicated_obs/","fine_mesh_",M,sep = "") 
    setwd(folder)
    
    lambda_post = read.table("Post_l_star.csv", sep = ",")
    plot(t(lambda_post)[1:n_iter], type = "l")
    
    Ls = read.table("Post_ls.csv", sep = ",")
    plot(t(Ls)[2,1:dim(Ls)[1]], type = "l")
    
    expon = read.table("Post_theta.csv", sep = ",")
    plot(t(expon)[2,1:dim(expon)[1]], type = "l")
    mean(t(expon)[1,1:n_iter])
    
    gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
    gs_m <- apply(gs_m, 2, as.numeric)
    pp <- ppp(x = vertices[,1], y = vertices[,2],
              marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
    post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
    plot(post_mean_m, zlim = range(c(0,45)))
    Map(function(c1,c2,l) points(c1[l], c2[l], cex = 0.7, col = "white", pch = "+"),
        cov1_list[1:M], cov2_list[1:M], loc2d_list_new[1:M])
    
    pp_uncert <- ppp(x = vertices[,1], y = vertices[,2],
                     marks = apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.95) - 
                       apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.05), window = owin(poly_win))
    
    post_uncert <- Smooth(pp_uncert, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
    
    
    dd = post_mean_m$dim[1]
    
    exp_cov <- expand.grid(x = seq(1.084202e-19, 1, length.out = 800), 
                           y = seq(1.084202e-19, 1, length.out = 800))
    
    rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                          sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
                   xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
    
    # d2_l2_dist_slope <- data.frame(N = c(2,12,60,300,1500), L2 = rep(NA,5))
    
    d2_l2_dist_slope$L2[i] <-
      sqrt(sum((post_mean_m - rho_imag)^2)/length(exp_cov$x))
    d2_l2_dist_slope$L2ratio[i] <-
      sqrt(sum((post_mean_m - rho_imag)^2)/length(exp_cov$x))/sqrt(sum((rho_imag)^2)/length(exp_cov$x))
  
    ####################################### ggplot ###############################################
    img_df <- data.frame(value = as.vector(post_mean_m$v))
    img_df$x <- exp_cov$y
    img_df$y <- exp_cov$x
    
    
    points <- data.frame(x = unlist(Map(function(cov, lo){cov[lo]}, cov1_list[1:M], loc2d_list_new[1:M])),
                         y = unlist(Map(function(cov, lo){cov[lo]}, cov2_list[1:M], loc2d_list_new[1:M])))
    
    plot_list2[[i]] <- ggplot() +
      geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
      # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
      scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"), name = expression(rho(z)), limits = c(0, 35)) +
      # scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  name = "Posterior \n width", 
      #                      limits = c(0, 6.1), 
      #                      breaks = c(0, 1, max(img_df$value), 3, 4, 5 ,6),  # Add max explicitly
      #                      labels = c("0","1",paste0("max = ", round(max(img_df$value),1)),"3","4", "5", "6")) +
      coord_fixed() +  # Keeps aspect ratio correct
      theme_minimal() + 
      labs(x = expression(Z[1](x)), y = expression(Z[2](x))) #+
      #geom_point(data = points, aes(x = x, y = y), shape = 3, color="white", size = 0.3, alpha = 0.5) #+
      # geom_polygon(data = triangles_df, aes(x = x, y = y, group = triangle_id), 
                   # fill = NA, color = "aquamarine2", linewidth = 0.1)
    
      
  }
  print(xtable(t(d2_l2_dist_slope)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = I)
  
  # d2_l2_dist_hole <- data.frame(N = c(2,12,60,300,1500), L2 = rep(NA,5))
  
  i = 0
  for(M in c(2, 12, 60, 300, 1500)){
    i = i+1
  
  kern_years <- lapply(c(1:M), function(i){rho2hat(loc2d_list_slope[[i]], cov1_list[[i]], cov2_list[[i]], method = "ratio", dimyx = c(200,200))$v})
  kern2_est$v <- apply(simplify2array(kern_years), c(1, 2), mean, na.rm = TRUE)
  plot(kern2_est)
  Map(function(c1,c2,l) points(c1[l], c2[l], ex = 0.7, col = "white", pch = "+"),
      cov1_list[1:M], cov2_list[1:M], loc2d_list_slope[1:M])
  
  d2_l2_dist_slope$L2[i] <-
    sqrt(sum((kern2_est$v - rho_imag)^2)/length(exp_cov$x))
  d2_l2_dist_slope$L2ratio[i] <-
    sqrt(sum((kern2_est$v - rho_imag)^2)/length(exp_cov$x))/sqrt(sum((rho_imag)^2)/length(exp_cov$x))
}

# plot_list2 <- list()
combined_plot <- wrap_plots(plot_list2, nrow = 3)
combined_plot
}

####################### diagonal plot ################# 
i <- 1
for(M in c(50, 250, 1000)){
  
  setwd("~/Research/IPP/")
  folder = paste("2dAdaptive_MCMC_post_mesh_slope_repl/", M, sep = "")
  setwd(folder)
  
  exp = 1
  
  folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_slope_repl/", M, "/", exp, sep = "")
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

  x = 1:800
  y = x #801 - x
  
  pp_upper <- ppp(x = vertices[,1], y = vertices[,2],
            marks = 1.05 * apply(gs_m[seq(2/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.9999, na.rm = TRUE), window = owin(poly_win))
  
  upper = Smooth(pp_upper, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
  
  pp_lower <- ppp(x = vertices[,1], y = vertices[,2],
                  marks = apply(gs_m[seq(2/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.01, na.rm = TRUE), window = owin(poly_win))
  
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
    geom_ribbon(data = data.frame(x = seq(0,1,length.out = 800), y2 = upper$v[cbind(x,y)], y1 = lower$v[cbind(x,y)], y=post_mean_m$v[cbind(x,y)]), mapping = aes(x = x, ymin = y1, y = y, ymax = y2),
                fill = "dodgerblue3", alpha = 0.2)
  # print(plot_M + ggtitle(paste(M,exp)))
  plot_list[[3]] <- plot_M
  i <- i+1
}

combined_plot <- wrap_plots(plot_list[1:3], nrow = 1)
combined_plot
