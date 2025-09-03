################## LOAD LIBRIARIES #############################################

library(locfit)
library(spatstat)
library(spatstat.explore)
library(spatstat.geom)
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
library(geometry)
library(alphahull)
library(RTriangle)
library(sp)
library(concaveman)
library(raster)
library(magrittr)
library(pracma)
library(patchwork)
library(viridis)
library(scico)
library(RColorBrewer)
library(pals)
library(fields)

################## DATA LOADING AND PREPRATION #################################

# import of replicated covariate processes and observations 
load("~/ontario_data_june_orig_scale.RData")

# map covariates to unit interval with ecdf transform  
rescaling <- function(x, list){
  pixels <- unlist(lapply(list, function(img) as.vector(img$v)))
  F_z <- ecdf(pixels) 
  return(eval.im(F_z(x)))
}

# WHEN CHANGING TO OTHER REGION SIMPLY FIND AND REPLACE "ontario" WITH REGION'S NAME

# subset to years from 2004 to 2022
ontario_temp_june <-ontario_temp_june[c(1:19)]
ontario_precip_june <-ontario_precip_june[c(1:19)]
ontario_wind_june <-ontario_wind_june[c(1:19)]

ontario_ppp_june <-ontario_ppp_june[c(1:19)]

cov_1_list <- lapply(ontario_temp_june, rescaling,ontario_temp_june) 
cov_2_list <- lapply(ontario_precip_june, rescaling,ontario_precip_june) 
cov_3_list <- lapply(ontario_wind_june, rescaling,ontario_wind_june) 

# define observation window and discretization grid  
ontario_mask = im(!is.na(ontario_precip_june$`2004`$v), xcol =ontario_precip_june$`2004`$xcol, yrow =ontario_precip_june$`2004`$yrow)
ontario_dom = owin(ontario_mask$xrange,ontario_mask$yrange)
grid_points <- expand.grid(x = cov_2_list$`2010`$xcol, y = cov_2_list$`2010`$yrow)

# plot covariate maps for selection of years: Figure 4

# define color palettes
my_col_map <- colorRampPalette(c("darkblue", "deeppink3", "darkgoldenrod1"))
temp_palette <- colorRampPalette(coolwarm(11))
spectral_colors <- brewer.pal(11, "Spectral")
wind_palette <- colorRampPalette(rev(spectral_colors))

# comment out line 68 and "[ontario_dom, drop = FALSE], add = TRUE," at line 69 to remove canada's shadow
par(mfrow = c(1,3),        
    mar = c(0, 1, 0, 2),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))  
plot(mask_canada, col = gray(c(1, 0.9)), main = "", ribbon = FALSE)
plot(ontario_temp_june$`2013`[ontario_dom, drop = FALSE], add = TRUE, main = "",
     col = temp_palette(256), ribbon = TRUE, ribargs = list(cex.axis = 1.5, cex.lab = 2)) #
points(ontario_ppp_june$`2013`$x,ontario_ppp_june$`2013`$y, pch = 18, col = "black")
plot(ontario_precip_june$`2013`, main = "", col = brewer.blues(256), ribargs = list(cex.axis = 1.5, cex.lab = 2))
points(ontario_ppp_june$`2013`$x,ontario_ppp_june$`2013`$y, pch = 18, col = "black")
plot(ontario_wind_june$`2013`, main = "", col = wind_palette(256), ribargs = list(cex.axis = 1.5, cex.lab = 2))
points(ontario_ppp_june$`2013`$x,ontario_ppp_june$`2013`$y, pch = 18, col = "black")

# pixel to geographical coordinates mapping 
# ONTARIO
pixel_to_lon <- function(x) scales::rescale(x, to = c(-95.16, -74.34), from = c(378.9401, 629.2343))
pixel_to_lat <- function(y) scales::rescale(y, to = c(41.68, 56.85), from = c(29.52635, 263.85878))
y_breaks <- seq(29.52635, 263.85878, length.out = 5)
x_breaks <- seq(378.9401, 629.2343, length.out = 5)
# ALBERTA
pixel_to_lon <- function(x) scales::rescale(x, to = c(-120, -110), from = c(134.66, 253.79))
pixel_to_lat <- function(y) scales::rescale(y, to = c(49, 60), from = c(148.7, 355.98))
y_breaks <- seq(148.7, 355.98, length.out = 5)
x_breaks <- seq(134.66, 253.79, length.out = 5)
# British columbia
pixel_to_lon <- function(x) scales::rescale(x, to = c(-139.03, -114.03), from = c(32.538, 172.7))
pixel_to_lat <- function(y) scales::rescale(y, to = c(48.18, 60), from = c(162.72, 438.09))
y_breaks <- seq(162.72, 438.09, length.out = 5)
x_breaks <- seq(32.538, 172.7, length.out = 3)
# Saskatchewan 
pixel_to_lon <- function(x) scales::rescale(x, to = c(-110, -101.3), from = c(212.75, 319.87))
pixel_to_lat <- function(y) scales::rescale(y, to = c(49, 60), from = c(132.68, 330.95))
y_breaks <- seq(132.68, 330.95, length.out = 5)
x_breaks <- seq(212.75, 319.87, length.out = 3)


################## ONE DIMENSIONAL PRELIMINARY ANALYSIS: Section 4.1 ###########

# average kernel estimator: change cov_i_list to switch between TEMP == 1, PRECIP == 2, WIND == 3
kern_list <- Map(function(ppp,cov){as.function(rhohat(ppp, cov, n = 200, method = "ratio"))},
                ontario_ppp_june, cov_1_list)
# match number of evaluation points (n = 300)
xgrid <- seq(0, 1, length.out = 300)
# weighted average 
fmat <- t(sapply(kern_list, function(f) f(xgrid))) * unlist(lapply(ontario_ppp_june, function(x){x$n}))
kern_est <- colSums(fmat/sum(unlist(lapply(ontario_ppp_june, function(x){x$n}))), na.rm = TRUE)
# smoothing step 
spline_fit <- smooth.spline(xgrid, kern_est, spar = 0.8)

# Bayesian analysis
setwd("~/")
source(file = "pCN_MCMC_GP_mult.R")

discr = xgrid
n_iter <- 20000

adapt_multi_MCMC("ontario_temp_1d", n_iter, beta = 0.1, K = length(discr), 
                 ontario_ppp_june, grid_points,  
                 cov_1_list, d_k = discr, squared_exponential_kernel, 
                 alpha1 = 1, alpha2 = 1, sigma_2 = 1, nu = 3/2, 
                 exp_param = c(2,2), shape = 1, rate = 2, "lambda")

gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs <- apply(gs, 2, as.numeric)
post_mean_tot <- colMeans(gs[seq(8/10 * dim(gs)[1],10/10*dim(gs)[1]),],na.rm = TRUE)
upper_tot = apply(gs[seq(8/10 * dim(gs)[1],10/10*dim(gs)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower_tot = apply(gs[seq(8/10 * dim(gs)[1],10/10*dim(gs)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)

# rescaling to original covariates for visualization purposes
# change ontario_cov_june to switch between temp == 1, precip == 2, wind == 3
all_pixels <- unlist(lapply(ontario_temp_june, function(img) as.vector(img$v)))
xx_grid = quantile(all_pixels, probs = discr, type = 1, na.rm = TRUE)

rug <- unlist(Map(function(c,l){return(c[l])},ontario_temp_june,ontario_ppp_june))

# Figure 5: change indexing of plot_1d list when changing covariate 
plot_1d <- list()
plot_1d[[1]] <- ggplot(data.frame(x=xx_grid, y=post_mean_tot[1:length(xgrid)]), aes(x = x, y = y)) + 
  geom_line(linetype = "solid", color = "gray30", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black"),
    axis.text  = element_text(size = 11)
  ) + labs(x = "Temperature", y = expression(rho(z))) + 
  geom_line(data = data.frame(x=xx_grid, y=predict(spline_fit, x = discr)$y), mapping = aes(x = x, y = y),
            color = "gray20", size = 1, linetype = "dashed") +
  geom_ribbon(data = data.frame(x = xx_grid, y2 = upper_tot[1:length(xgrid)], 
                                             y1 = lower_tot[1:length(xgrid)], 
                                             y=post_mean_tot[1:length(xgrid)]), 
              mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "gray20", alpha = 0.2) + 
  geom_rug(data = data.frame(x = rug, y = 0), sides = "b")

# Figure 6: change indexing of plot_lambda_1d list when changing year, 

plot_lambda_1d <- list()
# intensity on the spatial domain 
intensity_est <- function(rho, mask, Z1){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  Z1 <- Z1[!is.na(Z1)]
  res = matrix(NA, nrow = mask$dim[1], ncol = mask$dim[2])
  entries <- approx(seq(0,1,length.out = length(rho)), rho, xout = Z1)$y
  res[mask$v] <- entries
  return(res)
}

img_df2013 <- data.frame(value = as.vector(t((im(intensity_est(rho = post_mean_tot,ontario_mask, cov_1_list$`2013`),
                                                 xrange = cov_1_list$`2013`$xrange, yrange = cov_1_list$`2013`$yrange)$v)[dim(ontario_mask$v)[1]:1,])))
img_df2013$x <- grid_points$x
img_df2013$y <- rev(grid_points$y)
points <- data.frame(x =ontario_ppp_june$`2013`$x, y =ontario_ppp_june$`2013`$y)

bbox <- data.frame(
  xmin = min(img_df2013$x),
  xmax = max(img_df2013$x),
  ymin = min(img_df2013$y),
  ymax = max(img_df2013$y)
)

plot_lambda_1d[[1]] <- 
  ggplot(img_df2013, aes(x = x, y = y, fill = value)) +
  geom_raster() +  
  scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x))) +
  coord_fixed() + theme_minimal() + labs(x = "Longitude", y = "Latitude")+
  geom_rect(data = bbox,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            ,fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) + 
  scale_x_continuous(
    breaks = x_breaks,
    labels = round(pixel_to_lon(x_breaks), 2),
    name = "Longitude") +
  scale_y_continuous(
    breaks = y_breaks,
    labels = round(pixel_to_lat(y_breaks), 2),
    name = "Latitude") +
  geom_point(data = points, aes(x = x, y = y), shape = 16, color="turquoise1", size = 0.9, inherit.aes = FALSE) 

combined_plot <- wrap_plots(plot_1d2, nrow = 1)
combined_plot

################## TWO DIMENSIONAL ANALYSIS: Section D1 ########################

# triangular mesh on unit square 
points <- matrix(c(0,0,0,1,1,1,1,0), ncol = 2, byrow = T)
conv_hull <- chull(points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
tri_tess <- RTriangle::pslg(P = cbind(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y))
mesh <- RTriangle::triangulate(tri_tess, a = 0.0012, j = TRUE)  #Y = TRUE, D = TRUE # smaller 'a' = finer mesh
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
vertices <- mesh$P
extremes <- mesh$P[mesh$PB==1,]
colnames(extremes) <- c("x", "y")
colnames(vertices) <- c("x", "y")
# domain restriction to existing covatiarte values
x <- unlist(lapply(cov_1_list, function(img) as.vector(img$v)))
y <- unlist(lapply(cov_2_list, function(img) as.vector(img$v)))
all_points <- cbind(x,y)[complete.cases(cbind(x,y)),]
conv_hull <- chull(all_points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
poly_win <- owin(poly = list(x = all_points[,1][conv_hull], y = all_points[,2][conv_hull]))

# average kernel estimator
kern2_list <- Map(function(ppp,cov1, cov2){as.function(rho2hat(ppp, cov1, cov2,  dimyx = c(400,400), method = "ratio"))},
                 ontario_ppp_june, cov_1_list, cov_2_list)

fmat <- t(sapply(kern2_list, function(f) f(vertices[,1], vertices[,2]))) * unlist(lapply(ontario_ppp_june, function(x){x$n}))
pp_kern <- ppp(x = vertices[,1], y = vertices[,2],
               marks = colSums(fmat/sum(unlist(lapply(ontario_ppp_june, function(x){x$n}))),na.rm = TRUE), window = owin(c(0,1),c(0,1)))
kern_tot <- Smooth(pp_kern, kernel = "gaussian", dimyx = c(800,800), method = "pwl", pad = TRUE, edge = FALSE)
kern_tot$xcol <- seq(kern_tot$yrange[1], kern_tot$yrange[2], length.out = 800)
kern_tot$yrow <- seq(kern_tot$xrange[1], kern_tot$xrange[2], length.out = 800)


# Bayesian analysis 
setwd("~/")
source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
MCMC("2dAdaptive_MCMC_multi_ontario", 
     n_iter, beta = 0.1, K = dim(vertices)[1],
     ontario_ppp_june, 
     grid_points,ontario_dom,
     cov_1_list, cov_2_list,
     as.data.frame(vertices), triangle_centroids,
     exponential_cov, alpha1=1, alpha2=1, sigma_2 = 1, nu = NA,
     exp_param = c(2,2), shape = 1, rate = 2, 
     link = "lambda", LS = "anisotropic")


gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs_m <- apply(gs_m, 2, as.numeric)
pp_tot <- ppp(x = vertices[,1], y = vertices[,2],
              marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(c(0,1), c(0,1)))
post_mean_tot <- Smooth(pp_tot, kernel = "quartic",  dimyx = c(800,800), method = "pwl", pad = TRUE, edge = FALSE)
post_mean_tot$xcol <- seq(post_mean_tot$yrange[1], post_mean_tot$yrange[2], length.out = 800)
post_mean_tot$yrow <- seq(post_mean_tot$xrange[1], post_mean_tot$xrange[2], length.out = 800)

#Figure 19
gg_2D_plot <- list()

all_pixelsX <- unlist(lapply(ontario_temp_june, function(img) as.vector(img$v)))
all_pixelsY <- unlist(lapply(ontario_precip_june, function(img) as.vector(img$v)))

img_df <- data.frame(value = (as.vector(post_mean_tot$v))) 
exp_cov <- expand.grid(x = seq(0, 1, length.out = 800), 
                       y = seq(0, 1, length.out = 800))

img_df$x <- exp_cov$y
img_df$y <- exp_cov$x

points <- data.frame(x = unlist(Map(function(c, l){c[l]}, cov_1_list,ontario_ppp_june)),
                     y = unlist(Map(function(c, l){c[l]}, cov_2_list,ontario_ppp_june)))
img_df$value2 <- as.vector(kern_tot$v) 

df_bdry <- data.frame(
  x = poly_win$bdry[[1]]$x,
  y = poly_win$bdry[[1]]$y
)

{
gg_2D_plot[[2]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z)), limits = c(0, 0.0071)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = "Temperature", y = "Precipitation") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 

gg_2D_plot[[3]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value2)) +  
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + coord_fixed(ratio = 1) +
  labs(x = "Temperature", y = "Precipitation") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, 
                     probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, 
                     probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 

gg_2D_plot[[1]] <- ggplot() +
  theme_minimal() + coord_fixed(ratio = 1) +
  labs(x = "Temperature", y = "Precipitation") +
  geom_path(data = data.frame(x = c(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$x[1]), 
                              y = c(poly_win$bdry[[1]]$y, poly_win$bdry[[1]]$y[1])), 
            aes(x = x, y = y), color = "grey50", linewidth = 0.5) +
  geom_point(data = points, aes(x = x, y = y), shape = 19, size = 1, color="black") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, 
                     probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, 
                     probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2))
}

# intensity on the spatial domain, find and substitute "13" with other years for other spatial estimates 
intensity2d_est <- function(rho, mask, Z1, Z2){
  # browser()
  Z1 = Z1[!is.na(Z1)]
  Z2 = Z2[!is.na(Z2)]
  Z1 = Z1[!is.na(Z1)]
  Z2 = Z2[!is.na(Z2)]
  res = matrix(NA, nrow = mask$dim[1], ncol = mask$dim[2])
  entries <- interp2(rho$xcol, rho$yrow, rho$v, Z1, Z2, method = "nearest")
  res[mask$v] <- entries
  return(res)
}

est <- data.frame(mean21 = as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,ontario_mask,
                  Z1 = as.vector(cov_1_list$`2021`$v), Z2 = as.vector(cov_2_list$`2021`$v)), 
                  W =ontario_dom)$v)[dim(ontario_mask)[1]:1,])))
est$x <- grid_points$x
est$y <- rev(grid_points$y)
points21 <- data.frame(x =ontario_ppp_june$`2021`$x, y =ontario_ppp_june$`2021`$y)

ggplot(est, aes(x = x, y = y, fill = mean21)) +
  geom_raster() +  
  scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x))) +
  coord_fixed() + theme_minimal() + labs(x = "Longitude", y = "Latitude")+
  geom_rect(data = bbox,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            ,fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) + 
  scale_x_continuous(
    breaks = x_breaks,
    labels = round(pixel_to_lon(x_breaks), 2),
    name = "Longitude") +
  scale_y_continuous(
    breaks = y_breaks,
    labels = round(pixel_to_lat(y_breaks), 2),
    name = "Latitude") +
  geom_point(data = points21, aes(x = x, y = y), shape = 16, size = 1, color="turquoise1", inherit.aes = FALSE) 

################## THREE DIMENSIONAL ANALYSIS: Section 4.2 #####################

# mesh on unit cube using C++ program tetgen
{
cube <- as.matrix(expand.grid(rep(list(c(0, 1)), 3))[,1:3])
setwd(paste("~/","3d_mesh", sep = ""))
write_poly_file <- function(points, file = "covariate_domain.poly") {
  library(geometry)
  # Compute convex hull faces
  faces <- geometry::convhulln(points)
  
  # Open connection
  con <- file(file, open = "wt")
  
  ## Part 1: Point list
  cat(nrow(points), "3 0 0\n", file = con)
  for (i in 1:nrow(points)) {
    cat(i, points[i, ], "\n", file = con)
  }
  
  # Face list 
  cat(nrow(faces), "0\n", file = con)
  for (i in 1:nrow(faces)) {
    cat("1\n", file = con)                        
    cat(ncol(faces), faces[i, ], "\n", file = con) 
  }
  
  # Holes (none)
  cat("0\n", file = con)
  
  # Regions (none)
  cat("0\n", file = con)
  
  close(con)
  message("Written .poly file: ", file)
}
write_poly_file(cube)
run_tetgen <- function(poly_file, flags = "-pq1.2a0.0004") {
  cmd <- sprintf("./tetgen %s %s", flags, poly_file)
  system(cmd)
}
run_tetgen("covariate_domain.poly",flags = "-pq1.2a0.0004")
nodes <- read.table(paste("covariate_domain", "1", "node", sep = "."), skip = 1)[, 2:4]
tets <- read.table(paste("covariate_domain", "1", "ele", sep = "."), skip = 1)[, 2:5]
}

# Bayesian analysis 
setwd("~/")
source(file = "pCN_mbGP_3Dmesh_multi.R")
n_iter = 40000

MCMC("_ontario", 
     n_iter, beta = 0.1, K = dim(nodes)[1],
     ontario_ppp_june,
     ontario_dom,
     cov_1_list, cov_2_list, cov_3_list,
     as.data.frame(nodes),
     exponential_cov, alpha1 = 1, alpha2 = 1, sigma_2 = 1, nu = NA,
     exp_param = c(2,2), shape = 1, rate = 2, 
     "lambda", "anisotropic") 

gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs_m <- apply(gs_m, 2, as.numeric)
rho_est <- colMeans(gs_m[seq(7/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE)

intensity3d_est <- function(rho, dim1, dim2, mask, cov1, cov2, cov3, vert, year){
  Z1 = as.vector(cov1[[year]]$v)
  Z2 = as.vector(cov2[[year]]$v)
  Z3 = as.vector(cov3[[year]]$v)
  Z1 = Z1[!is.na(Z1)]
  Z2 = Z2[!is.na(Z2)]
  Z3 = Z3[!is.na(Z3)]
  res = matrix(NA, nrow = dim1, ncol = dim2)
  entries <- rho[unlist(apply(cbind(as.vector(Z1),as.vector(Z2), as.vector(Z3)),
             1, function(w){which.min(rowSums((vert - matrix(w, nrow = dim(vert)[1], ncol = 3, byrow = TRUE))^2))}))]
  res[(mask == TRUE)$v] <- entries
  return(res)
}

# Figure 8: change plot3d[[i]] index for every year and substitue "2013" with alternatives 
plot3d <- list()

img_df <- data.frame(value = as.vector(t((blur(as.im(intensity3d_est(rho = rho_est,ontario_mask$dim[1], 
                                                                    ontario_mask$dim[2],ontario_mask, 
                                                                     cov1 = cov_1_list, 
                                                                     cov2 = cov_2_list, 
                                                                     cov3 = cov_3_list,
                                                                     vert = nodes, year = 2013-2003), 
                                                     W =ontario_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE)$v)[ontario_mask$dim[1]:1,])))
img_df$x <- grid_points$x
img_df$y <- rev(grid_points$y)

bbox <- data.frame(
  xmin = min(img_df$x),
  xmax = max(img_df$x),
  ymin = min(img_df$y),
  ymax = max(img_df$y)
)

points2013 <- data.frame(x =ontario_ppp_june$`2013`$x, y =ontario_ppp_june$`2013`$y)

plot3d[[1]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  scale_fill_gradientn(na.value = "white", 
                   colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                   name = expression(lambda(x))) + coord_fixed() + theme_minimal() + 
  labs(x = "Longitude", y = "Latitude")+
  geom_rect(data = bbox,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            ,fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) + 
  geom_point(data = points2013, aes(x = x, y = y), shape = 19, 
             color="turquoise1", size = 0.7, inherit.aes = FALSE) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = round(pixel_to_lon(x_breaks), 2),
    name = "Longitude") +
  scale_y_continuous(
    breaks = y_breaks,
    labels = round(pixel_to_lat(y_breaks), 2),
    name = "Latitude") 
  
################## BIVARIATE MARGINAL PLOTS: Figure 7 ##########################

cloud_3d <- data.frame(x = nodes[,1], y = nodes[,2], z = nodes[,3], f = rho_est)

model <- loess(f ~ x + y + z, data = cloud_3d, span = 0.2)  # span controls smoothing

z0 = quantile(unlist(lapply(cov_3_list, function(img) as.vector(img$v))), prob = 0.05, na.rm = TRUE)
z1 = quantile(unlist(lapply(cov_3_list, function(img) as.vector(img$v))), prob = 0.5, na.rm = TRUE)
z2 = quantile(unlist(lapply(cov_3_list, function(img) as.vector(img$v))), prob = 0.95, na.rm = TRUE)

x_grid <- rep(seq(0, 1, length.out = 200), each = 200)
y_grid <- rep(seq(0, 1, length.out = 200),200)
pred_df0 <- data.frame(x = x_grid, y = y_grid, z = rep(z0,40000))
pred_df1 <- data.frame(x = x_grid, y = y_grid, z = rep(z1,40000))
pred_df2 <- data.frame(x = x_grid, y = y_grid, z = rep(z2,40000))
q1_slice <- predict(model, newdata = pred_df0)
q2_slice <- predict(model, newdata = pred_df1)
q3_slice <- predict(model, newdata = pred_df2)

all_pixelsZ <-  unlist(lapply(ontario_wind_june, function(img) as.vector(img$v)))
round(quantile(all_pixelsZ, probs = c(c(0, 0.05, 0.5, 0.95, 1)), type = 1, na.rm = TRUE),2)

plot_marginal <- list()

all_pixelsX <- unlist(lapply(ontario_temp_june, function(img) as.vector(img$v)))
all_pixelsY <- unlist(lapply(ontario_precip_june, function(img) as.vector(img$v)))

img_df <- data.frame(value = q1_slice) 

img_df$x <- x_grid
img_df$y <- y_grid

points <- data.frame(x = unlist(Map(function(c, l){c[l]}, cov_1_list,ontario_ppp_june)),
                     y = unlist(Map(function(c, l){c[l]}, cov_2_list,ontario_ppp_june)))
img_df$value2 <- q2_slice
img_df$value3 <- q3_slice

# change index of plot_marginal and fill = value consistently: 1<->value, 2<->value2, 3<->value3
plot_marginal[[3]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value3)) + 
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) + #, limits = c(0, 0.0068)
  coord_fixed() + theme_minimal() + labs(x = "Temperature", y = "Precipitation") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 

combined_plot <- wrap_plots(plot_marginal, nrow = 1)
combined_plot







  
