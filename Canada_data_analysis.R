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

# import of replicated covariate processes and observations 
{
setwd(paste("~/Research/IPP/","Replicated_data_canada_fires", sep = ""))
load("~/Research/IPP/Replicated_data_canada_fires/data_june.RData")

rescaling <- function(x,list){
  all_pixels <- unlist(lapply(list, function(img) as.vector(img$v)))
  return((x - min(all_pixels,na.rm = TRUE))/(max(all_pixels,na.rm = TRUE) - min(all_pixels,na.rm = TRUE)))
}
#max(all_pixels,na.rm = TRUE)

rescaling <- function(x, list){
  all_pixels <- unlist(lapply(list, function(img) as.vector(img$v)))
  return(eval.im(pnorm((x - mean(all_pixels,na.rm = TRUE))/sd(all_pixels, na.rm = TRUE), sd = 1.5)))
}

rescaling <- function(x, list){
  pixels <- unlist(lapply(list, function(img) as.vector(img$v)))
  F_z <- ecdf(pixels) 
  return(eval.im(F_z(x)))
}

saska_temp_june <-saska_temp_june[c(1:19)]
saska_precip_june <-saska_precip_june[c(1:19)]
saska_wind_june <-saska_wind_june[c(1:19)]

saska_ppp_june <-saska_ppp_june[c(1:19)]

cov_1_list <- lapply(saska_temp_june, rescaling,saska_temp_june) #saska_temp_june #
cov_2_list <- lapply(saska_precip_june, rescaling,saska_precip_june) #saska_precip_june #
cov_3_list <- lapply(saska_wind_june, rescaling,saska_wind_june) #saska_wind_june #

pixels <- unlist(lapply(cov_1_list, function(img) as.vector(img$v)))
hist(pixels)
rug(unlist(Map(function(c,l){return(c[l])} , cov_1_list,saska_ppp_june)))

Map(function(pix, dat, i){hist(pix, freq = FALSE, main = 2003+i); rug(pix[dat]); hist(pix[dat], add = TRUE, freq = FALSE, col = rgb(0, 0, 1, 0.5))}, 
    cov_1_list,saska_ppp_june, seq_along(saska_ppp_june))

my_col_map <- colorRampPalette(c("darkblue", "deeppink3", "darkgoldenrod1"))

temp_palette <- colorRampPalette(coolwarm(11))
spectral_colors <- brewer.pal(11, "Spectral")
wind_palette <- colorRampPalette(rev(spectral_colors))
#brewer.blues

par(mfrow = c(1,3),        
    mar = c(0, 1, 0, 2),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))  
plot(mask_canada, col = gray(c(1, 0.9)), main = "", ribbon = FALSE)
plot(saska_temp_june$`2012`[saska_dom, drop = FALSE], add = TRUE, main = "",
     col = temp_palette(256), ribbon = TRUE, ribargs = list(cex.axis = 1.5, cex.lab = 2)) #[saska_dom, drop = FALSE], add = TRUE,
points(saska_ppp_june$`2012`$x,saska_ppp_june$`2012`$y, pch = 18, col = "black")
plot(saska_precip_june$`2012`, main = "", col = brewer.blues(256), ribargs = list(cex.axis = 1.5, cex.lab = 2))
points(saska_ppp_june$`2012`$x,saska_ppp_june$`2012`$y, pch = 18, col = "black")
plot(saska_wind_june$`2012`, main = "", col = wind_palette(256), ribargs = list(cex.axis = 1.5, cex.lab = 2))
points(saska_ppp_june$`2012`$x,saska_ppp_june$`2012`$y, pch = 18, col = "black")


saska_mask = im(!is.na(saska_precip_june$`2004`$v), xcol =saska_precip_june$`2004`$xcol, yrow =saska_precip_june$`2004`$yrow)
saska_dom = owin(saska_mask$xrange,saska_mask$yrange)

par(mfrow = c(1,3), 
    mar = c(0, 1, 0, 2),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))
plot(saska_mask, col = gray(c(1, 0.9)), main = "", ribbon = FALSE)
points(saska_ppp_june$`2006`$x,saska_ppp_june$`2006`$y, pch = 18)
plot(saska_mask, col = gray(c(1, 0.9)), main = "", ribbon = FALSE)
points(saska_ppp_june$`2013`$x,saska_ppp_june$`2013`$y, pch = 18)
plot(saska_mask, col = gray(c(1, 0.9)), main = "", ribbon = FALSE)
points(saska_ppp_june$`2013`$x,saska_ppp_june$`2013`$y, pch = 18)

grid_points <- expand.grid(x = cov_2_list$`2010`$xcol, y = cov_2_list$`2010`$yrow)

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
}
# mesh on the unit cube 
{
cube <- as.matrix(expand.grid(rep(list(c(0, 1)), 3))[,1:3])
setwd(paste("~/Research/IPP/","3d_mesh", sep = ""))
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
  
  ## Part 2: Face list (facets)
  cat(nrow(faces), "0\n", file = con)
  for (i in 1:nrow(faces)) {
    cat("1\n", file = con)                        # One polygon per facet
    cat(ncol(faces), faces[i, ], "\n", file = con) # Triangle with vertex indices
  }
  
  ## Part 3: Holes (none)
  cat("0\n", file = con)
  
  ## Part 4: Regions (none)
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
# posterior inference
{
  setwd("~/Research/IPP/")
  source(file = "pCN_mbGP_3Dmesh_multi.R")
  # unlink(paste("~/Research/IPP/", "3d_mesh_canada_fires_replicated_prec", sep = ""), recursive = TRUE)
  
  n_iter = 40000
  
  MCMC("_saska_ecdf_natural", n_iter, beta = 0.1, 
       K = dim(nodes)[1],saska_ppp_june,saska_dom,
       cov_1_list, cov_2_list, cov_3_list,
       as.data.frame(nodes),
       exponential_cov, alpha1 = 1, alpha2 = 2, sigma_2 = 1, nu = NA,
       exp_param = c(2,2), shape = 2, rate = 1, 
       "lambda", "anisotropic") #,Index3d_over_dom_prec, Index3d_over_loc_prec
  
  lambda_post = read.table("Post_l_star.csv", sep = ",")
  plot(t(lambda_post)[1:dim(lambda_post)[1]], type = "l")
  
  Ls = read.table("Post_ls.csv", sep = ",")
  plot(t(Ls)[3,1:dim(Ls)[1]], type = "l")
  
  expon = read.table("Post_theta.csv", sep = ",")
  plot(t(expon)[1,1:dim(expon)[1]], type = "l")
  
  gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs_m <- apply(gs_m, 2, as.numeric)
  rho_est <- colMeans(gs_m[seq(7/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE)
  
  intensity3d_est <- function(rho, dim1, dim2, mask, cov1, cov2, cov3, vert, year){
    # browser()
    Z1 = as.vector(cov1[[year]]$v)
    Z2 = as.vector(cov2[[year]]$v)
    Z3 = as.vector(cov3[[year]]$v)
    Z1 = Z1[!is.na(Z1)]
    Z2 = Z2[!is.na(Z2)]
    Z3 = Z3[!is.na(Z3)]
    res = matrix(NA, nrow = dim1, ncol = dim2)
    entries <- #rho[Index3d_over_dom_prec[[year]]]
      rho[unlist(apply(cbind(as.vector(Z1),as.vector(Z2), as.vector(Z3)),
                                 1, function(w){which.min(rowSums((vert - matrix(w, nrow = dim(vert)[1], ncol = 3, byrow = TRUE))^2))}))]
    res[(mask == TRUE)$v] <- entries
    return(res)
  }
  for(aa in 1:1){
  plot(blur(as.im(intensity3d_est(rho = rho_est,saska_mask$dim[1], 
                            saska_mask$dim[2],saska_mask, 
                             cov1 = cov_1_list, 
                             cov2 = cov_2_list, 
                             cov3 = cov_3_list,
                             vert = nodes, year = aa), 
             W =saska_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE),
       main = paste("Posterior mean of the intensity function with observations", 2003+aa))
  plot(saska_ppp_june[[aa]], add = TRUE, pch = "+", cex = .7, col = "green")
  }
  
  #### ggplots for 3 covariates ###
  
  img_df <- data.frame(value = as.vector(t((blur(as.im(intensity3d_est(rho = rho_est,saska_mask$dim[1], 
                                                                      saska_mask$dim[2],saska_mask, 
                                                                       cov1 = cov_1_list, 
                                                                       cov2 = cov_2_list, 
                                                                       cov3 = cov_3_list,
                                                                       vert = nodes, year = 2013-2003), 
                                                       W =saska_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE)$v)[saska_mask$dim[1]:1,])))
  img_df$x <- grid_points$x
  img_df$y <- rev(grid_points$y)
  
  
  bbox <- data.frame(
    xmin = min(img_df$x),
    xmax = max(img_df$x),
    ymin = min(img_df$y),
    ymax = max(img_df$y)
  )
  
  img_df$value2 <- as.vector(t((blur(as.im(intensity3d_est(rho = rho_est,saska_mask$dim[1], 
                                                                      saska_mask$dim[2],saska_mask, 
                                                                       cov1 = cov_1_list, 
                                                                       cov2 = cov_2_list, 
                                                                       cov3 = cov_3_list,
                                                                       vert = nodes, year = 2013-2003), 
                                                       W =saska_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE)$v)[saska_mask$dim[1]:1,]))
  
  img_df$value3 <- as.vector(t((blur(as.im(intensity3d_est(rho = rho_est,saska_mask$dim[1], 
                                                          saska_mask$dim[2],saska_mask, 
                                                           cov1 = cov_1_list, 
                                                           cov2 = cov_2_list, 
                                                           cov3 = cov_3_list,
                                                           vert = nodes, year = 2013-2003), 
                                           W =saska_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE)$v)[saska_mask$dim[1]:1,]))

  # img_df$value4 <- as.vector(t((blur(as.im(intensity3d_est(rho = rho_est,saska_mask$dim[1], 
  #                                                         saska_mask$dim[2],saska_mask, 
  #                                                          cov1 = cov_1_list, 
  #                                                          cov2 = cov_2_list, 
  #                                                          cov3 = cov_3_list,
  #                                                          vert = nodes, year = 2022-2003), 
  #                                          W =saska_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE)$v)[saska_mask$dim[1]:1,]))
  
  
  points1 <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  points2 <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  points3 <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  points4 <- data.frame(x =saska_ppp_june$`2022`$x, y =saska_ppp_june$`2022`$y)
  
  plot3d <- list()
  
  for(i in 1:19){
    
    img_df <- data.frame(value = as.vector(t((blur(as.im(intensity3d_est(rho = rho_est,saska_mask$dim[1], 
                                                                        saska_mask$dim[2],saska_mask, 
                                                                         cov1 = cov_1_list, 
                                                                         cov2 = cov_2_list, 
                                                                         cov3 = cov_3_list,
                                                                         vert = nodes, year = i), 
                                                         W =saska_dom),sigma = 2, NA, bleed=FALSE, normalise = TRUE)$v)[saska_mask$dim[1]:1,])))
    img_df$x <- grid_points$x
    img_df$y <- rev(grid_points$y)
    
    
    bbox <- data.frame(
      xmin = min(img_df$x),
      xmax = max(img_df$x),
      ymin = min(img_df$y),
      ymax = max(img_df$y)
    )
    points1 <- data.frame(x =saska_ppp_june[[i]]$x, y =saska_ppp_june[[i]]$y)
    
  plot3d[[i]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x))) + # , limits = c(0, 0.011)
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points1, aes(x = x, y = y), shape = 19, color="turquoise1", size = 0.7, inherit.aes = FALSE) 
  }
  plot3d[[2]] <- ggplot(img_df, aes(x = x, y = y, fill = value2)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x))) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points2, aes(x = x, y = y), shape = 19, color="turquoise1", inherit.aes = FALSE) 
  
  
  plot3d[[3]] <- ggplot(img_df, aes(x = x, y = y, fill = value3)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x))) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points3, aes(x = x, y = y), shape = 19, color="turquoise1", inherit.aes = FALSE) 
  
  plot3d[[4]] <- ggplot(img_df, aes(x = x, y = y, fill = value4)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x))) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points4, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  
  combined_plot <- wrap_plots(plot3d[c(3,5,7,9,12,13,15,18,19)], nrow = 3)
  combined_plot
  
  ################################## marginal plots ##################################
  dom <- seq(0,1,length.out = 20)
  m3 <- sapply(nodes[,1], function(x){which.min((x - dom)^2)})
  f3 <- tapply(rho_est, m3, mean)
  plot(dom, f3, main = "Polynomial Regression")
  model <- lm(f3 ~ poly(dom, 3, raw = TRUE))
  p3 <- predict(model)
  lines(dom, p3, col = "blue", lwd = 2)
  
  plot_marginal <- list()
  
  rug1 <- unlist(Map(function(c,l){return(c[l])} , cov_1_list,saska_ppp_june))
  rug2 <- unlist(Map(function(c,l){return(c[l])} , cov_2_list,saska_ppp_june))
  rug3 <- unlist(Map(function(c,l){return(c[l])} , cov_3_list,saska_ppp_june))


  # ggplots for univariate intensity function 
  plot_marginal[[1]] <- ggplot(data.frame(x=xpred, y=p1), aes(x = x, y = y)) +
    geom_line(color = "gray30", size = 1) +
    theme_minimal(base_family = "sans") +  # Clean white background
    theme(
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      text = element_text(color = "black")
    ) + labs(x = expression(Z[1] * (x)), y = expression(rho[1] * (z))) + ylim(0.001, 0.004) + 
    geom_rug(data = data.frame(x = rug1, y = 0), sides = "b")
  
  plot_marginal[[2]] <- ggplot(data.frame(x=dom, y=p2), aes(x = x, y = y)) +
    geom_line(color = "gray30", size = 1) +
    theme_minimal(base_family = "sans") +  # Clean white background
    theme(
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      text = element_text(color = "black")
    ) + labs(x = expression(Z[2] * (x)), y = expression(rho[2] * (z))) + ylim(0.001, 0.004) + 
    geom_rug(data = data.frame(x = rug2, y = 0), sides = "b")
  
  plot_marginal[[3]] <- ggplot(data.frame(x=dom, y=p3), aes(x = x, y = y)) +
    geom_line(color = "gray30", size = 1) +
    theme_minimal(base_family = "sans") +  # Clean white background
    theme(
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      text = element_text(color = "black")
    ) + labs(x = expression(Z[3] * (x)), y = expression(rho[3] * (z))) + ylim(0.001, 0.004) + 
    geom_rug(data = data.frame(x = rug3, y = 0), sides = "b")
  
  combined_plot <- wrap_plots(plot_marginal, nrow = 1)
  combined_plot
}  

################################ 1 dimensional slicing ##########################

cloud_3d <- data.frame(x = nodes[,1], y = nodes[,2], z = nodes[,3], f = rho_est)

model <- loess(f ~ x + y + z, data = cloud_3d, span = 0.2)  # span controls smoothing

x_pred <- rep(0,100)
for(y0 in seq(0,1,0.01)){
  for(z0 in seq (0,1,0.01)){
    y0 = median(unlist(lapply(cov_2_list, function(img) as.vector(img$v))), na.rm = TRUE)
    z0 = median(unlist(lapply(cov_3_list, function(img) as.vector(img$v))), na.rm = TRUE)
    x_grid <- seq(min(cloud_3d$x), max(cloud_3d$x), length.out = 100)
    pred_df <- data.frame(x = x_grid, y = y0, z = z0)
    x_pred <- x_pred + predict(model, newdata = pred_df)
  }
}
# Plot conditional curve
all_pixels <- unlist(lapply(saska_temp_june, function(img) as.vector(img$v)))
x_grid_orig = quantile(all_pixels, probs = x_grid, type = 1, na.rm = TRUE)

plot(x_grid_orig, x_pred, type = "l", col = "blue", main = "f(x | y ≈ y0, z ≈ z0)")
rug_x <- unlist(Map(function(c,l){return(c[l])} ,saska_temp_june,saska_ppp_june))
rug(rug_x, ticksize=0.03, side=1, lwd=0.5)


y_pred <- rep(0,100)
for(x0 in seq(0,1,0.01)){
  for(z0 in seq (0,1,0.01)){
    x0 = median(unlist(lapply(cov_1_list, function(img) as.vector(img$v))), na.rm = TRUE)
    z0 = median(unlist(lapply(cov_3_list, function(img) as.vector(img$v))), na.rm = TRUE)
    y_grid <- seq(min(cloud_3d$x), max(cloud_3d$x), length.out = 100)
    pred_df <- data.frame(x = x0, y = y_grid, z = z0)
    y_pred <- y_pred + predict(model, newdata = pred_df)
  }
}
# Plot conditional curve
all_pixels <- unlist(lapply(saska_precip_june, function(img) as.vector(img$v)))
y_grid_orig = quantile(all_pixels, probs = y_grid, type = 1, na.rm = TRUE)

plot(y_grid_orig, y_pred, type = "l", col = "blue", main = "f(y | x ≈ x0, z ≈ z0)")
rug_y <- unlist(Map(function(c,l){return(c[l])} ,saska_precip_june,saska_ppp_june))
rug(rug_y, ticksize=0.03, side=1, lwd=0.5)


z_pred <- rep(0,100)
for(y0 in seq(0,1,0.01)){
  for(x0 in seq (0,1,0.01)){
    x0 = median(unlist(lapply(cov_1_list, function(img) as.vector(img$v))), na.rm = TRUE)
    y0 = median(unlist(lapply(cov_2_list, function(img) as.vector(img$v))), na.rm = TRUE)
    z_grid <- seq(min(cloud_3d$x), max(cloud_3d$x), length.out = 100)
    pred_df <- data.frame(x = x0, y = y0, z = z_grid)
    z_pred <- z_pred + predict(model, newdata = pred_df)
  }
}
# Plot conditional curve
all_pixels <- unlist(lapply(saska_wind_june, function(img) as.vector(img$v)))
z_grid_orig = quantile(all_pixels, probs = z_grid, type = 1, na.rm = TRUE)

plot(z_grid_orig, z_pred, type = "l", col = "blue", main = "f(z | x ≈ x0, y ≈ y0)")
rug_z <- unlist(Map(function(c,l){return(c[l])} ,saska_wind_june,saska_ppp_june))
rug(rug_z, ticksize=0.03, side=1, lwd=0.5)

## ggplot

plot_marginal <- list()

rug1 <- unlist(Map(function(c,l){return(c[l])} , cov_1_list,saska_ppp_june))
rug2 <- unlist(Map(function(c,l){return(c[l])} , cov_2_list,saska_ppp_june))
rug3 <- unlist(Map(function(c,l){return(c[l])} , cov_3_list,saska_ppp_june))


# ggplots for univariate intensity function 
plot_marginal[[1]] <- ggplot(data.frame(x=x_grid_orig, y=x_pred), aes(x = x, y = y)) +
  geom_line(color = "gray30", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) + labs(x = "Temperatures", y = expression(rho[1] * (z))) + ylim(0, 0.003) + 
  geom_rug(data = data.frame(x = rug_x, y = 0), sides = "b")

plot_marginal[[2]] <- ggplot(data.frame(x=y_grid_orig, y=y_pred), aes(x = x, y = y)) +
  geom_line(color = "gray30", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) + labs(x = "Precipitation", y = expression(rho[2] * (z))) + ylim(0, 0.005) + 
  geom_rug(data = data.frame(x = rug_y, y = 0), sides = "b")

plot_marginal[[3]] <- ggplot(data.frame(x=z_grid_orig, y=z_pred), aes(x = x, y = y)) +
  geom_line(color = "gray30", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) + labs(x = "Wind speed", y = expression(rho[3] * (z))) + ylim(0, 0.0015) + 
  geom_rug(data = data.frame(x = rug_z, y = 0), sides = "b")

combined_plot <- wrap_plots(plot_marginal, nrow = 1)
combined_plot

################################ 2 dimensional slicing ##########################

cloud_3d <- data.frame(x = nodes[,1], y = nodes[,2], z = nodes[,3], f = rho_est)

model <- loess(f ~ x + y + z, data = cloud_3d, span = 0.2)  # span controls smoothing

z0 = quantile(unlist(lapply(cov_3_list, function(img) as.vector(img$v))), prob = 0.05, na.rm = TRUE)
x_grid <- rep(seq(0, 1, length.out = 200), each = 200)
y_grid <- rep(seq(0, 1, length.out = 200),200)
pred_df <- data.frame(x = x_grid, y = y_grid, z = rep(z0,40000))
q1_slice <- predict(model, newdata = pred_df)

plot(im(matrix(q1_slice,200)))

## ggplot

all_pixelsZ <-  unlist(lapply(saska_wind_june, function(img) as.vector(img$v)))
round(quantile(all_pixelsZ, probs = c(c(0, 0.05, 0.5, 0.95, 1)), type = 1, na.rm = TRUE),2)

plot_marginal <- list()

all_pixelsX <- unlist(lapply(saska_temp_june, function(img) as.vector(img$v)))
all_pixelsY <- unlist(lapply(saska_precip_june, function(img) as.vector(img$v)))

img_df <- data.frame(value = q1_slice) 

img_df$x <- x_grid
img_df$y <- y_grid

points <- data.frame(x = unlist(Map(function(c, l){c[l]}, cov_1_list,saska_ppp_june)),
                     y = unlist(Map(function(c, l){c[l]}, cov_2_list,saska_ppp_june)))
img_df$value2 <- q2_slice
img_df$value3 <- q3_slice


plot_marginal[[1]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) + #, limits = c(0, 0.0068)
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = "Temperature", y = "Precipitation") +
  #geom_point(data = points, aes(x = x, y = y), shape = 19, size = 0.5, color="turquoise1")+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 

plot_marginal[[2]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value2)) +  # or geom_tile() if you prefer
  scale_fill_gradientn(na.value = "darkgoldenrod1", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) + #,limits = c(0, 0.0068)
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + coord_fixed(ratio = 1) +
  labs(x = "Temperature", y = "Precipitation") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 


plot_marginal[[3]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value3)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
  scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) + # , limits = c(0, 0.0068)
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + coord_fixed(ratio = 1) +
  labs(x = "Temperature", y = "Precipitation") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 


combined_plot <- wrap_plots(plot_marginal, nrow = 1)
combined_plot


################################ 1D on temperature ################################

setwd("~/Research/IPP")
source("pCN_adaptive_MCMC_GP.R")
# unlink(paste("~/Research/IPP/", "Adaptive_MCMC_post_saska_2013", sep = ""), recursive = TRUE)

discr <- seq(0, 1, length.out = 300)
grid_points <- expand.grid(x =saska_temp_june$`2013`$xcol, y =saska_temp_june$`2013`$yrow)

n_iter = 30000

pixels <- as.vector(saska_temp_june$`2013`$v)
F_z <- ecdf(pixels)
x =saska_temp_june$`2013`

adapt_MCMC("saska_2013_NATecdf", n_iter, beta = 0.5, 
           K = length(discr),saska_ppp_june$`2013`, grid_points, 
           eval.im(F_z(x)), discr, squared_exponential_kernel, 
           alpha1 = 1, alpha2 = 3, sigma_2 = 1, nu = 3/2, 
           exp_param = c(2,2), shape = 1, rate = 1, "lambda")


lambda_post = read.table("Post_l_star.csv", sep = ",")
# lambda_post[lambda_post>10] <- 0
plot(t(lambda_post)[1:n_iter], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
plot(t(Ls)[1:n_iter], type = "l")

expon = read.table("Post_theta.csv", sep = ",")
plot(t(expon)[1:n_iter], type = "l")

gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs <- apply(gs, 2, as.numeric)
post_mean_dom <- colMeans(gs[seq(8/10 * dim(gs)[1],dim(gs)[1]),],na.rm = TRUE)

rho_hat <- rhohat(saska_ppp_june$`2013`, eval.im(F_z(x)), n = 300)

discr = quantile(x, prob = discr, type = 1)

plot(discr, rho_hat$rho,
     xlab = "Covariate values (with observations Rug)",
     ylab = "Intensity", main = "Comparison of freq kernel estimator and BNP posterior inference",
     legend = FALSE, type = "l")
lines(discr, post_mean_dom, col = "red", lwd = 2)
rug(saska_temp_june$`2013`[saska_ppp_june$`2013`], ticksize = 0.03, side = 1)

upper = apply(gs[seq(n_iter/2,n_iter),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower = apply(gs[seq(n_iter/2,n_iter),], 2, quantile, probs = 0.025, na.rm = TRUE)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(1, 0, 0, alpha = 0.2), border = NA)

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

plot(im(intensity_est(rho = rho_hat$rho,saska_mask, eval.im(F_z(x))),
        xrange =saska_temp_june$`2013`$xrange, yrange =saska_temp_june$`2013`$yrange),
     main = "Posterior mean of the intensity function with observations",
     zlim = c(0,0.006))
plot(saska_ppp_june$`2013`, add = TRUE, 
       cex = 0.5, col = "green", pch = "+")

plot(im(intensity_est(rho = post_mean_dom,saska_mask, eval.im(F_z(x))),
        xrange =saska_temp_june$`2013`$xrange, yrange =saska_temp_june$`2013`$yrange),
     main = "Posterior mean of the intensity function with observations",
     zlim = c(0,0.006))
points(saska_ppp_june$`2013`$x, 
      saska_ppp_june$`2013`$y, 
       col = "green", pch = "+")

############################### GGplots for single obs ######################################
rug <-saska_temp_june$`2013`[saska_ppp_june$`2013`]

plot_saska <- list()

# ggplots for univariate intensity function 
plot_saska[[1]] <- ggplot(data.frame(x=discr, y=rho_hat$rho), aes(x = x, y = y)) +
  geom_line(linetype = "dashed", color = "gray30", size = 1) +
  geom_line(data = data.frame(x=discr, y=rho_hat$ave), aes(x = x, y = y),
            linetype = "dotted", color = "gray20", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) + labs(x = expression(z^(2013)), y = expression(rho(z))) + 
  geom_line(data = data.frame(x=discr, y=post_mean_dom), mapping = aes(x = x, y = y), 
            color = "gray20", size = 1) + 
  geom_ribbon(data = data.frame(x = discr, y2 = upper, y1 = lower, y=post_mean_dom), mapping = aes(x = x, ymin = y1, y = post_mean_dom, ymax = y2),
              fill = "gray20", alpha = 0.2) + ylim(0, 0.01) + 
  geom_rug(data = data.frame(x = rug, y = 0), sides = "b")

#  gpglots for two dimensional spatial intensity function
img_df <- data.frame(value = as.vector(t((im(intensity_est(rho = rho_hat$rho,saska_mask, eval.im(F_z(x))),
                                          xrange =saska_temp_june$`2013`$xrange, yrange =saska_temp_june$`2013`$yrange)$v)[246:1,])))
img_df$x <- grid_points$x
img_df$y <- rev(grid_points$y)

img_df2 <- data.frame(value = as.vector(t((im(intensity_est(rho = post_mean_dom,saska_mask, eval.im(F_z(x))),
                                             xrange =saska_temp_june$`2013`$xrange, yrange =saska_temp_june$`2013`$yrange)$v)[246:1,])))
img_df2$x <- grid_points$x
img_df2$y <- rev(grid_points$y)


bbox <- data.frame(
  xmin = min(img_df$x),
  xmax = max(img_df$x),
  ymin = min(img_df$y),
  ymax = max(img_df$y)
)


points <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)


plot_saska[[2]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", 
                       colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x)), limits = c(0, 0.006)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = "Longitude", y = "Latitude")+
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
  geom_point(data = points, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 


plot_saska[[3]] <- ggplot(img_df2, aes(x = x, y = y, fill = value)) +
  geom_raster() +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x)), limits = c(0, 0.006)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = "Longitude", y = "Latitude")+
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
  geom_point(data = points, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 


combined_plot <- wrap_plots(plot_saska, nrow = 1)
combined_plot

####################################### when using replicated observations ###################################

rho_list <- Map(function(ppp,cov){as.function(rhohat(ppp, cov, n = 200, method = "ratio"))},
               saska_ppp_june, cov_2_list)
xgrid <- seq(0, 1, length.out = 200)
fmat <- t(sapply(rho_list, function(f) f(xgrid))) * unlist(lapply(saska_ppp_june, function(x){x$n}))
kern_est <- colSums(fmat/sum(unlist(lapply(saska_ppp_june, function(x){x$n}))), na.rm = TRUE)

spline_fit <- smooth.spline(xgrid, kern_est, spar = 0.8)

setwd("~/Research/IPP")
# unlink(paste("~/Research/IPP/", "Adaptive_MCMC_post_multi_saska_exp", sep = ""), recursive = TRUE)
source(file = "pCN_MCMC_GP_mult.R")

discr = seq(0, 1, length.out = 300)

n_iter <- 20000
adapt_multi_MCMC("saska_wind_orig_full", n_iter, beta = 0.1, K = length(discr), 
                saska_ppp_june, grid_points,  
                 cov_3_list, d_k = discr, squared_exponential_kernel, 
                 alpha1 = 2, alpha2 = 2, sigma_2 = 1, nu = 3/2, 
                 exp_param = c(2,2), shape = 1, rate = 2, "lambda")

lambda_post = read.table("Post_l_star.csv", sep = ",")
# lambda_post[lambda_post>10] <- 0
plot(t(lambda_post)[1:dim(lambda_post)[1]], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
plot(t(Ls)[1:dim(Ls)[1]], type = "l")

gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs <- apply(gs, 2, as.numeric)
post_mean_tot <- colMeans(gs[seq(8/10 * dim(gs)[1],10/10*dim(gs)[1]),],na.rm = TRUE)

all_pixels <- unlist(lapply(saska_wind_june, function(img) as.vector(img$v)))
# xx = qnorm(discr)*sd(all_pixels, na.rm = TRUE) + mean(all_pixels, na.rm = TRUE)
# xx = discr * (max(all_pixels, na.rm = TRUE) - min(all_pixels, na.rm = TRUE)) + min(all_pixels, na.rm = TRUE) 
xx = quantile(all_pixels, probs = discr, type = 1, na.rm = TRUE)

rug <- unlist(Map(function(c,l){return(c[l])},saska_wind_june,saska_ppp_june))

plot(xx, kern_est, type = "l")
lines(xx, predict(spline_fit, x = discr)$y, col = "red", lwd = 2, type= "l")
lines(xx, post_mean_tot[1:200], col = "blue", lwd = 2, type= "l")
rug(rug, ticksize=0.03, side=1, lwd=0.5)

upper_tot = apply(gs[seq(dim(gs)[1]/2,dim(gs)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower_tot = apply(gs[seq(dim(gs)[1]/2,dim(gs)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)
polygon(c(xx, rev(xx)), c(upper_tot, rev(lower_tot)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)


plot(im(intensity_est(rho = kern_est,saska_mask, cov_1_list$`2013`),
        xrange = cov_1_list$`2013`$xrange, yrange = cov_1_list$`2013`$yrange),
     main = "Posterior mean of the intensity function with observations",
     zlim = c(0, 0.01))
plot(saska_ppp_june$`2013`, add = TRUE, 
     cex = 0.5, col = "green", pch = "+")
for(i in c(1:19)){
plot(im(intensity_est(rho = post_mean_tot,saska_mask, cov_1_list[[i]]),
        xrange = cov_1_list[[i]]$xrange, yrange = cov_1_list[[i]]$yrange),
     main = paste("Posterior mean for yy", 2003+i))
points(saska_ppp_june[[i]]$x, 
      saska_ppp_june[[i]]$y, 
       col = "green", pch = "+")
}
############################### GGplots for multiple obs ######################################

# ggplots
plot_1d2 <- list()
plot_1d2[[3]] <- ggplot(data.frame(x=xx, y=post_mean_tot[1:300]), aes(x = x, y = y)) + #predict(spline_fit, x = discr)$y
  geom_line(linetype = "solid", color = "gray30", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black"),
    axis.text  = element_text(size = 11)
  ) + labs(x = "Temperature", y = expression(rho(z))) + 
  # geom_line(data = data.frame(x=xx, y=post_mean_tot[1:300]), mapping = aes(x = x, y = y), 
  #           color = "gray20", size = 1) + 
  geom_ribbon(data = data.frame(x = xx, y2 = upper_tot[1:300], y1 = lower_tot[1:300], y=post_mean_tot[1:300]), mapping = aes(x = x, ymin = y1, y = y, ymax = y2),
              fill = "gray20", alpha = 0.2) + 
  geom_rug(data = data.frame(x = rug, y = 0), sides = "b")


img_df2013 <- data.frame(value = as.vector(t((im(intensity_est(rho = post_mean_tot,saska_mask, cov_1_list$`2013`),
                                              xrange = cov_1_list$`2013`$xrange, yrange = cov_1_list$`2013`$yrange)$v)[234:1,])))
img_df2013$x <- grid_points$x
img_df2013$y <- rev(grid_points$y)
points <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)

bbox <- data.frame(
  xmin = min(img_df2013$x),
  xmax = max(img_df2013$x),
  ymin = min(img_df2013$y),
  ymax = max(img_df2013$y)
)

plot_1d[[3]] <- 
  ggplot(img_df2013, aes(x = x, y = y, fill = value)) +
  geom_raster() +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x)), limits = c(0, 0.004)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = "Longitude", y = "Latitude")+
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
  geom_point(data = points, aes(x = x, y = y), shape = 16, color="turquoise1", inherit.aes = FALSE) 


combined_plot <- wrap_plots(plot_1d2, nrow = 1)
combined_plot

######################## Adding precipitation as covariate ##############################

  points <- matrix(c(0,0,0,1,1,1,1,0), ncol = 2, byrow = T)
  conv_hull <- chull(points)
  conv_hull <- rev(c(conv_hull, conv_hull[1]))
  poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
  
  plot(poly_win, main = "Covariate domain", xlab = "Slope", ylab = "Elevetion")
  
  # triangular tesseletion of the convex hull, independent of points 
  
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
  
  plot(mesh, asp = 1, main = "Uniform Triangular Tessellation of Convex Hull")
  points(triangle_centroids, col = "blue", pch = 16, cex = .3)
  vertices <- mesh$P
  extremes <- mesh$P[mesh$PB==1,]
  colnames(extremes) <- c("x", "y")
  colnames(vertices) <- c("x", "y")
  
  inside.owin(x = extremes[,1], y = extremes[,2], w = poly_win)
  
  xrange <- c(0,1)
  yrange <- c(0,1)
  
  # Define a rectangular window
  W_dom <- owin(xrange = xrange, yrange = yrange)
  
  x <- unlist(lapply(cov_1_list, function(img) as.vector(img$v)))
  y <- unlist(lapply(cov_2_list, function(img) as.vector(img$v)))
  points <- cbind(x,y)[complete.cases(cbind(x,y)),]
  # points <- cbind((x-min(elev))/(max(elev)- min(elev)),(y-min(slope))/(max(slope)- min(slope)))
  conv_hull <- chull(points)
  conv_hull <- rev(c(conv_hull, conv_hull[1]))
  poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
  plot(W_dom, main = "Covariate domain", xlab = "Slope", ylab = "Elevetion", asp = 1)
  lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "orange", cex = 2, pch = 16)
  
################################ kernel estimator #######################################
  
  pixels1 <- as.vector(saska_temp_june$`2013`$v)
  F_z1 <- ecdf(pixels1)
  x1 =saska_temp_june$`2013`
  pixels2 <- as.vector(saska_precip_june$`2013`$v)
  F_z2 <- ecdf(pixels2)
  x2 =saska_precip_june$`2013`
  
  rho_kernel <- rho2hat(saska_ppp_june$`2013`, eval.im(F_z1(x1)), eval.im(F_z2(x2)), W = W_dom, dimyx = 400)
  plot(rho_kernel)
  points(eval.im(F_z1(x1))[saska_ppp_june$`2013`], eval.im(F_z2(x2))[saska_ppp_june$`2013`], pch = "+", col = "green")
  
  scale = diff(range(pixels1, na.rm = TRUE))/diff(range(pixels2, na.rm = TRUE))
  
  xgrid_orig <- seq(range(pixels1, na.rm = TRUE)[1], range(pixels1, na.rm = TRUE)[2], length.out = 400)
  ygrid_orig <- seq(range(pixels2, na.rm = TRUE)[1], range(pixels2, na.rm = TRUE)[2], length.out = 400)
  
  
  f_interp <- interp.surface(
    obj = list(x =  rho_kernel$yrow, y = rho_kernel$xcol, z = rho_kernel$v),
    loc = expand.grid(F_z1(xgrid_orig), F_z2(ygrid_orig)) 
  )
  
  f_original_scale <- matrix(f_interp, nrow = 400, byrow = TRUE)
  
  # Plot
  image.plot(xgrid_orig, ygrid_orig, f_original_scale,
        main = "Function on Original Scale", col = temp_pa(100), asp = scale)
  
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
  plot(as.im(intensity2d_est(rho = rho_kernel,saska_mask,
                             Z1 = as.vector(eval.im(F_z1(x1))$v), Z2 = as.vector(eval.im(F_z2(x2))$v)), 
             W =saska_dom),
       main = "Posterior mean of the intensity function with observations")
  plot(saska_ppp_june$`2013`, add = TRUE, pch = "+", cex = .7, col = "green")
  
########################################## single observation ##################################
  
  setwd("~/Research/IPP/")
  source(file = "pCN_adaptive_MCMC_3mesh_GP.R")
  # unlink("~/Research/IPP/2dAdaptive_MCMC_post_saska_2013", recursive = TRUE)
  
  n_iter = 40000
  
  MCMC("_saska_2013_ecdf", n_iter, beta = 0.2, K = dim(vertices)[1], 
      saska_ppp_june$`2013`, grid_points,saska_dom,
       eval.im(F_z1(x1)), eval.im(F_z2(x2)),
       as.data.frame(vertices), triangle_centroids,
       exponential_cov, alpha1=1, alpha2=3, sigma_2 = 1, nu = NA,
       exp_param = c(2,2), shape = 2, rate = 1, 
       link = "lambda", LS = "anisotropic")
  
  lambda_post = read.table("Post_l_star.csv", sep = ",")
  plot(t(lambda_post)[1:dim(lambda_post)[1]], type = "l")
  
  Ls = read.table("Post_ls.csv", sep = ",")
  plot(t(Ls)[2,1:dim(Ls)[1]], type = "l")
  
  theta = read.table("Post_theta.csv", sep = ",", fill = TRUE)
  plot(t(theta)[1,1:dim(theta)[1]], type = "l")
  
  gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs_m <- apply(gs_m, 2, as.numeric)
  pp <- ppp(x = vertices[,1], y = vertices[,2],
            marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(W_dom))
  post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(400,400), method = "pwl")
  plot(post_mean_m, main = "Posterior for single year data")
  
  points(vertices, col = "black", cex = 1, pch = 16)
  points(cov_1_list$`2013`[saska_ppp_june$`2013`], cov_2_list$`2013`[saska_ppp_june$`2013`], pch = "+", col = "white", cex = .8)
  
  pp_uncert <- ppp(x = vertices[,1], y = vertices[,2],
                   marks = apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.95) - 
                     apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.05), window = owin(poly_win))
  
  post_uncert <- Smooth(pp_uncert, kernel = "gaussian", dimyx = c(dim(mesh$P)[1], dim(mesh$P)[1]), method = "pwl")
  
  plot(post_uncert, 
       main = "Posterior bandwidth on unit square")
  
  plot(as.im(intensity2d_est(rho = post_mean_m,saska_mask,
                             Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
             W =saska_dom),
       main = "Posterior mean of the intensity function with observations")
  
  plot(saska_ppp_june$`2013`, add = TRUE, pch = "+", cex = .7, col = "white")
  
  
  
####################################### ggplot single year ######################################
  ggcube <- list()
  
  img_df <- data.frame(value = as.vector(rho_kernel$v))
  img_df$x <- rep(rho_kernel$yrow, 400)
  img_df$y <- rep(rho_kernel$xcol, each = 400) 
  
  points <- data.frame(x = cov_1_list$`2013`[saska_ppp_june$`2013`],
                       y = cov_2_list$`2013`[saska_ppp_june$`2013`])
  img_df$value2 <- as.vector(post_mean_m$v)
  
  ggcube[[1]] <- ggplot() +
    geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.007)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Temperature", y = "Precipitation") + 
    geom_point(data = points, aes(x = x, y = y), shape = 3, color="white")
  
  ggcube[[2]] <- ggplot() +
    geom_raster(data = img_df, aes(x = x, y = y, fill = value2)) +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.007)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Temperature", y = "Precipitation") + 
    geom_point(data = points, aes(x = x, y = y), shape = 3, color="white")
  
  est2013 <- data.frame(kern = as.vector(t((as.im(intensity2d_est(rho = rho_kernel,saska_mask,
                                                  Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
                                                  W =saska_dom)$v)[246:1,])))
  est2013$x <- grid_points$x
  est2013$y <- rev(grid_points$y)
  points <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  
  ggcube[[3]] <- ggplot(est2013, aes(x = x, y = y, fill = kern)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.006)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  est2013$gp <- as.vector(t((as.im(intensity2d_est(rho = post_mean_m,saska_mask,
                                   Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
                                   W =saska_dom)$v)[246:1,]))
  
  
  ggcube[[4]] <- ggplot(est2013, aes(x = x, y = y, fill = gp)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.006)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  combined_plot <- wrap_plots(ggcube, nrow = 2)
  combined_plot
  
################## when adding data from other years, like our method can naturally account for ###############
  
  kern_years <- lapply(c(1:length(saska_ppp_june)), function(i){rho2hat(saska_ppp_june[[i]], 
                       cov_1_list[[i]], cov_2_list[[i]], method = "ratio", dimyx = c(800,800))$v})
  kern2_est = rho_kernel
  kern2_est$v <- apply(simplify2array(kern_years), c(1, 2), mean, na.rm = TRUE)
  plot(kern2_est)
  
  plot(as.im(intensity2d_est(rho = kern2_est,saska_mask,
                             Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
             W =saska_dom),
       main = "Posterior mean of the intensity function with observations")
  
  plot(saska_ppp_june$`2013`, add = TRUE, pch = "+", cex = .7, col = "white")
  
  rho2_list <- Map(function(ppp,cov1, cov2){as.function(rho2hat(ppp, cov1, cov2,  dimyx = c(400,400), method = "ratio"))},
                  saska_ppp_june, cov_1_list, cov_2_list)
  
  fmat <- t(sapply(rho2_list, function(f) f(vertices[,1], vertices[,2]))) * unlist(lapply(saska_ppp_june, function(x){x$n}))
  pp_kern <- ppp(x = vertices[,1], y = vertices[,2],
                marks = colSums(fmat/sum(unlist(lapply(saska_ppp_june, function(x){x$n}))),na.rm = TRUE), window = owin(W_dom))
  kern_tot <- Smooth(pp_kern, kernel = "gaussian", dimyx = c(800,800), method = "pwl", pad = TRUE, edge = FALSE)
  kern_tot$xcol <- seq(kern_tot$yrange[1], kern_tot$yrange[2], length.out = 800)
  kern_tot$yrow <- seq(kern_tot$xrange[1], kern_tot$xrange[2], length.out = 800)
  plot(kern_tot, main = "Kernel intensity estimator over aggregated data")
  
  lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "green", lwd = 3)
  
  pixels1 <- unlist(lapply(saska_temp_june, function(x){as.vector(x$v)}))
  F_z1 <- ecdf(pixels1)
  pixels2 <-  unlist(lapply(saska_precip_june, function(x){as.vector(x$v)}))
  F_z2 <- ecdf(pixels2)
  
  scale = diff(range(pixels1, na.rm = TRUE))/diff(range(pixels2, na.rm = TRUE))
  
  xgrid_orig <- seq(range(pixels1, na.rm = TRUE)[1]+1e-06, range(pixels1, na.rm = TRUE)[2]-1e-06, length.out = 400)
  ygrid_orig <- seq(range(pixels2, na.rm = TRUE)[1]+1e-06, range(pixels2, na.rm = TRUE)[2]-1e-06, length.out = 400)
  
  
  f_interp <- interp.surface(
    obj = list(x =  kern_tot$yrow, y = kern_tot$xcol, z = kern_tot$v),
    loc = expand.grid(F_z1(xgrid_orig), F_z2(ygrid_orig)) 
  )
  
  f_original_scale <- matrix(f_interp, nrow = 400, byrow = TRUE)
  
  # Plot
  image.plot(xgrid_orig, ygrid_orig, f_original_scale,
        main = "Function on Original Scale", col = my_col_map(256), asp = scale)
  lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "green", lwd = 3)
  
##################################### multiple observations ###########################################
  setwd("~/Research/IPP/")
  source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
  # unlink(paste("~/Research/IPP/", "2dAdaptive_MCMC_post_multi_canadian_fires_multiple_years", sep = ""), recursive = TRUE)
  
  n_iter  = 30000
  
  MCMC("2dAdaptive_MCMC_multi_saska_ecdf", n_iter, beta = 0.1, 
       K = dim(vertices)[1],saska_ppp_june, grid_points,saska_dom,
       cov_1_list, cov_2_list,
       as.data.frame(vertices), triangle_centroids,
       exponential_cov, alpha1=1, alpha2=2, sigma_2 = 1, nu = NA,
       exp_param = c(2,2), shape = 2, rate = 1, 
       link = "lambda", LS = "anisotropic",
       Index_over_dom_hole_19, Index_over_loc_prec_hole_19)
  
  
  lambda_post = read.table("Post_l_star.csv", sep = ",")
  plot(t(lambda_post)[1:dim(lambda_post)[1]], type = "l")
  
  Ls = read.table("Post_ls.csv", sep = ",")
  plot(t(Ls)[2,1:dim(Ls)[1]], type = "l")
  
  theta = read.table("Post_theta.csv", sep = ",", fill = TRUE)
  plot(t(theta)[2,1:dim(theta)[1]], type = "l")
  
  
  gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs_m <- apply(gs_m, 2, as.numeric)
  pp_tot <- ppp(x = vertices[,1], y = vertices[,2],
            marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(W_dom))
  post_mean_tot <- Smooth(pp_tot, kernel = "quartic",  dimyx = c(800,800), method = "pwl", pad = TRUE, edge = FALSE)
  post_mean_tot$xcol <- seq(post_mean_tot$yrange[1], post_mean_tot$yrange[2], length.out = 800)
  post_mean_tot$yrow <- seq(post_mean_tot$xrange[1], post_mean_tot$xrange[2], length.out = 800)
  
  plot(post_mean_tot, main = "Posterior mean intensity over aggregated data")
  
  points(vertices, col = "black", cex = 1, pch = 16)
  
  Map(function(c1,c2,l) points(c1[l], c2[l], ex = 0.7, col = "white", pch = "+"),
      cov_1_list, cov_2_list,saska_ppp_june)
  
  
  post_orig <- interp.surface(
    obj = list(x =  post_mean_tot$yrow, y = post_mean_tot$xcol, z = post_mean_tot$v),
    loc = expand.grid(F_z1(xgrid_orig), F_z2(ygrid_orig)) 
  )
  
  post_original_scale <- matrix(post_orig, nrow = 400, byrow = TRUE)
  
  # Plot
  image.plot(xgrid_orig, ygrid_orig, post_original_scale,
             main = "Function on Original Scale", col = my_col_map(256), asp = scale)
  lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "green", lwd = 3)
  
  
  pp_uncert <- ppp(x = vertices[,1], y = vertices[,2],
                   marks = apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.95) - 
                     apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, quantile, probs = 0.05), window = owin(W_dom))
  
  post_uncert <- Smooth(pp_uncert, kernel = "quartic",  dimyx = c(800,800), method = "pwl", pad = TRUE, edge = FALSE)
  post_uncert$xcol <- seq(post_uncert$yrange[1], post_uncert$yrange[2], length.out = 800)
  post_uncert$yrow <- seq(post_uncert$xrange[1], post_uncert$xrange[2], length.out = 800)
  
  
  post_orig_uncert <- interp.surface(
    obj = list(x =  post_uncert$yrow, y = post_uncert$xcol, z = post_uncert$v),
    loc = expand.grid(F_z1(xgrid_orig), F_z2(ygrid_orig)) 
  )

  uncert_original_scale <- matrix(post_orig_uncert, nrow = 400, byrow = TRUE)
  # Plot
  image.plot(xgrid_orig, ygrid_orig, uncert_original_scale,
             main = "Function on Original Scale", col = my_col_map(256), asp = scale)
  lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "green", lwd = 3)
  
  
  pp_var <- ppp(x = vertices[,1], y = vertices[,2],
                   marks = apply(gs_m[seq(8/12 * (dim(gs_m)[1]),(dim(gs_m)[1])),], 2, function(x){var(x, na.rm = TRUE)/1}), window = owin(W_dom))
  
  post_var <- Smooth(pp_var, kernel = "quartic",  dimyx = c(800,800), method = "pwl", pad = TRUE, edge = FALSE)
  post_var$xcol <- seq(post_var$yrange[1], post_var$yrange[2], length.out = 800)
  post_var$yrow <- seq(post_var$xrange[1], post_var$xrange[2], length.out = 800)
  
  
  post_orig_var <- interp.surface(
    obj = list(x =  post_var$yrow, y = post_var$xcol, z = post_var$v),
    loc = expand.grid(F_z1(xgrid_orig), F_z2(ygrid_orig)) 
  )
  
  var_original_scale <- matrix(post_orig_var, nrow = 400, byrow = TRUE)
  
  # Plot
  image.plot(xgrid_orig, ygrid_orig, var_original_scale,
             main = "Function on Original Scale", col = my_col_map(256), asp = scale)
  lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "green", lwd = 3)
  
  plot(as.im(intensity2d_est(rho = post_uncert,saska_mask,
                             Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
             W =saska_dom),
       main = "Posterior mean of the intensity function with observations")
  
  plot(saska_ppp_june$`2013`, add = TRUE, pch = "+", cex = .7, col = "white")
  
  for(i in 1:21){
    plot(as.im(intensity2d_est(rho = post_mean_tot,saska_mask, 
                               Z1 = as.vector(cov_1_list[[i]]), Z2 = as.vector(cov_2_list[[i]])), 
               W =saska_dom),
         main = paste("Posterior mean for yy", 2003+i))
    points(saska_ppp_june[[i]]$x, 
          saska_ppp_june[[i]]$y, 
           col = "green", pch = "+")
  }
  
####################################### ggplot multiple years ######################################
  ggcube <- list()
  
  all_pixelsX <- unlist(lapply(saska_temp_june, function(img) as.vector(img$v)))
  all_pixelsY <- unlist(lapply(saska_precip_june, function(img) as.vector(img$v)))
  
  img_df <- data.frame(value = (as.vector(post_mean_tot$v))) # data.frame(value = (as.vector(t(post_original_scale))))
  
  exp_cov <- expand.grid(x = seq(0, 1, length.out = 800), 
                         y = seq(0, 1, length.out = 800))
  
  img_df$x <- exp_cov$y # rep(xgrid_orig, each = 400)
  img_df$y <- exp_cov$x #rep(ygrid_orig) 
  
  points <- data.frame(x = unlist(Map(function(c, l){c[l]}, cov_1_list,saska_ppp_june)),
                       y = unlist(Map(function(c, l){c[l]}, cov_2_list,saska_ppp_june)))
  img_df$value2 <- as.vector(kern2_est$v) # as.vector(t(f_original_scale))
  img_df$value3 <- as.vector(t(var_original_scale))
  
  
  df_bdry <- data.frame(
    x = poly_win$bdry[[1]]$x,
    y = poly_win$bdry[[1]]$y
  )

  ggcube[[2]] <- ggplot() +
    geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.0071)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Temperature", y = "Precipitation") + 
    # geom_path(data = df_bdry, aes(x = x, y = y), color = "turquoise2", linewidth = 1.5) +
    #geom_point(data = points, aes(x = x, y = y), shape = 19, size = 0.5, color="turquoise1")+
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 
  
  ggcube[[3]] <- ggplot() +
    geom_raster(data = img_df, aes(x = x, y = y, fill = value2)) +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "viridis", name = expression(rho(z))) +
    # scale_fill_gradientn(colors = colors, name = expression(Var(z))) +
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
                       probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) #+ 
    geom_path(data = df_bdry, aes(x = x, y = y), color = "turquoise2", linewidth = 1.5) +
    geom_point(data = points, aes(x = x, y = y), shape = 19, size = 0.5, color="turquoise1")
  
  ggcube[[1]] <- ggplot() +
    # geom_raster(data = img_df, aes(x = x, y = y, fill = value2)) +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    # scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
    #                      name = expression(rho(z))) +
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
  
  est <- data.frame(mean13 = as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                                  Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
                                                  W =saska_dom)$v)[234:1,])))
  est$x <- grid_points$x
  est$y <- rev(grid_points$y)
  points13 <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  
  ggcube[[1]] <- ggplot(est, aes(x = x, y = y, fill = mean13)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x)), limits = c(0, 0.0065)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points13, aes(x = x, y = y), shape = 16, size = 1, color="turquoise1", inherit.aes = FALSE) 
  
  est$mean15 <- as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                   Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
                                   W =saska_dom)$v)[234:1,]))
  
  
  points15 <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  
  ggcube[[2]] <- ggplot(est, aes(x = x, y = y, fill = mean15)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x))) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points15, aes(x = x, y = y),  shape = 16, size = 1, color="turquoise1", inherit.aes = FALSE) 
  
  est$mean21 <- as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                   Z1 = as.vector(cov_1_list$`2013`$v), Z2 = as.vector(cov_2_list$`2013`$v)), 
                                   W =saska_dom)$v)[234:1,]))
  
  
  points21 <- data.frame(x =saska_ppp_june$`2013`$x, y =saska_ppp_june$`2013`$y)
  
  ggcube[[3]] <- ggplot(est, aes(x = x, y = y, fill = mean21)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(lambda(x))) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points21, aes(x = x, y = y),  shape = 16, size = 1, color="turquoise1", inherit.aes = FALSE) 
  
  
  combined_plot <- wrap_plots(ggcube, nrow = 1)
  combined_plot
  
  #### addtional ggplots for 2 covariates ###
  
  img_df <- data.frame(value = as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                                  Z1 = as.vector(cov_1_list$`2006`$v), Z2 = as.vector(cov_2_list$`2006`$v)), 
                                                  W =saska_dom)$v)[246:1,])))
  img_df$x <- grid_points$x
  img_df$y <- rev(grid_points$y)
  
  img_df$value2 <- as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                      Z1 = as.vector(cov_1_list$`2011`$v), Z2 = as.vector(cov_2_list$`2011`$v)), 
                                      W =saska_dom)$v)[246:1,]))
  
  img_df$value3 <- as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                      Z1 = as.vector(cov_1_list$`2017`$v), Z2 = as.vector(cov_2_list$`2017`$v)), 
                                      W =saska_dom)$v)[246:1,]))
  
  img_df$value4 <- as.vector(t((as.im(intensity2d_est(rho = post_mean_tot,saska_mask,
                                                      Z1 = as.vector(cov_1_list$`2019`$v), Z2 = as.vector(cov_2_list$`2019`$v)), 
                                      W =saska_dom)$v)[246:1,]))
  
  
  points1 <- data.frame(x =saska_ppp_june$`2006`$x, y =saska_ppp_june$`2006`$y)
  points2 <- data.frame(x =saska_ppp_june$`2011`$x, y =saska_ppp_june$`2011`$y)
  points3 <- data.frame(x =saska_ppp_june$`2017`$x, y =saska_ppp_june$`2017`$y)
  points4 <- data.frame(x =saska_ppp_june$`2019`$x, y =saska_ppp_june$`2019`$y)
  
  plot3d <- list()
  
  plot3d[[1]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.007)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points1, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  plot3d[[2]] <- ggplot(img_df, aes(x = x, y = y, fill = value2)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.007)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points2, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  
  plot3d[[3]] <- ggplot(img_df, aes(x = x, y = y, fill = value3)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.007)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points3, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  plot3d[[4]] <- ggplot(img_df, aes(x = x, y = y, fill = value4)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(na.value = "white", 
                         colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 0.007)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = "Longitude", y = "Latitude")+
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
    geom_point(data = points4, aes(x = x, y = y), shape = 3, color="white", inherit.aes = FALSE) 
  
  
  combined_plot <- wrap_plots(plot3d, nrow = 2)
  combined_plot
  
