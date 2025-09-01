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

################################# SEMPRE: CARICAMENTO ESPERIMENTI ##############################################
{
  setwd("~/Research/IPP/2dAdaptive_MCMC_post_mesh_replicated_obs")
  load("~/Research/IPP/2dAdaptive_MCMC_post_mesh_replicated_obs/ed_exp_kernel.RData")
}

########################################### LANCIARE SIMULAZIONI ###############################################
{
  # mettere estremi forl-loop "M = 2,10,50,250,500,1000,2000" ed "exp = 1:50/1:100"
  # set parametro "beta" tra 0 e 1. maggiore M -> minore beta
  i <- 1
  bet <- c(0.1)
  for(M in c(1000)){
    setwd("~/Research/IPP/")
    source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
    folder = paste("2dAdaptive_MCMC_post_mesh_slope_repl/", M, sep = "")
    dir.create(paste("~/Research/IPP/", folder, sep = ""))
    setwd(folder)
    for(exp in c(1:5)){
      setwd("~/Research/IPP/")
      source(file = "Adaptive_pCN_MCMC_2dGP_mult_mesh.R")
      n_iter  = 20000
      MCMC(paste("2dAdaptive_MCMC_post_mesh_slope_repl/",M, "/", exp,sep = ""), n_iter,
           beta = 0.08, K = dim(vertices)[1], 
           loc2d_list_slope[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
           grid_points, window,
           cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
           cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
           as.data.frame(vertices), triangle_centroids,
           exponential_cov, alpha1=1, alpha2=2, sigma_2 = 1, nu = NA,
           exp_param = c(2,2), shape = 2, rate = 0.5, 
           link = "lambda", LS = "anisotropic")
    }
  i <- i+1
  }
}

########################################### PRODURRE TABELLE ###################################################
{
  # mettere estremi forl-loop "M = 2,10,50,250,500,1000,2000" ed "exp = 1:50/1:100"
  L2_2dbay_exp = matrix(NA,5,10)
  L2_2dbay_rel_exp = matrix(NA,5,10)
  i = 1
  exp_cov <- expand.grid(x = seq(1.084202e-19, 1, length.out = 800), 
                         y = seq(1.084202e-19, 1, length.out = 800))
  
  rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                        sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
                 xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
  for(M in c(1000)){
    setwd("~/Research/IPP/")
    folder = paste("2dAdaptive_MCMC_post_mesh_hole_repl/", M, sep = "")
    setwd(folder)
    for(exp in c(1:4)){
      folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_hole_repl/", M, "/", exp, sep = "")
      setwd(folder)
      gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
      gs_m <- apply(gs_m, 2, as.numeric)
      pp <- ppp(x = vertices[,1], y = vertices[,2],
                marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
      post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(800, 800), method = "pwl")
      print(paste(M,exp))
      L2_2dbay_exp[i,exp]<- sqrt(sum((post_mean_m - rho_imag)^2)/length(exp_cov$x))
      L2_2dbay_rel_exp[i,exp]<-  sqrt(sum((post_mean_m - rho_imag)^2))/sqrt(sum((rho_imag)^2))
    }
    i<-i+1
  }
  
  funcs <- list(mean = mean, sd = sd)
  
  result_rel <- lapply(funcs, function(f) apply(L2_2dbay_rel_exp, 1, f, na.rm = TRUE))
  print(paste(round(result_rel$mean,2), "(", round(result_rel$sd,2), ")&"))
  
  result <- lapply(funcs, function(f) apply(L2_2dbay_exp, 1, f, na.rm = TRUE))
  print(paste(round(result$mean,2), "(", round(result$sd,2), ")&"))
  
  save.image("~/Research/IPP/Adaptive_MCMC_post_multi/salvataggio_albania.RData")
}

########################################### PRODURRE IMMAGINI ###############################################
{
  plot_list2 <- list()
  i = 1
  for(M in c(50, 250, 1000)){
    
    setwd("~/Research/IPP/")
    folder = paste("2dAdaptive_MCMC_post_mesh_replicated_obs/repexp/", M, sep = "")
    setwd(folder)
    
    exp = 3
    
    folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_replicated_obs/repexp/", M, "/", exp, sep = "")
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
    
    plot_list2[[1]] <-
      ggplot() +
      geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  
      scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"), 
                           name = expression(rho(z)), limits = c(0, 55)) +
      coord_fixed() +  # Keeps aspect ratio correct
      theme_minimal() + 
      labs(x = expression(z[1]), y = expression(z[2])) + theme(legend.position = "none")
    
    i<-i+1
  }
  
  # add ground truth 
  
  exp_cov <- expand.grid(x = seq(0,1, length.out = 200), 
                         y = seq(0,1, length.out = 200))
  
  rho_imag <- im(matrix(rho2d(exp_cov$x, exp_cov$y),
                        sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
                 xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
  img_df <- data.frame(value = as.vector(rho_imag$v))
  img_df$x <- gridcentres(window, nx = 200, ny = 200)$y
  img_df$y <- gridcentres(window, nx = 200, ny = 200)$x
  
  plot_list2[[4]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +  # or geom_tile() if you prefer
    # scale_fill_viridis(option = "magma", name = expression(rho(z)), limits = c(0,40)) +
    # scale_fill_gradientn(colors = colors, name = expression(rho(z)), limits = c(0, 120)) +
    scale_fill_gradientn(colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                         name = expression(rho(z)), limits = c(0, 55)) +
    coord_fixed() +  # Keeps aspect ratio correct
    theme_minimal() + 
    labs(x = expression(z[1]), y = expression(z[2]))
  
  
  combined_plot <- wrap_plots(plot_list2, nrow = 1)
  combined_plot
}  

############################################ convergence diagnostics #######################################
col = brewer.pal(8, "Dark2")
n_iter <- 20000
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


exp_cov <- expand.grid(x = seq(0,1, length.out = 200), 
                       y = seq(0,1, length.out = 200))

which.min(rowSums((sweep(vertices, 2, c(0.8,0.3), "-"))^2))
which.min(rowSums((sweep(vertices, 2, c(0.3,0.8), "-"))^2))
which.min(rowSums((sweep(vertices, 2, c(0.5,1), "-"))^2))

RHO2 <- rho2d(vertices[,1], vertices[,2])

for(exp in c(1:5)){
  setwd("~/Research/IPP/2dAdaptive_MCMC_post_mesh_hole_repl/1000")
  folder = paste("~/Research/IPP/2dAdaptive_MCMC_post_mesh_hole_repl/1000", "/", exp, sep = "")
  setwd(folder)
  
  # lambda_post = read.table("Post_l_star.csv", sep = ",")
  # 
  # theta = read.table("Post_theta.csv", sep = ",")
  # 
  # Ls = read.table("Post_ls.csv", sep = ",")
  
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs <- apply(gs, 2, as.numeric)
  # 
  # Z1 <-Map(function(cov, lo){cov[lo]},
  #         cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
  #         loc2d_list_bumphole[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  # Z2 <-Map(function(cov, lo){cov[lo]},
  #          cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
  #          loc2d_list_bumphole[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  # 
  # list_index_dom <- Map(function(cov1, cov2){unlist(apply(cbind(na.omit(as.vector(cov1$v)),na.omit(as.vector(cov2$v))),
  #                                                         1, function(w){which.min(rowSums((vertices - matrix(w, nrow = dim(vertices)[1], ncol = 2, byrow = TRUE))^2))}))},
  #                       cov1_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], cov2_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  # list_index_loc <- Map(function(z1, z2){apply(cbind(z1,z2), 1,
  #                       function(w){which.min(rowSums((vertices - matrix(w, nrow = dim(vertices)[1], ncol = 2, byrow = TRUE))^2))})},
  #                       Z1, Z2)

  # log_lik <- rbind(log_lik,apply(as.data.frame(gs), 1, function(r){
  #   rho_obs <- log(r[unlist(list_index_loc)])
  #   I <- lapply(list_index_dom, function(i){sum(r[i]-1, na.rm = TRUE)/2500})
  #   (sum(rho_obs) - sum(unlist(I)))/1000
  # }))

  # true_lik[exp] <- (sum(log(RHO2[unlist(list_index_loc)])) -
  #   sum(unlist(lapply(list_index_dom, function(i){sum(RHO2[i]-1, na.rm = TRUE)/2500}))))/1000

  # trace_plot1 <- trace_plot1 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(lambda_post)[1:n_iter]), mapping = aes(x = x, y = y),
  #                     linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho^"*"))
  # trace_plot2 <- trace_plot2 + geom_line(data = data.frame(x=seq(1:n_iter), y=Ls[1:n_iter,1]), mapping = aes(x = x, y = y),
  #                     linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(ℓ[1]))
  # trace_plot4 <- trace_plot4 + geom_line(data = data.frame(x=seq(1:n_iter), y=Ls[1:n_iter,2]), mapping = aes(x = x, y = y),
                                         # linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(ℓ[2]))
  # trace_plot3 <- trace_plot3 + geom_line(data = data.frame(x=seq(2:n_iter), y=log_lik[exp+6,2:n_iter]), mapping = aes(x = x, y = y),
  #                                        linetype = "solid", color = col[exp], size = 1) + labs(x = "Number of iterations", y = "log-likelihood") +
  #                              geom_hline(yintercept = 0.995*true_lik[exp], color = col[exp], linetype = "dashed", size = 1)
  # trace_plot5 <- trace_plot5 + geom_line(data = data.frame(x=seq(1:n_iter), y=theta[1:n_iter,1]), mapping = aes(x = x, y = y),
  #                                        linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(theta[1]))
  # trace_plot6 <- trace_plot6 + geom_line(data = data.frame(x=seq(1:n_iter), y=theta[1:n_iter,2]), mapping = aes(x = x, y = y),
  #                                        linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(theta[2]))

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
