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
library(latex2exp)

################################# SEMPRE: CARICAMENTO ESPERIMENTI ##############################################
{
  setwd("~/Research/IPP/Adaptive_MCMC_post_multi")
  load("~/Research/IPP/Adaptive_MCMC_post_multi/plots_1d_rep_exp.RData")
}

########################################### LANCIARE SIMULAZIONI ###############################################
{
  # mettere estremi forl-loop "M = 2,10,50,250,500,1000,2000" ed "exp = 1:50/1:100"
  # set parametro "beta" tra 0 e 1. maggiore M -> minore beta
  beta = c(0.8,0.5,0.15,0.1,0.08,0.06) 
  i <- 5
  for(M in c(500,1000)){
    setwd("~/Research/IPP/")
    source(file = "pCN_MCMC_GP_mult.R")
    folder = paste("Adaptive_MCMC_post_multi", "/min_max/", M, sep = "")
    dir.create(paste("~/Research/IPP/", folder, sep = ""))
    setwd(folder)
    for(exp in c(1:5)){
      setwd("~/Research/IPP/")
      source(file = "pCN_MCMC_GP_mult.R")
      discr <- seq(0,1, length.out = 200)
      n_iter <- 20000
      adapt_multi_MCMC(paste0("min_max/", M, "/", exp, sep = ""), n_iter, beta = beta[i], K = length(discr), 
                       loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], grid_points, 
                       cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], d_k = discr,
                       squared_exponential_kernel, alpha1 = 1, alpha2 = 1, 
                       sigma_2 = 1, nu = 3/2, exp_param = c(2,2),
                       shape = 3, rate = 1/3, "lambda")
    }
    i <- i+1
  }   
}

########################################### PRODURRE TABELLE ###################################################
{
  # mettere estremi forl-loop "M = 2,10,50,250,500,1000,2000" ed "exp = 1:50/1:100"
  L2_bay_exp = matrix(NA,5,100)
  L2_bay_rel_exp = matrix(NA,5,100)
  i = 1
  for(M in c(500)){
    setwd("~/Research/IPP/")
    folder = paste("Adaptive_MCMC_post_multi/neg_exp/", M, sep = "")
    setwd(folder)
    for(exp in c(1:5)){
      folder = paste("~/Research/IPP/Adaptive_MCMC_post_multi/neg_exp/", M, "/", exp, sep = "")
      setwd(folder)
      gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
      gs <- apply(gs, 2, as.numeric)
      post_mean <- colMeans(gs[seq(n_iter/4,n_iter),],na.rm = TRUE)
      print(exp)
    
      L2_bay_exp[i,exp]<- sqrt(sum((post_mean - rho(discr))^2)/length(discr))
      L2_bay_rel_exp[i,exp]<- sqrt(sum((post_mean - rho(discr))^2)/(sum(rho(discr)^2)))
    }
    i<-i+1
  }
  
  funcs <- list(mean = mean, sd = sd)
  
  result_rel <- lapply(funcs, function(f) apply(L2_bay_rel_exp, 1, f, na.rm = TRUE))
  print(paste(round(result_rel$mean,2), "(", round(result_rel$sd,3), ")&"))
  
  result <- lapply(funcs, function(f) apply(L2_bay_exp, 1, f, na.rm = TRUE))
  print(paste(round(result$mean,2), "(", round(result$sd,3), ")&"))

  save.image("~/Research/IPP/Adaptive_MCMC_post_multi/salvataggio_albania.RData")
}

########################################### PRODURRE IMMAGINI ###############################################
{
  i <- 1
  for(M in c(250, 500, 1000)){
  setwd("~/Research/IPP/")
  folder = paste("Adaptive_MCMC_post_multi/neg_exp/", M, sep = "")
  setwd(folder)
  for(exp in c(1:5)){
  # exp <- c(1,1,4)[i]
  folder = paste("~/Research/IPP/Adaptive_MCMC_post_multi/neg_exp/", M, "/", exp, sep = "")
  setwd(folder)
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs <- apply(gs, 2, as.numeric)
  post_mean <- colMeans(gs[seq(5/10*dim(gs)[1],dim(gs)[1]),],na.rm = TRUE)
    
  upper = apply(gs[seq(3/10*dim(gs)[1],dim(gs)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
  lower = apply(gs[seq(3/10*dim(gs)[1],dim(gs)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)
  
  rug <- unlist(Map(function(c,l){return(c[l])} , cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)], 
                    loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)]))
  
  plot_list[[i]] <-
    plot_list[[i]] +
  # print(ggplot(data.frame(x=discr, y=rho(discr)), aes(x = x, y = y)) +
    # geom_line(linetype = "dashed", color = "gray30", size = 1) +
      geom_line(data = data.frame(x=discr, y=post_mean), mapping = aes(x = x, y = y),
                color = "dodgerblue3", size = 1) +
      geom_ribbon(data = data.frame(x = discr, y2 = upper, y1 = lower, y=post_mean), mapping = aes(x = x, ymin = y1, y = post_mean, ymax = y2),
                  fill = "dodgerblue3", alpha = 0.2) #+ ylim(0, 15) #+
      geom_rug(data = data.frame(x = rug, y = 0), sides = "b")
    
  i<-i+1
  }
  }

combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot
save.image("~/Research/IPP/Adaptive_MCMC_post_multi/salvataggio_albania.RData")
}


############################################ convergence diagnostics #######################################
col = brewer.pal(12,"Paired")
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

for(exp in c(2,4,6,8,10,12)){
  setwd("~/Research/IPP/Adaptive_MCMC_post_multi/1000")
  folder = paste("~/Research/IPP/Adaptive_MCMC_post_multi/1000", "/", exp, sep = "")
  setwd(folder)
  
  # lambda_post = read.table("Post_l_star.csv", sep = ",")
  # plot(t(lambda_post)[1:n_iter], type = "l")
  
  # theta = read.table("Post_theta.csv", sep = ",")
  
  # Ls = read.table("Post_ls.csv", sep = ",")
  # plot(t(Ls)[1:n_iter], type = "l")
  
  gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
  gs <- apply(gs, 2, as.numeric)
  
  # Z <-Map(function(cov, lo){cov[lo]},
  #         cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
  #         loc_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)])
  # index_dom <- lapply(cov_list[((exp-1)*1000 + 1):((exp-1)*1000 + M)],
  #                     function(c){pmin(round(c$v*200), 200)})
  # index_loc <- lapply(Z, function(c){pmin(round(c*200), 200)})
  # log_lik <- apply(as.data.frame(gs[seq(1,n_iter,1),]), 1, function(rho){
  #   rho_obs <- log(rho[unlist(index_loc)])
  #   I <- lapply(index_dom, function(i){sum(rho[i]-1, na.rm = TRUE)/2500})
  #   sum(rho_obs) - sum(unlist(I))
  # })
  # true_lik <- sum(log(rho(discr)[unlist(index_loc)])) -
  #   sum(unlist(lapply(index_dom, function(i){sum(rho(discr)[i]-1, na.rm = TRUE)/2500})))

  # trace_plot1 <- trace_plot1 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(lambda_post)[1:n_iter]), mapping = aes(x = x, y = y, color = col[exp]),
  # linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(rho^"*"))
  # trace_plot2 <- trace_plot2 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(Ls)[1:n_iter]), mapping = aes(x = x, y = y),
  #                                        linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(l))
  # trace_plot3 <- trace_plot3 + geom_line(data = data.frame(x=seq(2:n_iter), y=t(log_lik)[2:n_iter]), mapping = aes(x = x, y = y),
  #                                        linetype = "solid", color = col[exp], size = 0.3) + labs(x = "Number of iterations", y = "log-likelihood") +
  #                              geom_hline(yintercept = true_lik, color = col[exp], linetype = "dashed", size = 0.5)
  # 
  # trace_plot4 <- trace_plot4 + geom_line(data = data.frame(x=seq(1:n_iter), y=t(theta)[1:n_iter]), mapping = aes(x = x, y = y),
  #                                        linetype = "solid", color = col[exp], size = 0.3, alpha = 0.6) + labs(x = "Number of iterations", y = expression(theta))
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
