CODE and DATA supporting "A nonparametric Bayesian analysis of independent and identically distributed observations of covariate-driven Poisson processes" by Patric Dolmeta and Matteo Giordano: https://arxiv.org/abs/2509.02299

Files named CODE_ are R scripts to repeat the experiments/analysis in the paper. 
Follow comments for 1to1 referencing with Sections/Tables and Figures in the paper. 

RData files are the datasets as used in the Canadian wildfires application. 
  
In CODE_DataGeneration_Analysis.R:  univariate covariate process generation 
                                    point pattern simulation
                                    1D frequentist and Bayesian analysis
                                    bivariate covariate process generation 
                                    triangukar mesh construction on domain
                                    2D frequentist and Bayesian analysis

In CODE_Canada_processing.R:        download and processing of Canadian weather and wildfire location maps from https://cwfis.cfs.nrcan.gc.ca/maps/wx
                                    single region data extraction 

CODE_Canada_analysis.R:             preparation of canadian covariate maps for analysis
                                    1D preliminary analysis with temperature/precipitation or wind speed, separately
                                    2D joint analysis
                                    3D analysis 

Adaptive_Anisotropic_GP.R:          Markov Chain Monte Carlo sampler for Adaptive Anisotropic Gaussian Process prior
                                    independent on number of covariates (input as list) and corresponding discretization (input as matrix)
                                    working on parallel architecture (multisession on Windows, multicore otherwise)
