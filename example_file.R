rm(list = ls())

library(scalreg)
library(glmnet)
library(doParallel)
library(foreach)
library(MASS)
require(corpcor)
library(dplyr)



ncores <- 11
cl <- makeCluster(ncores)   
registerDoParallel(cl)




#source codes
source("all_source.R")


#file path
#my_path <- "/scratch/sosarfo/logerror_new_codes/" 





sample_size<-200
taxa_size<-300
conf_size<-4
sig_u <- 1
sig_factor <- 0.6
nz_alp<-7
alpha_truth<-c(1,-0.8,1.5,0.6,-0.9,1.2, 0.4, rep(0, taxa_size-nz_alp-1) )
beta_truth<- rnorm(conf_size)
nsim <- 100
seednum<-2025


#file_name = sprintf("%ssim_res_%s_%s_%s_%s.RData", my_path, sample_size, taxa_size, conf_size, "checks")



results <- sim_run(
  sample_size=sample_size,
  taxa_size =taxa_size, 
  conf_size = conf_size, 
  alpha_truth =alpha_truth, 
  beta_truth  = beta_truth, 
  nsim = nsim,
  seednum= seednum, 
  sig_factor = sig_factor, 
  keep_raw = F,
  export_funs = export_funs_vec,
  export_packs = export_packs_vec
)



#save(results, file = file_name)



stopCluster(cl)



