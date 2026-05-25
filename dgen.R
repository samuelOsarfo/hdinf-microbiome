 

dgen_fn <- function(n, p, q, sig_u = 1, sig_eps = 1, sig_factor = 1, 
                    alpha0 = NULL, beta0 = NULL, nz_alp = 3, reps=3){
  
  # 1. simulate the log-ratio of true compositions from the factor model with hidden confounders
  H   <- matrix(rnorm(n * q), nrow = n, ncol = q) # confounders (factors)
  Psi <- matrix(rnorm((p - 1) * q), nrow = p - 1, ncol = q) # matrix of loadings
  true_lr <- tcrossprod(H, Psi) + matrix(rnorm(n * (p - 1), sd = sig_factor), nrow = n, ncol = p - 1) 
  
  # 2. simulate y using both true_lr and H
  if(is.null(alpha0)){
    alpha0 <- c(rnorm(nz_alp), rep(0, p - 1 - nz_alp)) #alpha0 <- runif(p - 1)
  }
  
  if(is.null(beta0)){
    beta0 <- rnorm(q) #beta0 <- runif(p - 1)
  }
  
  y <- true_lr %*% alpha0 + H%*%beta0 + rnorm(n, sd = sig_eps)
  
  # 3. generate the observed log-ratio from the model (4) in Zhao and Wang (2024)
  obs_lr <- true_lr + matrix(rnorm(n * (p - 1), sd = sqrt(2*sig_u^2)), nrow = n, ncol = p - 1)
  
  
  #replicates for measurement error
  obs_lr_reps <- vector("list", reps)
  for (r in 1:reps) {
    obs_lr_reps[[r]] <- true_lr + matrix(rnorm(n * (p - 1), sd = sqrt(2*sig_u^2)), nrow = n, ncol = p - 1)
  }
  
  
  # 4. compute the true conditional expectation (denoted as M(V) in (8) in Zhao and Wang)
  #mu_ztilde <- 0
  cov_ztilde <- tcrossprod(Psi) + diag(sig_factor^2, nrow = p - 1, ncol = p - 1)
  A     <- solve(cov_ztilde + diag(2 * sig_u^2, nrow = p - 1, ncol = p - 1))%*%cov_ztilde
  MV_mat <- obs_lr%*%A
  cov_MV <- cov_ztilde %*% (diag(1, p - 1) - A)
  Sig_mat <- cov_ztilde%*%A
  
  return(list(y = y, true_lr = true_lr, obs_lr = obs_lr, H = H, MV_mat = MV_mat, Sig_mat = Sig_mat, 
              cov_MV = cov_MV, alpha0 = alpha0, beta0 = beta0, obs_lr_reps=obs_lr_reps))
}


