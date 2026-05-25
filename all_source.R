

########################################################################################################################
## CODES FOR DATA GENERATION
#########################################################################################################################
# Arguments:
#   n          - sample size (number of observations)
#   p          - number of taxa (compositions); log-ratios are (p-1)-dimensional
#   q          - number of latent confounders (hidden factors)
#   sig_u      - measurement error SD for observed log-ratios (default 1)
#   sig_eps    - residual SD for the outcome y (default 1)
#   sig_factor - SD of idiosyncratic noise in the factor model for true_lr (default 1)
#   alpha0     - true coefficient vector for log-ratios in the outcome model;
#                if NULL, generated with nz_alp non-zero entries ~ N(0,1)
#   beta0      - true coefficient vector for confounders in the outcome model;
#                if NULL, generated as rnorm(q)
#   nz_alp     - number of non-zero entries in alpha0 when auto-generated (default 3)
#   reps       - number of replicate measurements of obs_lr for SIMEX/RC (default 3)
#
# Returns a list:
#   y            - n-vector of simulated outcomes
#   true_lr      - n x (p-1) matrix of true log-ratios 
#   obs_lr       - n x (p-1) matrix of observed log-ratios 
#   H            - n x q matrix of latent confounder values
#   MV_mat       - n x (p-1) matrix: conditional expectation E[true_lr | obs_lr]
#                  i.e. M(V) from eq. (8) in Zhao & Wang (2024)
#   Sig_mat      - (p-1) x (p-1) covariance of the signal component; = Cov_Z * A
#   cov_MV       - (p-1) x (p-1) residual covariance: Cov(true_lr - MV_mat)
#   alpha0       - true alpha used in simulation (useful when auto-generated)
#   beta0        - true beta used in simulation (useful when auto-generated)
#   obs_lr_reps  - list of length reps, each an n x (p-1) noisy replicate of true_lr

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














##########################################################################################################################
## CODES FOR DDL
##########################################################################################################################

### estimate_coefficients
### FUNCTION: using a single data set, creates unbiased estimators
###           for betahat and bhat using single value decomposition
###           and a lasso fit on trimmed X and Y for betahat then
###           a lasso fit on U and the residuals for bhat
###
### INPUT:  X, an nxp matrix; represents the measurements for each subject and covariate
###         Y, an nx1 matrix; the measured response for each subject
###         rho, a double; the trim parameter for the Q transform, default is 0.5(median)
###
### OUTPUT: a list (a) betahat, an unbiased estimator for Beta_init
###                (b) bhat, an unbiased estimator for b
estimate_coefficients <- function(X,Y,rho=0.5){
  #perform Single Val Decomp on X
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  
  #forms QX and QY(STEP 1)
  tau = quantile(D, rho)
  Dtilde = pmin(D, tau)
  Q = diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*%  t(U)
  Xtilde = Q %*% X
  Ytilde = Q %*% Y
  
  #perform a lasso fit for trimmed QX and QY(STEP 2)
  fit = glmnet::cv.glmnet(x=Xtilde, y=Ytilde)
  
  #B-init
  betahat = as.matrix((coef(fit,S =fit$lambda.min)[-1]))
  
  #determine bhat
  res = Y-X %*% betahat
  fit = glmnet::cv.glmnet(x=U %*% diag(D), y=res, alpha=1)
  bhat = t(V) %*% coef(fit, s=fit$lambda.min)[-1]
  
  return_listcoeff = list("betahat" = betahat,"bhat"=bhat)
  return(return_listcoeff)
}

### estimate_sigma
### FUNCTION: Creates an unbiased estimator of noise level based on either the default method
###           of the proposed estimator. If alt is chosen to be true, the sigmahat is found using
###           the residuals of normal lasso regression
###
### INPUT: X, an nxp matrix; represents the measurements for each subject and covariate
###        Y, an nx1 matrix; the measured response for each subject
###        rho, a double; the trim parameter for the Q transform, default is 0.5(median)
###        alt, a boolean; determines which method to use for computing sigmahat, default false
###        active_set_scaling: the correction for the size of the active set, scale estimate by n/(n-k)
###
### OUTPUT: a double, a point approximation for the true noise error in Y
estimate_sigma <- function(X,Y,rho=0.5,alt=FALSE,active_set_scaling=FALSE){
  
  #uses both fitting estimators to create an unbiased estimator of sigma
  est_coef=estimate_coefficients(X,Y,rho)
  betahat=est_coef$betahat
  bhat=est_coef$bhat
  
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  
  tau = quantile(D,rho)
  Dtilde = pmin(D, tau)
  # The length of D has length min(n,p). In order to avoid the over-shrinkage to the last n-p eigenvector
  # when n>p, the following is necessary.
  Q =diag(nrow(X))-U %*% diag(1-Dtilde / D) %*%  t(U)
  Xtilde = Q %*% X
  Ytilde = Q %*% Y
  
  #Step 7, two methods of sigmahat computation, the first has shown more robust results
  divisor=sum(diag(Q%*%Q))
  error=(norm(Ytilde-Xtilde%*%betahat,type='2'))^2
  sigmahat=(error/divisor)^0.5
  if(active_set_scaling){
    sigmahat = (nrow(X)/(nrow(X)-min(nrow(X)/2, Matrix::nnzero(betahat))))^0.5*sigmahat
  }
  if(alt){
    residuals = Y - X %*% betahat - X %*% bhat
    sigmahat = mean(residuals^2)^0.5
  }
  
  return (sigmahat)
}

### find_z
### FUNCTION: computes the normalized projection direction used to construct CI;
###           P transform must be applies to X before performing this function
###
### INPUT: X, a nxp matrix with the sample observations
###        index, an integer; the index of X(the subject) which should be used for projection
###
### OUTPUT: z, a nx1 matrix; the projection direction vector
find_z <- function(X,index){
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  #Xj and X-j
  X_j = X[,index]
  X_negj = X[,-index]
  
  #regress X-j on xj, use least min lambda to estimate gamma(Step 4)
  cvfit = glmnet::cv.glmnet(x=X_negj, y=X_j)
  gamma = coef(cvfit, s=cvfit$lambda.min)[-1]
  #eq. (8), residuals(Step 5)
  z = n^-0.5 *(X_j - X_negj %*% gamma)
  
  #variation from eq. (23) with 25% increase(read 3.6)
  V = 1.25*n^0.5*norm(z, type ="2") / (t(z) %*% X_j)
  
  #take first z whose variance is at most 25% larger than for the CV lambda
  for (lam in cvfit$glmnet.fit$lambda){
    gamma = coef(cvfit,s =lam)[-1]
    z = n^(-0.5)*(X_j - X_negj %*% matrix(gamma,p-1,1))
    if (n^0.5 * (norm(z,type="2")/(t(z) %*% X_j)) > V){
      break
    }
  }
  
  #normalize Z with 2 norm
  z = z/norm(z,type = "2")
  return(z)
}



#coverage estimation
ci.dd_lasso <- function(x, alpha = 0.05){
  se = x$se
  est_ddl = x$est_ddl
  index = x$index
  
  
  output.ci = cbind(est_ddl - qnorm(1 - alpha / 2)*se, est_ddl + qnorm(1 - alpha / 2)*se)
  
  output.ci = data.frame(cbind(index, output.ci))
  output.est = est_ddl / se
  output.pval = 2 * pnorm(-abs(output.est))
  
  output.ci$'p_value'=output.pval
  return(output.ci)
  
}



dd_lasso <- function(X,Y,index,rho=0.5,rhop=0.5){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]
  dblasso = rep(NA,length(index))
  stddev = rep(NA,length(index))
  lower = rep(NA,length(index))
  upper = rep(NA,length(index))
  est = estimate_coefficients(X,Y,rho)
  betahat = est$betahat
  for (i in seq(length(index))){
    X_negj=X[,-index[i]]
    
    #single value decomposition of X (Trim transform)
    UDV_list = svd(X_negj)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)
    
    #Reduce D, then make P trim
    Dtilde = pmin(D, quantile(D,rhop))
    P = diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P_X = P %*% X
    
    #determine projection direction then estimate betahat and bhat(the z here has been multiplied with P)
    z = find_z(P_X, index[i])
    # bhat = est$bhat
    
    # eq. (12) for a point estimation of Bj(Step 6)
    dblasso[i] = t(z) %*% P %*% (Y - X_negj %*% betahat[-index[i]]) / (t(z) %*% P %*% X[,index[i]])
    
    #eq. (23) for
    Variance = (t(z) %*% (P^2) %*% z)/(t(z) %*% P %*% X[,index[i]]) ^ 2
    
    #(Step 7)
    sigmahat = estimate_sigma(X,Y,rho)
    stddev[i] = sigmahat*Variance^.5
    
    #2-sided CI critical value for alpha
    # quantile = qnorm(1-alpha/2,mean=0,sd = 1,lower.tail=TRUE)
    # #determine bounds of interval, from eq. (13) (Step 8)
    # lower[i] = dblasso[i] - quantile * stddev[i]
    # upper[i] = dblasso[i] + quantile * stddev[i]
    # B_b=t(z) %*% P_X %*% b/(t(z) %*% P_X[,index]*stddev)
    # B_beta= t(z) %*% P_X[,-index] %*% (beta[-index]-betahat[-index])/(t(z) %*% P_X[,index]*stddev)
  }
  obj = list(index = index,
             est_init = betahat[index],
             est_ddl= dblasso,
             se = stddev)
  
  obj
}



##########################################################################################################################
#  CODES FOR APPROX. ORTHOGONALIZATION
##########################################################################################################################
# est_sigma_app <-function(y, M,rhop, ahat=NULL,plug_in=T){
#   
#   UDV_list = svd(M)
#   U = UDV_list$u
#   D = UDV_list$d
#   V = t(UDV_list$v)
#   
#   
#   tau = quantile(D,rhop)
#   Dtilde = pmin(D, tau)
#   # The length of D has length min(n,p). In order to avoid the over-shrinkage to the last n-p eigenvector
#   # when n>p, the following is necessary.
#   Q =diag(nrow(M))-U %*% diag(1-Dtilde / D) %*%  t(U)
#   Mtilde = Q %*% M
#   Ytilde = Q %*% y
#   
#   #selectiveInference::estimateSigma(Mtilde, Ytilde)$sigmahat 
#   sig_hat <- selectiveInference::estimateSigma(Mtilde, Ytilde)$sigmahat 
#   selectiveInference::estimateSigma(M, y)$sigmahat 
#   sig_hat<-       scalreg(M,y)$hsigma
#   
#   return(sig_hat)
#   
# }




app_orth <- function(y, M, rhop=0.5,k = 1,alpha=0.05 , index=NULL){
  
  
  n <- length(y) ;
  q1 <- ncol(M) ; 
  ts<-ts1 <- pval<-ahat <-se<-lower_ci <- upper_ci <-sigma_vv<-res <- rep(NA, q1)
  
  
  
  for(j in 1:q1){
    M_negj=M[,-j]
    UDV_list = svd(M_negj)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)
    
    #Reduce D, then make P trim
    Dtilde = pmin(D, quantile(D,rhop))
    P = diag(nrow(M_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P_M = P%*%M
    P_y  <- P %*% y
    
    
    
    
    # the proposed projection direction in Battey and Reid (2023)
    proj_vec <- MASS::ginv(P)%*%solve(k*diag(n) + tcrossprod(M[,-j]))%*%M[,j] # dimension : n by 1
    
    
    # compute the test statistic which asymptotically follows a standard normal dist'n
    ahat[j] <-sum(proj_vec*P_y)/sum(proj_vec*P_M[, j])
    
    
    # compute sigma_vv from Equation (5)
    sigma_vv<- (t(proj_vec) %*% P%*%P%*%proj_vec) /(t(proj_vec)%*%P_M[,j])^2
    
    #residuals
    r_j <- as.numeric(P_y - ahat[j] * P_M[, j])
    
    # estimate sigma^2 from this r_j using the formula in the notes:
    #   \hat sigma_j^2 = (1/(n-1)) sum_i r_{j,i}^2
    sigma2_j <- sum(r_j^2) / (n-1)
    
    # Standard error
    se[j] <- sqrt(sigma2_j*sigma_vv)
    
    
    
    ts[j] <- ahat[j]/se[j]
    
    
    # compute the p-value
    pval[j] <- 2*pnorm(abs(ts[j]), lower.tail = F)
    
    z_critical <- qnorm(1 - alpha/2)
    
    # Confidence interval
    lower_ci[j] <- ahat[j] - z_critical * se[j]
    upper_ci[j] <- ahat[j] + z_critical * se[j]
  }
  
  #covered<-(alpha0[1:7] >= lower_ci[1:7]) & (alpha0[1:7] <= upper_ci[1:7])
  #cbind(Alpha=alpha0[1:7], ahat=ahat[1:7],L=lower_ci[1:7], U=upper_ci[1:7], Covered=covered,Pval= pval[1:7])
  
  
  return(list(ts = ts, pval = pval, ahat=ahat, se=se,lower_ci = lower_ci, upper_ci = upper_ci))
  #return(list(ahat=ahat, se=se,lower_ci = lower_ci, upper_ci = upper_ci))
}




estimate_logratio_ME <- function(V_array, V_lr,center_barV = TRUE, eps=1e-6) {
  
  if (length(dim(V_array)) != 3L) {
    stop("V_array must be a 3D array: n x p1 x R.")
  }
  
  n  <- dim(V_array)[1L]
  p1 <- dim(V_array)[2L]   # p-1
  R  <- dim(V_array)[3L]
  
  ## 1. Subject-level mean log-ratios: bar{V}_i
  # barV[i, j] = average over replicates r of V_array[i, j, r]
  barV <- apply(V_array, c(1, 2), mean)   # n x p1
  
  
  ## 2. Estimate mu_ztilde = E[Ztilde]
  mu_ztilde_hat <- colMeans(barV)         # length p1
  
  
  ## 3. Estimate sigma_u^2 using within-subject differences
  Sigma_within <- matrix(0, p1, p1)
  
  for (i in seq_len(n)) {
    Vi <- V_array[i, , , drop = FALSE]    # 1 x p1 x R
    Vi <- matrix(Vi, nrow = p1, ncol = R) # p1 x R
    for (r in 1:(R - 1)) {
      for (s in (r + 1):R) {
        d_rs <- Vi[, r] - Vi[, s]         # length p1
        Sigma_within <- Sigma_within + d_rs %*% t(d_rs)
      }
    }
  }
  
  num_pairs <- R * (R - 1) / 2
  Sigma_within <- Sigma_within / (n * num_pairs)
  
  # E[Sigma_within] = 4 * sigma_u^2 * I_p1  =>  use trace to solve
  sigma_u2_hat <- sum(diag(Sigma_within)) / (4 * p1)
  
  ## 4. Estimate Cov(barV) with chosen method
  
  # Optionally center barV for covariance estimation
  barV_for_cov <- if (center_barV) {
    scale(barV, center = TRUE, scale = FALSE)
  } else {
    barV
  }
  
  
  #require  library(corpcor)
  Sigma_barV_hat <- corpcor::cov.shrink(barV_for_cov)
  
  
  ## 5. De-noise Cov(barV) to get Sigma_ztilde
  # Cov(barV) = Sigma_ztilde + (2 sigma_u^2 / R) * I  => subtract
  Sigma_ztilde_hat <- Sigma_barV_hat - (2 * sigma_u2_hat / R) * diag(p1)
  
  
  A <- Sigma_ztilde_hat %*% solve(Sigma_ztilde_hat + 2 * sigma_u2_hat * diag(p1))
  
  ## 6. Compute mu(V_i) for each subject
  M_V_hat <- matrix(NA, n, p1)
  for (i in 1:n) {
    diff_i <- V_lr[i, ] - mu_ztilde_hat     # length p-1
    M_V_hat[i, ] <- mu_ztilde_hat + A %*% diff_i
  }
  
  ## 7. Return everything
  list(
    mu_ztilde   = as.numeric(mu_ztilde_hat),
    sigma_u2    = as.numeric(sigma_u2_hat),
    Sigma_ztilde = Sigma_ztilde_hat,
    barV        = barV,
    MV_hat=M_V_hat
  )
}



#summary function

summarize_simulation <- function(sim_results) {
  # sim_results: list of data.frames, each with columns like
  # index, alpha_true, estimate, std_error, covered (0/1)
  
  # combine all results into one long data frame
  df <- do.call(rbind, sim_results)
  
  # compute per-index summary
  summary <- df %>%
    dplyr::group_by(index) %>%
    dplyr::summarise(
      alpha_true   = unique(alpha_true),
      mean_est     = mean(estimate, na.rm=TRUE),
      coverage_pct = mean(covered, na.rm = TRUE),
      
      # bias and variability of the estimator
      avg_bias   = mean(estimate - alpha_true, na.rm = TRUE),
      emp_sd     = sd(estimate, na.rm = TRUE),
      
      # standard error behaviour
      avg_se     = mean(std_error, na.rm = TRUE),
      sd_se      = sd(std_error, na.rm = TRUE),
      
      
      reps       = dplyr::n(),
      .groups    = "drop"
    )
  
  
  
  summary
  
}




##########################################################################################################################
## SIMULATION FUNCTIONS
##########################################################################################################################


# sim_fn: Single simulation replicate
#
# Generates data via dgen_fn, then fits three estimators:
#   1. DD-Lasso       (ddl)
#   2. Approximate orthogonalization (approximate orthogonalization)
#   3. DR-C-Lasso from Zhao & Wang / biometric paper (biom)
#
# Arguments:
#   n, p, q        - sample size, taxa count, confounder count
#   alpha0, beta0  - true coefficient vectors (passed through to dgen_fn)
#   sig_u          - measurement error SD
#   sig_eps        - outcome residual SD
#   sig_factor     - factor model noise SD
#   nz_alp         - number of non-zero entries in alpha0
#   reps           - number of obs_lr replicates
#   seednum        - random seed for reproducibility
#
# Returns a named list: app, biom, ddl — each a data frame with columns:
#   index, alpha_true, estimate, std_error, p_value, lower_ci, upper_ci, covered (coverage probability)


sim_fn <- function(n, p, q, alpha0, beta0, sig_u = 1, sig_eps = 1, sig_factor = 1, nz_alp = 3, reps=3, seednum){
  
  set.seed(seednum)
  
  dd <- dgen_fn(
    n = n,
    p = p, 
    q = q, 
    alpha0 = alpha0, 
    beta0 = beta0, 
    sig_u = sig_u, 
    sig_eps = sig_eps, 
    sig_factor = sig_factor,
    nz_alp = nz_alp,
    reps=reps
  )
  
  
  y <- dd$y ; MV <- dd$MV_mat
  
  
  
  #ddl
  ddlasso.obj <-dd_lasso(X=MV, Y=y, index=c(1:10), rho=0.5, rhop=0.5)
  ci.obj<-ci.dd_lasso(ddlasso.obj)
  colnames(ci.obj)<-c('index', 'lower_ci', 'upper_ci','p_value')
  covered <- (alpha0[1:10] >= ci.obj$lower_ci) & (alpha0[1:10] <= ci.obj$upper_ci)
  ci.obj$covered <-covered
  ci.obj$alpha_true <-alpha0[1:10]
  ci.obj$estimate <-ddlasso.obj$est_ddl
  ci.obj$std_error <-ddlasso.obj$se
  #
  #
  results_ddl<-ci.obj
  
  
  
  # approx. orth..
  app_orth.obj<-app_orth(y=y,M= MV,rhop = .5)
  results_app<-data.frame(
    index = 1:10,
    alpha_true = alpha0[1:10],
    estimate = app_orth.obj$ahat[1:10],
    std_error = app_orth.obj$se[1:10],
    p_value = app_orth.obj$pval[1:10],
    lower_ci = app_orth.obj$lower_ci[1:10],
    upper_ci = app_orth.obj$upper_ci[1:10],
    covered = (alpha0[1:10] >= app_orth.obj$lower_ci[1:10]) & (alpha0[1:10] <= app_orth.obj$upper_ci[1:10])
  )
  
  
  
  #biometric paper...
  y = y - mean(y)
  MV <- scale(dd$MV_mat, scale = F)
  
  rclasso_model <- cv.glmnet(MV, y, alpha = 1, nfolds = 5, family = "gaussian", intercept = F)
  alpha_hat = as.vector(coef(rclasso_model, s = rclasso_model$lambda.min))[-1]
  Sinv_mat <- solve(dd$Sig_mat)
  alpha_drclasso = alpha_hat + (1/n)*tcrossprod(Sinv_mat, dd$MV_mat)%*%(y - MV%*%alpha_hat)
  
  eg <- eigen(dd$Sig_mat)
  Sigma_inv_half <- eg$vectors %*% diag(1/sqrt(eg$values)) %*% t(eg$vectors)
  fit = scalreg(MV,y) # scaled lasso estimator for sigma
  cov_drclasso = (fit$hsigma^2/n)*Sigma_inv_half%*%(crossprod(dd$MV_mat)/n)%*%Sigma_inv_half
  se_drclasso = sqrt(diag(cov_drclasso)) 
  
  
  # Confidence interval
  lower_ci <- alpha_drclasso - qnorm(1-0.05/2) * se_drclasso
  upper_ci <- alpha_drclasso + qnorm(1-0.05/2) * se_drclasso
  
  cover_vec <- as.vector((dd$alpha0 >= lower_ci) & (dd$alpha0 <= upper_ci))
  p_drclasso <- as.vector(2*pnorm(abs(alpha_drclasso)/se_drclasso, lower.tail = F))
  
  results_biom<-data.frame(
    index = 1:10,
    alpha_true = alpha0[1:10],
    estimate =  alpha_drclasso[1:10],
    std_error = se_drclasso[1:10],
    p_value =    p_drclasso[1:10],
    lower_ci =   lower_ci[1:10],
    upper_ci =   upper_ci[1:10],
    covered = cover_vec[1:10]
  )
  
  list(app=results_app, biom=results_biom, ddl=results_ddl)
  
  
}





###########################################################################################################################
# PARALLEL SIMULATION SETUP
###########################################################################################################################


# ─────────────────────────────────────────────────────────────────────────────
# Exports for parallel workers
# ─────────────────────────────────────────────────────────────────────────────

# Functions that must be exported to each foreach worker
export_funs_vec <- c(
  "sim_fn",
  "dgen_fn",
  "estimate_coefficients",
  "estimate_sigma",
  "find_z",
  "dd_lasso",
  "ci.dd_lasso",
  "app_orth"
)

# Packages that must be loaded on each foreach worker
export_packs_vec <- c(
  "scalreg",
  "glmnet",
  "MASS",
  "dplyr",
  "corpcor"
)


# Function for parrallel simulation setup
# Arguments:
#   sample_size, taxa_size, conf_size - n, p, q passed to sim_fn
#   alpha_truth, beta_truth           - true coefficient vectors
#   nsim                              - number of iterations
#   seednum                           - base seed; replicate i uses seednum + i
#   sig_u, sig_eps, sig_factor        - noise parameters forwarded to sim_fn
#   nz_alp                            - sparsity of alpha_truth
#   reps                              - obs_lr replicates per sim_fn call passed to data generation functon
#   keep_raw                          - if TRUE, return raw list of nsim results
#                                       instead of summarized output
#   export_packs, export_funs         - packages/functions forwarded to workers
#
# Returns:
#   If keep_raw = TRUE  : raw list of length nsim, each element as from sim_fn
#   If keep_raw = FALSE : named list (approx. ortho, biometric paper, ddl) of summarized results



sim_run <-function(sample_size, taxa_size, conf_size, alpha_truth, 
                   beta_truth, nsim, seednum,sig_u = 1, sig_eps = 1, 
                   sig_factor = 1, nz_alp = 3, reps=3, keep_raw=F, 
                   export_packs, export_funs){
  
  all_res <- foreach(i = 1:nsim,
                     .packages = export_packs,
                     .export =export_funs) %dopar% {
                       
                       seed_num = seednum + i   
                       sim_fn(n=sample_size, 
                              p=taxa_size, 
                              q=conf_size, 
                              alpha0= alpha_truth, 
                              beta0= beta_truth,
                              sig_u = sig_u,
                              sig_eps = sig_eps,
                              sig_factor = sig_factor,
                              nz_alp = nz_alp,
                              reps   =reps,
                              seednum=seed_num
                       )
                     } 
  
  
  if(keep_raw) return(all_res)
  
  
  all_data_app <- lapply(all_res, function(x) x$app)
  all_data_biom <- lapply(all_res, function(x) x$biom)
  all_data_ddl <- lapply(all_res, function(x){ x$ddl})
  
  
  #approx orthogo.
  summary_res_app <-summarize_simulation(all_data_app)
  
  #biometric paper
  summary_res_biom <-summarize_simulation(all_data_biom)
  
  #ddl
  summary_res_ddl <-summarize_simulation(all_data_ddl)
  
  list(
    app=summary_res_app, 
    biom=  summary_res_biom, 
    ddl=summary_res_ddl
  )
  
  
}







