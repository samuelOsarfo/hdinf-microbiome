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




#App_Orth
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

