library(MASS)


# Define two functions to generate covariance matrices for comparison
# n = sample_size
# p = taxa size

corfactor_fn <-function(n, p, q, sig_factor = 1){
  H   <- matrix(rnorm(n * q), nrow = n, ncol = q) # confounders (factors)
  Psi <- matrix(rnorm((p - 1) * q), nrow = p - 1, ncol = q) # matrix of loadings
  true_lr <- tcrossprod(H, Psi) + matrix(rnorm(n * (p - 1), sd = sig_factor), nrow = n, ncol = p - 1) 
  
  
  return(cor(true_lr))
}


corBiom_fn <-function(n,p){
  
  mu  <- c(rep(log(round(p/2)), 5), rep(0, p-5))
  
  sgma <-toeplitz(0.2^(0: (p-1)))
  
  true_lx <- mvrnorm(n, mu, sgma)
  
  true_lr <-true_lx[, -p] - true_lx[,p]
  
  
  return(cor(true_lr))
  
}




set.seed(123)
n_sim <- 100      # number of simulation replicates
n <- 300          # sample size
p <- 300           # number of taxa
q <- 3            # number of factors
sig_factor <- 0.6  # noise sd in factor model (can experiment)

mean_abs_cor_factor <- numeric(n_sim)
mean_abs_cor_biom <- numeric(n_sim)

for (i in 1:n_sim) {
  cor_factor <- corfactor_fn(n, p, q, sig_factor)
  cor_biom   <- corBiom_fn(n, p)
  
  # Mean absolute off-diagonal correlation (excluding diagonal)
  off_diag <- function(cor_mat) {
    cor_mat[upper.tri(cor_mat)]
  }
  mean_abs_cor_factor[i] <- mean(abs(off_diag(cor_factor)))
  mean_abs_cor_biom[i]   <- mean(abs(off_diag(cor_biom)))
}

# Compare distributions
summary(mean_abs_cor_factor)
summary(mean_abs_cor_biom)


# Proportion of replicates where factor model has higher mean absolute correlation
mean(mean_abs_cor_factor > mean_abs_cor_biom)