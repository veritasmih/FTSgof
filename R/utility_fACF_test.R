# V_WS_quantile computes the 1-alpha qunatile of the beta * chi-squared distribution with nu
#   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#   of the test statistic V_K. This quantile is used to conduct an approximate size alpha test
#   of the hypothesis H'_0_K.
# Input: f_data = J x N functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for the test statistic V_K
#        alpha = the significance level to be used in the hypothesis test
#        M = the number of Monte-Carlo simulations for Welch-Satterthwaite approximation under the WWN assumption
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates V_K)
V_WS_quantile <- function(f_data, K, alpha=0.05, M=NULL) {
  mean_V_K <- mean_hat_V_K(f_data, K)
  var_V_K <- variance_hat_V_K(f_data, K, M=M)
  beta <- var_V_K / (2 * mean_V_K)
  nu <- 2 * (mean_V_K^2) / var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_V(f_data, K)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# t_statistic_V computes the statistic V_{T,K} = T*sum_h(||y^hat_h||^2) or h in 1:K and for T
#   inferred from the functional data f_data that is passed to the function.
# Input: f_data = J x N functional data with functions in columns
#        K = the max value in the range of time lags (1:K) used
# Output: scalar value of the statistic V_{T,K} to test the hypothesis
#         H'_{0,K} : for all h in 1:K y_h(t,s) = 0.
t_statistic_V <- function(f_data, K) {
  V_T_K <- 0
  for (h in 1:K) {
    V_T_K <- V_T_K + t_statistic_Q(f_data, h)
  }
  V_T_K
}


# mean_hat_V_K computes the approximation of the mean defined in (15) which is used in the Welch-
#   Satterthwaite approximation as mean of the chi-squared random variable approximating V_K.
# Input: f_data = J x N functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for for the test statistic V_K
# Output: scalar approximation of the mean of the test statistic V_K.
mean_hat_V_K <- function(f_data, K) {
  J <- NROW(f_data)
  sum1 <- 0
  store <- covariance_diag_store(f_data, K)
  for (i in 1:K) {
    sum1 <- sum1 + sum(store[[i]])
  }
  mu_hat_V_K <- (1 / (J^2)) * sum1
  mu_hat_V_K
}

# covariance_diag_store returns a list storage of the approximate covariances c^hat_i_i(t,s,t,s),
#   for all i in 1:K, for each encoding all values of t,s in U_J X U_J.
# Input: f_data = J x N functional data matrix with functions in columns
#        K = the maximum value of i in the range 1:K for which to compute c^hat_i_i(t,s,t,s)
# Output: a list containing K 2-D arrays encoding c^hat_i_j(t,s,t,s) evaluated at all (t,s) in
#         U_JxU_J, for i in 1:K
#
covariance_diag_store <- function(f_data, K) {
  cov_i_store <- list()
  for (j in 1:K) {
    cov_i_store[[j]] <- diagonal_covariance_i(f_data, j)
  }
  cov_i_store
}

# variance_hat_V_K computes the approximation of the variance defined in (15) which is used in
#   the Welch- Satterthwaite approximation as the variance of the chi-squared random variable
#   approximating V_K.
# Input: f_data = J x N functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for the test statistic V_K
#        M = optional argument specifying the sampling size in the related Monte Carlo method
# Output: scalar approximation of the variance of the test statistic V_K.
variance_hat_V_K <- function(f_data, K, M=NULL) {
  N <- NCOL(f_data)
  sum1 <- 0
  for (i in 1:K) {
    sum1 <- sum1 + MCint_eta_approx_i_j(f_data, i, i, M=M)
  }
  bandwidth <- ceiling(0.25 * (N ^ (1/3)))
  if (K > 1) {
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        if (abs(i-j) > bandwidth) { # empirically, past a lag of 15, error is less than 1%
          next
        }
        sum1 <- sum1 + (2 * MCint_eta_approx_i_j(f_data, i, j, M=M))
      }
    }
  }
  variance_V_K <- sum1
  variance_V_K
}



V_WS_quantile_iid <- function(f_data, K, alpha=0.05) {
  mean_V_K <- mean_hat_V_K_iid(f_data, K)
  var_V_K <- variance_hat_V_K_iid(f_data, K)
  beta <- var_V_K / (2 * mean_V_K)
  nu <- 2 * (mean_V_K^2) / var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_V(f_data, K)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# mean_hat_V_K_iid computes the approximation of the mean defined in (15) which is used in the
#   Welch-Satterthwaite approximation under the assumption that the functional data follows a
#   strong white noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for for the test statistic V_K
# Output: scalar approximation of the mean of the test statistic V_K under a strong white noise
#         assumption.
mean_hat_V_K_iid <- function(f_data, K) {
  J <- NROW(f_data)
  cov <- iid_covariance(f_data)
  mu_hat_Q_h <- K * ((sum(diag(cov)) / J)^2)
  mu_hat_Q_h
}


# variance_hat_V_K_iid computes the approximation of the variance defined in (15) which is used
#   in the Welch- Satterthwaite approximation under the assumption that the functional data
#   follows a strong white noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        K = specifies the range of lags 1:K for the test statistic V_K
# Output: scalar approximation of the variance of the test statistic V_K
variance_hat_V_K_iid <- function(f_data, K) {
  J <- NROW(f_data)
  cov_iid <- iid_covariance(f_data)
  variance_V_K_iid <- K * 2 * ( sum(cov_iid^2) / (J^2) )^2
  variance_V_K_iid
}
