# Store all utility functions used in the 'fACF' function

# autocorrelation_coff_h computes the approximate functional autocorrelation coefficient rho^hat_h at lag h
# Input: f_data = J x N functional data matrix with functions in columns
#        lag = lag for which to compute the coefficient
# Output: scalar value of the approximate functional autocorrelation coefficient at lag h.
#
autocorrelation_coeff_h <- function(f_data, lag) {
  N <- NCOL(f_data)
  num <- sqrt(t_statistic_Q(f_data, lag))
  denom <- sqrt(N) * diagonal_autocov_approx_0(f_data)
  coefficient <- num / denom
  coefficient
}

# t_statistic_Q computes the test statistic Q_{T,h} = T*||y^hat_h||^2 for fixed h and for T
# Input: f_data = J x N functional data matrix with observed functions in columns
#        lag = the fixed time lag used in the computation of the statistic
# Output: scalar value of the statistic Q_{T,h} to test the hypothesis H_{0,h} : y_h(t,s) = 0.
t_statistic_Q <- function(f_data, lag) {
  N = NCOL(f_data)
  J = NROW(f_data)
  gamma_hat <- autocov_approx_h(f_data, lag)
  Q_T_h <- N * sum(gamma_hat^2) / (J^2)
  Q_T_h
}

# diagonal_autocov_approx_0 computes the intergral of y^hat_0(t,t) with respect to \mu(dt).
# Input: f_data = J x M functional data matrix with functions in columns
# Output: scalar value of the integral of y^hat_0(t,t) with respect to \mu(dt).
#
diagonal_autocov_approx_0 <- function(f_data) {
  J <- NROW(f_data)
  gamma_hat_0 <- autocov_approx_h(f_data, 0)
  sum(diag(gamma_hat_0)) / J
}

# autocov_approx_h computes the approximate autocovariance gamma^hat_h(t,s) for a given lag (h)
#   for every (t,s) in U_J X U_J.
# Input: f_data = J x M functional data matrix with observed functions in columns
#        lag = the fixed lag for which to compute gamma^hat_h(t,s)
# Output: a 2-D array encoding the values of gamma^hat_h(t,s) for every (t,s) in U_J X U_J.
#
autocov_approx_h <- function(f_data, lag) {
  N = NCOL(f_data)
  c_f_data <- center(f_data)
  gamma_hat <- c_f_data[,(1+lag):N]%*%t(c_f_data[,1:(N-lag)])/ N
  gamma_hat
}

# center centers the functional data f_data by substracting the row means from the data.
# Input: f_data = J x N functional data matrix with observed functions in columns
# Output: a matrix containing the centered functional data
center <- function(f_data) {
  f_data - rowMeans(f_data)
}

# B_h_bound returns an approximate asymptotic upper 1-alpha confidence bound for the functional
#   autocorrelation coefficient at lag h under the assumption that f_data forms a weak white
#   noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        lag = the lag for which to ccmpute the bound
#        alpha = significance level of the bound
#        M = A positive Integer value. The number of Monte-Carlo simulations used under the WWN assumption.
# Output: scalar value of the 1-alpha confidence bound for the functional autocorrelation
#         coefficient at lag h under a weak white noise assumption.
#
B_h_bound <- function(f_data, lag, alpha=0.05, M=NULL) {
  N <- NCOL(f_data)
  quantile = Q_WS_quantile(f_data, lag, alpha=alpha, M=M)$quantile
  num <- sqrt(quantile)
  denom <- sqrt(N) * diagonal_autocov_approx_0(f_data)
  bound <- num / denom
  bound
}


# Q_WS_quantile computes the 1-alpha quantile of the beta * chi-squared distribution with nu
#   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#   of the test statistic Q_h. This quantile is used to conduct an approximate size alpha test
#   of the hypothesis H_0_h under the assumption that the data follows a weak  white noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        lag = specifies the lag used for the test statistic Q_h
#        alpha = the significance level to be used in the hypothesis test
#        M = A positive Integer value. The number of Monte-Carlo simulations used under the WWN assumption.
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates Q_h).
Q_WS_quantile <- function(f_data, lag, alpha=0.05, M=NULL) {
  mean_Q_h <- mean_hat_Q_h(f_data, lag)
  var_Q_h <- variance_hat_Q_h(f_data, lag, M=M)
  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_Q(f_data, lag)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# mean_hat_Q_h computes the approximation ofthe mean of the test statistic Q_h which is used in the
# Welch-Satterthwaite approximation as mean of the chi-squared random variable approximating Q_h.
# Input: f_data = J x N functional data matrix with functions in columns
#        lag = specifies the lag used in the test statistic Q_h
# Output: scalar approximation of the mean of the test statistic Q_h.
mean_hat_Q_h <- function(f_data, lag) {
  J <- NROW(f_data)
  cov <- diagonal_covariance_i(f_data, lag)
  mu_hat_Q_h <- (1 / (J^2)) * sum(cov)
  mu_hat_Q_h
}

# diagonal_covariance_i returns the approximate covariance c^hat_i_i(t,s,t,s), encoding all
#   values of t,s in U_J X U_J, i in 1:T.
# Input: f_data = the functional data matrix with functions in columns
#        i = the index i in 1:T that we are computing the covariance diagonal for
# Output: a 2-D array encoding c^hat_i_j(t,s,t,s) evaluated at all (t,s) in U_JxU_J.
diagonal_covariance_i <- function(f_data, i) {
  N = NCOL(f_data)
  J = NROW(f_data)
  c_f_data <- center(f_data)
  sum1 <- array(0, c(J, J))
  for (k in (1+i):N) {
    sum1 <- sum1 + ((c_f_data[,k-i])^2 %o% (c_f_data[,k])^2)
  }
  cov <- (1 / N) * sum1
  cov
}

# variance_hat_Q_h computes the approximation of the variance of the test statistic Q_h which is used in
#   the Welch-Satterthwaite approximation as variance of the chi-squared random variable
#   approximating Q_h.
# Input: f_data = J x N functional data matrix with functions in columns
#        lag = specifies the lag used in the test statistic Q_h (lag = h in paper)
#        M = optional argument specifying the sampling size in the related Monte Carlo method
# Output: scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h <- function(f_data, lag, M=NULL) {
  variance_Q_h <- MCint_eta_approx_i_j(f_data, lag, lag, M=M)
  variance_Q_h
}

# MCint_eta_approx_i_j computes an approximation of eta_i_j using the second
#   Monte Carlo integration method "MCint" defined on page 8.
# Input: f_data = J x M functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing eta^hat_i_j for
#        M = number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
# Output: scalar value of eta^_hat_i_j computed using the MCint method.
MCint_eta_approx_i_j <- function(f_data, i, j, M=NULL) {
  J <- NROW(f_data)
  T <- NCOL(f_data)
  if (is.null(M)) {
    M = floor((max(150 - T, 0) + max(100-J,0) + (J / sqrt(2))))
  }
  rand_samp_mat <- matrix(nrow=M, ncol=4)
  rand_samp_mat <- cbind(sample(1:J, M, replace = TRUE),sample(1:J, M, replace = TRUE),sample(1:J, M, replace = TRUE),sample(1:J, M, replace = TRUE))

  eta_hat_i_j_sum <- 0
  for (k in 1:M) {
    cov <- scalar_covariance_i_j(f_data, i, j, rand_samp_mat[k,])
    eta_hat_i_j_sum <- eta_hat_i_j_sum + (cov^2)
  }
  eta_hat_i_j <- (2/M) * eta_hat_i_j_sum
  eta_hat_i_j
}

# scalar_covariance_i_j returns the approximate covariance c^hat_i_j(t,s,u,v) evaluated at a
#   given t,s,u,v in U_J X U_J X U_J X U_J (for use in MCint method).
# Input: f_data = J x M functional data matrix with functions in columns
#        i,j = the indices i,j in 1:T that we are computing the covariance for
#        times = a 4-element vector representing the values (t,s,u,v)
# Output: scalar value of the computed covariance c^hat_i_j(t,s,u,v).
#
scalar_covariance_i_j <- function(f_data, i, j, times) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum1 <- 0
  for (k in (1+max(i,j)):N) {
    sum1 <- sum1 + c_f_data[times[1],k-i] * c_f_data[times[2],k] * c_f_data[times[3],k-j]  *
      c_f_data[times[4],k]
  }
  cov <- (1/N) * sum1
  cov
}

# B_iid_bound returns an approximate asymptotic upper 1-alpha confidence bound for the functional
#   autocorrelation coefficient at lag h under the assumption that f_data forms a strong
#   white noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        alpha = significance level of the bound
# Output: scalar value of the 1-alpha confidence bound for the functional autocorrelation
#         coefficient at lag h under a strong white noise assumption.
#
B_iid_bound <- function(f_data, alpha=0.05) {
  N <- NCOL(f_data)
  quantile_iid = Q_WS_quantile_iid(f_data, alpha=alpha)$quantile
  num <- sqrt(quantile_iid)
  denom <- sqrt(N) * diagonal_autocov_approx_0(f_data)
  bound <- num / denom
  bound
}

# Q_WS_quantile_iid computes the size alpha test of the hypothesis H_0_h using the Welch-Satterthwaite
#   approximation under the assumption that the data follows a strong white noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        alpha = the significance level to be used in the hypothesis test
# Output: scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#         degrees of freedom (which approximates Q_h) (computed under a strong white noise
#         assumption).
Q_WS_quantile_iid <- function(f_data, alpha=0.05) {
  mean_Q_h <- mean_hat_Q_h_iid(f_data)
  var_Q_h <- variance_hat_Q_h_iid(f_data)
  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * qchisq(1 - alpha, nu)
  statistic <- t_statistic_Q(f_data, lag = 1)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}

# mean_hat_Q_h_iid computes the approximation of the test statistic Q_h which is used in the
#   Welch-Satterthwaite approximation under the assumption that the functional data follows a
#   strong white noise.
# Input: f_data = J x N functional data matrix with functions in columns
# Output: scalar approximation of the mean of the test statistic Q_h under a strong white noise
#         assumption.
mean_hat_Q_h_iid <- function(f_data) {
  J <- NROW(f_data)
  cov <- iid_covariance(f_data)
  mu_hat_Q_h_iid <- (sum(diag(cov)) / J)^2
  mu_hat_Q_h_iid
}

# iid_covariance returns one of the two independent sum terms in the approximate covariance
#   c^*_0(t,s,u,v) definition, encoding all values of (t,s) in U_J X U_J, i,j in 1:T.
# Input: f_data = J x N functional data matrix with functions in columns
# Output: returns a 2-D tensor of c^*(t,s), one of the two independent sums in the computation
#         of c^*(t,s,u,v).
#
iid_covariance <- function(f_data) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum1 <- 0
  for (i in 1:N) {
    sum1 <- sum1 + c_f_data[,i] %o% c_f_data[,i]
  }
  sum1 / N
}

# variance_hat_Q_h_iid computes the approximationof the test statistic Q_h which is used
#   in the Welch- Satterthwaite approximation under the assumption that the functional data
#   follows a strong white noise.
# Input: f_data = J x N functional data matrix with functions in columns
#        lag = specifies the lag used in the test statistic Q_h (lag = h in paper)
# Output: scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h_iid <- function(f_data) {
  J <- NROW(f_data)
  cov_iid <- iid_covariance(f_data)
  variance_Q_h_iid <-2 * ( sum(cov_iid^2) / (J^2) )^2
  variance_Q_h_iid
}
