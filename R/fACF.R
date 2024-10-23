#' Functional Autocorrelation Function (fACF) Plot
#'
#' @description This function provides a graphical summary of the fACF of a functional time series (FTS) across different time lags \eqn{h = 1:H}.
#' It also plots the \eqn{100 (1-\alpha)\%} confidence bounds, developed under both weak white noise (WWN) and strong white noise (SWN) assumptions for all lags \eqn{h = 1:H}.
#'
#' @param f_data A functional object composed by \eqn{N} functional time series, given \eqn{J} grid points in the domains; or a \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param H A positive integer value. The maximum lag for which to compute the coefficients and confidence bounds.
#' @param alpha A numeric value between 0 and 1 specifying the significance level to be used for the confidence bounds.
#' @param wwn_bound A Boolean value allowing the user to turn on the WWN bound. FALSE by default. Speeds down computation when TRUE.
#' @param M A positive integer value. The number of Monte-Carlo simulations used to compute the confidence bounds under the WWN assumption. 
#' If \eqn{M = NULL, M = \text{floor}((\max(150 - N, 0) + \max(100 - J, 0) + (J / \sqrt{2})))},
#' ensuring that the number of Monte Carlo simulations is adequate based on the dataset size.
#' @param J A positive integer value. Evaluate the functional objects at a pre-specified number of points on a grid J. The default value is NULL, indicating no change on the number of grid points.
#' 
#' @details This function computes and plots functional autocorrelation coefficients at lag \eqn{h}, for \eqn{h \in 1:H}. Given functional observations, \eqn{X_1,\ldots, X_N}, the sample autocovariance kernel at lag \eqn{h} can be computed by
#' \deqn{
#' \hat{\gamma}_{N,h}(t,s)=\frac{1}{N}\sum_{i=1}^{N-h} (X_i(t)-\bar{X}_N(t))(X_{i+h}(s)-\bar{X}_N(s)),\ \ \ \ 0 \le h < N,
#' }
#' where \eqn{\bar{X}_N(t) = \frac{1}{N} \sum_{i=1}^N X_i(t)}. Then, the fACF at lag \eqn{h} is defined by measuring
#' the magnitude (\eqn{L^2}-norm) of the lagged autocovariance kernel \eqn{\hat\gamma_{N,h}}:
#' \deqn{
#' \hat\rho_h =\frac{\|\hat{\gamma}_{N,h}\|}{\int \hat{\gamma}_{N,0}(t,t)dt}, \ \ \ \ \|\hat{\gamma}_{N,h}\|^2=\iint  \hat{\gamma}_{N,h}^2(t,s) dtds.
#' }
#' This function plots estimated asymptotic \eqn{100 (1-\alpha)\%} confidence bounds under the WWN assumption.
#' Additionally, it computes similar (constant) bounds under the SWN assumption.
#'
#'
#' @return Plot of the estimated functional autocorrelation coefficients for lags \eqn{h \in 1:H} with the WWN
#' \eqn{100 (1-\alpha)\%} upper confidence bound for each lag, as well as the constant SWN
#' \eqn{100 (1-\alpha)\%} upper confidence bound.
#'
#' @references
#' [1] Kokoszka P., Rice G., Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' [2] Mestre G., Portela J., Rice G., Roque A. M. S., Alonso E. (2021). Functional time series model identification
#'  and diagnosis by means of auto-and partial autocorrelation analysis. Computational Statistics & Data Analysis, 155, 107108.
#'
#' @examples
#' \donttest{
#' data(Spanish_elec) # Daily Spanish electricity price profiles
#' fACF(Spanish_elec)
#' fACF(Spanish_elec, H=10, wwn_bound=TRUE)
#' }
#'
#' @export
#' @import stats
#'
fACF <- function(f_data, H=20, alpha=0.05, wwn_bound=FALSE, M=NULL, J=NULL) {

  if ((H < 1) | (H %% 1 != 0)) {
    stop("The parameter 'H' must be a positive integer.")
  }
  if ((alpha > 1) | (alpha < 0)) {
    stop("The 'alpha' parameter must be a value between 0 and 1.")
  }
  
  data_class <- class(f_data)[[1]]
  if (data_class=="funData" ) {
    
    if (is.null(J) ){
      f_data <- t(f_data@X)
    } else {
      f_data <- t(f_data@X)
      J_raw <- NROW(f_data)
      basis <- create.bspline.basis(rangeval = c(0,1), nbasis = 25)
      fd_data <- smooth.basis(1:J_raw/J_raw, y = f_data, fdParobj = basis)$fd
      f_data <- eval.fd(fd_data, 1:J/J)
    }
    
  } else if (data_class=="matrix") {
  
    if (is.null(J) ){
      f_data <- f_data
    } else {
      J_raw <- NROW(f_data)
      basis <- create.bspline.basis(rangeval = c(0,1), nbasis = 25)
      fd_data <- smooth.basis(1:J_raw/J_raw, y = f_data, fdParobj = basis)$fd
      f_data <- eval.fd(fd_data, 1:J/J)
    }
    
  } else if (data_class=="fd" ) {
    
    if (is.null(J) ){
      tempFun <- fd2funData(f_data, argvals = seq(0, 1, length.out = length(f_data[["fdnames"]][["time"]])  ) )
      f_data <- t(tempFun@X) # realization can be different from the original discrete data given the smoothing
    } else {
      tempFun <- fd2funData(f_data, argvals = seq(0, 1, length.out = J)  ) 
      f_data <- t(tempFun@X)
    }
    
  } else {
    stop("The input must be either a matrix or a funData object.")
  }

  
  J = NROW(f_data)
  coefficients = array(0, H)
  B_iid_bounds = array(0,H)
  lags = 1:H
  if (wwn_bound == TRUE) {
    B_h_bounds = array(0,H)
    for (h in lags){
      coefficients[h] <- autocorrelation_coeff_h(f_data, h)
      B_h_bounds[h] <- B_h_bound(f_data, h, M=M)
    }
  } else {
    for (h in lags){
      coefficients[h] <- autocorrelation_coeff_h(f_data, h)
    }
  }
  m<-max(coefficients, B_iid_bounds)
  #par(mar = c(4, 4, 2, 1))
  plot(lags, coefficients, ylim=c(0,m+0.3),
       type='h', xlab='Lag', ylab='Coefficient', main = 'fACF plot')
  lines(rep(B_iid_bound(f_data), H), col='red', lty='solid')
  if (wwn_bound == TRUE) {
    lines(B_h_bounds, col='blue', lty='dotted')
    legend("topleft",
           legend=c('Estimated Autocorrelation Coefficients',
                    paste('WWN ', (1 - alpha) * 100,'% Confidence Bound', sep=' '),
                    paste('SWN ', (1 - alpha) * 100,'% Confidence Bound', sep=' ')),
           col=c('black', 'blue', 'red'), lty=c('solid', 'dotted', 'solid')
           , cex=0.9, bty = "n",y.intersp = 1)
  } else {
    legend('topleft',
           legend=c('Estimated Autocorrelation Coefficients', paste('SWN ', (1 - alpha) * 100,'% Confidence Bound', sep=' ')),
           col=c('black', 'red'), lty=c('solid', 'solid'), cex=0.9, bty = "n",y.intersp = 1)
  }
}


