#' Functional Spherical Autocorrelation Function (fSACF) Plot
#'
#' @description This function offers a graphical summary of the fSACF of a functional time series (FTS) across different time lags \eqn{h = 1:H}.
#' It also plots \eqn{100 \times (1-\alpha)\%}  confidence bounds developed under strong white noise (SWN) assumption for all lags \eqn{h = 1:H}.
#'
#' @param f_data A functional object composed by \eqn{N} functional time series, given \eqn{J} grid points in the domains; or a \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param H A positive integer value. The maximum lag for which to compute the coefficients and confidence bounds.
#' @param alpha A numeric value between 0 and 1 specifying the significance level to be used for the confidence bounds.
#' @param J A positive integer value. Evaluate the functional objects at a pre-specified number of points on a grid J. The default value is NULL, indicating no change on the number of grid points.
#' @details This function computes and plots functional spherical autocorrelation coefficients at lag \eqn{h}, for \eqn{h = 1:H}.
#' The fSACF at lag \eqn{h} is computed by the average of the inner product of lagged pairs of
#' the series \eqn{X_i} and \eqn{X_{i+h}} that have been centered and scaled:
#' \deqn{
#' \tilde\rho_h=\frac{1}{N}\sum_{i=1}^{N-h} \langle \frac{X_i - \tilde{\mu}}{\|X_i - \tilde{\mu}\|}, \frac{X_{i+h} - \tilde{\mu}}{\|X_{i+h} - \tilde{\mu}\|} \rangle,\ \ \ \ 0 \le h < N,
#' }
#' where \eqn{\tilde{\mu}} is the estimated spatial median of the series.
#' It also computes estimated asymptotic \eqn{(1-\alpha)100 \%} confidence lower and upper bounds, under the SWN assumption.
#'
#' @return Plot of the estimated autocorrelation coefficients for lags \eqn{h} in \eqn{1:H} with the SWN
#' \eqn{(1-\alpha)100 \%} upper and lower confidence bounds for each lag.
#'
#' @references
#' [1] Yeh C.K., Rice G., Dubin J.A. (2023). Functional spherical autocorrelation: A robust estimate of
#' the autocorrelation of a functional time series. Electronic Journal of Statistics, 17, 650–687.
#'
#' @examples
#' \donttest{
#' data(Spanish_elec) # Daily Spanish electricity price profiles
#' fSACF(Spanish_elec)
#' }
#' @export
fSACF <- function (f_data, H = 20, alpha = 0.05, J = NULL) {
  if ((H < 1) | (H%%1 != 0)) {
    stop("The parameter 'H' must be a positive integer.")
  }
  if ((alpha > 1) | (alpha < 0)) {
    stop("The 'alpha' parameter must be a value between 0 and 1.")
  }
  
  data_class <- class(f_data)[[1]]
  if (data_class=="funData" ) {
    if (is.null(J) == TRUE){
      f_data <- t(f_data@X)
    } else {
      f_data <- t(f_data@X)
      J_raw <- NROW(f_data)
      basis <- create.bspline.basis(rangeval = c(0,1), nbasis = 25)
      fd_data <- smooth.basis(1:J_raw/J_raw, y = f_data, fdParobj = basis)$fd
      f_data <- eval.fd(fd_data, 1:J/J)
    }
    
  } else if (data_class=="matrix" ) {
    if (is.null(J) == TRUE){
      f_data <- f_data
    } else {
      J_raw <- NROW(f_data)
      basis <- create.bspline.basis(rangeval = c(0,1), nbasis = 25)
      fd_data <- smooth.basis(1:J_raw/J_raw, y = f_data, fdParobj = basis)$fd
      f_data <- eval.fd(fd_data, 1:J/J)
    }
    
  } else if (data_class=="fd" ) {
    if (is.null(J) == TRUE){
      tempFun <- fd2funData(f_data, argvals = seq(0, 1, length.out = length(f_data[["fdnames"]][["time"]])  ) )
      f_data <- t(tempFun@X) # realization can be different from the original discrete data given the smoothing
    } else {
      tempFun <- fd2funData(f_data, argvals = seq(0, 1, length.out = J)  ) 
      f_data <- t(tempFun@X)
    }
    
  } else {
    stop("The input must be either a matrix or a funData object.")
  }
  
  
  res_raw <- my_new_receipt(t(f_data), H)
  J = NROW(f_data)
  N = NCOL(f_data)
  coefficients = rep(0, H)
  B_iid_bounds_L = rep(0, H)
  B_iid_bounds_U = rep(0, H)
  lags = 1:H

  res_raw = as.matrix(res_raw)

  coefficients =  res_raw[,3]
  B_iid_bounds_L = qt(alpha/2, N - 1) * res_raw[,5]/sqrt(N)
  B_iid_bounds_U = qt(1 - alpha/2, N - 1) * res_raw[,5]/sqrt(N)

  mu<-max(coefficients, B_iid_bounds_U[1])
  ml<-min(coefficients, B_iid_bounds_L[1])
  plot(lags, coefficients, ylim = c(ml-0.05, mu+0.3),
       type = 'h', xlab = "Lag", ylab = "Coefficient",
       main = "fSACF plot")
  lines(B_iid_bounds_L, col = 'red', lty = 'solid')
  lines(B_iid_bounds_U, col = 'red', lty = 'solid')
  lines(rep(0,H), col = 'black', lty = 'solid')

  legend('topleft', legend = c('Estimated Spherical Autocorrelation Coefficients',
                               paste('SWN ', (1 - alpha) * 100,'% Confidence Bound', sep=' ')),
         col = c('black', 'red'),
         lty = c("solid", "solid"), cex=0.9, y.intersp = 1, bty = "n")
}
