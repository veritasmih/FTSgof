#' Exploratory Data Analysis for Functional Time Series.
#'
#' @description This function sequentially displays the fACF plot, the fSACF and the rainbow plot of a functional time series (FTS) for comprehensive exploratory data analysis.
#'
#' @param f_data A \eqn{J \times N} matrix of FTS data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param H A positive integer representing the maximum lag for computing the coefficients and confidence bounds. This value determines the range of lags included in the fACF and fSACF plots.
#' @param alpha A numeric value between 0 and 1 indicating the significance level for the confidence bounds in the fACF and fSACF plots.
#' @param wwn_bound A Boolean value allowing the user to turn on the WWN bound in the fACF plot. FALSE by default. Speeds down computation when TRUE.
#' @param M A positive integer value. The number of Monte-Carlo simulations used to compute the confidence bounds under the WWN assumption. 
#' If \eqn{M = NULL, M = \text{floor}((\max(150 - N, 0) + \max(100 - J, 0) + (J / \sqrt{2})))},
#' ensuring that the number of Monte Carlo simulations is adequate based on the dataset size.
#' @details This function sequentially displays the rainbow plot, the fACF plot, and the fSACF of an FTS for comprehensive exploratory data analysis.
#' See the help page of \code{\link{rainbow3D}}, \code{\link{fACF}}, \code{\link{fSACF}}, for more details.
#'
#' @return A 3D rainbow plot, a fACF plot for lags \eqn{h \in 1:H} with the WWN
#' \eqn{(1-\alpha)100 \%} upper confidence bound and the constant strong white noise (SWN)
#' \eqn{(1-\alpha)100 \%} upper confidence bound, and a fSACF plot for lags \eqn{h \in 1:H} with the SWN
#' \eqn{(1-\alpha)100 \%} upper and lower confidence bounds.
#'
#' @references
#' [1] Kokoszka P., Rice G., Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' [2] Mestre G., Portela J., Rice G., Roque A. M. S., Alonso E. (2021). Functional time series model identification
#'  and diagnosis by means of auto-and partial autocorrelation analysis. Computational Statistics & Data Analysis, 155, 107108.
#'
#' [3] Yeh CK, Rice G, Dubin JA (2023). “Functional spherical autocorrelation: A robust estimate of
#' the autocorrelation of a functional time series.” Electronic Journal of Statistics, 17, 650–687.
#'
#' @examples
#' \donttest{
#' data(Spanish_elec) # Daily Spanish electricity price profiles
#' fport_eda(Spanish_elec)
#' }
#'
#' @export
#' @importFrom rgl open3d plot3d lines3d axes3d title3d rglwidget
#' @import stats
#'
fport_eda <- function(f_data, H = 20, alpha = 0.05, wwn_bound = FALSE, M = NULL) {
  
  if ((H < 1) | (H %% 1 != 0)) {
    stop("The parameter 'H' must be a positive integer.")
  }
  if ((alpha > 1) | (alpha < 0)) {
    stop("The 'alpha' parameter must be a value between 0 and 1.")
  }

  fACF(f_data, H, alpha, wwn_bound, M)
  readline(prompt = "Hit <Return> to see next plot: ")

  fSACF(f_data, H, alpha)
  readline(prompt = "Hit <Return> to see next plot: ")
  
  rainbow3D(f_data)
  return(rglwidget())
}


