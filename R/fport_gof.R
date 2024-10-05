#' Goodness-of-fit Tests for Functional Times Series
#'
#' @description It computes three different goodness-of-fit tests for functional time series. All goodness-of-fit tests in this package are accessible through this function.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param test A string specifying the goodness-of-fit test. Currently available tests are referred
#' to by their string handles: "far", "arch" and "garch". Please see the Details section of the documentation.
#' @param H A positive integer specifying the maximum lag for which test statistics are computed.
#' @param M A positive integer specifying the number of Monte Carlo simulations used to approximate the null distribution in the "far" test, and the number of basis functions used in the "arch" and "garch" tests. 
#' If \eqn{M = NULL}, then for the "far" test, \eqn{M} is calculated as \eqn{\text{floor}((\max(150 - N, 0) + \max(100 - J, 0) + (J / \sqrt{2})))}, ensuring the number of Monte Carlo simulations is adequate based on the dataset size.
#' For the "arch" or "garch" test, \eqn{M} is set to 1.
#' @param pplot A Boolean value. If TRUE, the function will produce a plot of p-values of the test
#' as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.
#' @param residual A data frame. If TRUE, the function will provide the residuals obtained from fitting the FAR(1) model.
#'
#' @details
#' This function conducts goodness-of-fit tests for functional time series data,
#' providing several testing options:
#'
#' 1. When test = "far", it tests the goodness-of-fit of the FAR(1) model.
#' The function fits the FAR(1) model to the input data and then applies the test statistic \eqn{KRS_{N,H}}, as described in \code{\link{fport_wn}}, to the residuals.
#' The null distribution of the test statistic accounts for the dependence structure present in the residuals.
#' The optional parameters for this test are 'fdata', 'test', 'H', 'M', 'pplot' and 'residual'.
#'
#' 2. When test = "arch" or "garch", it tests the goodness-of-fit of the fARCH(1) or fGARCH(1,1) models.
#' It fits the model to the input data and applies the test \eqn{M_{N,H}} in \code{\link{fport_wn}} to the model residuals.
#' The asymptotic distribution is adjusted to account for the estimation effect,
#' because the model residual depends on the joint asymptotics of the innovation process and
#' the estimated parameters. We assume that the kernel parameters are consistently estimated
#' by the Least Squares method proposed in Aue et al. (2017).
#' Then, the asymptotic distribution of the statistic \eqn{M_{N,H}} is given in Theorem 3.1
#' in Rice et al. (2020).
#' The optional parameters for this test are 'fdata', 'test', 'H', 'M', and 'pplot'.
#'
#'
#'
#' @return
#' A summary is printed with a brief explanation of the test and the p-value.
#'
#'
#' @references
#' [1] Kim, M., Kokoszka, P., & Rice, G. (2023). White noise testing for functional time series. Statistic Surveys, 17, 119-168.
#'
#' [2] Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#'
#' [3] Rice, G., Wirjanto, T., Zhao, Y. (2020). Tests for conditional heteroscedasticity of functional data. Journal of Time Series Analysis. 41(6), 733-758. <doi:10.1111/jtsa.12532>.\cr
#'
#'
#' @examples
#' \donttest{
#' data(Spanish_elec)
#' fport_gof(Spanish_elec, test = "far", H = 20, pplot=TRUE)
#'
#' data(sp500)
#' fport_gof(OCIDR(sp500), test = "arch", M = 1, H = 5)
#' fport_gof(OCIDR(sp500), test = "garch", M = 1, H = 10)
#' }
#' @export
#'
#'
#' @import stats
#' @importFrom graphics abline plot
#' @importFrom nloptr cobylac
#' @importFrom fda create.bspline.basis smooth.basis eval.fd pca.fd fdPar
#'
fport_gof <- function(f_data, test = "far", H=10, M=NULL, 
                      pplot=FALSE, residual = FALSE) {

  tests = c("far", "arch", "garch")
  if (!(test %in% tests)) {
    stop("Please see the documentation for available tests.")
  }
  if (!is.null(H)) {
    if (!all.equal(H, as.integer(H)) | H <= 0) {
      stop("Invalid arguments, lag must be a positive integer")
    }
  }
  if (!is.matrix(f_data)) {
    stop("Invalid arguments, functional data f_data must be passed in matrix form.")
  }
  if (!is.null(M)) {
    if (!all.equal(M, as.integer(M)) | M < 0) {
      stop("Invalid arguments, M must be a positive integer or NULL.")
    }
  }
  if (is.null(M) && (test %in% c("arch", "garch"))) {
    M=1
  }
  if (test == "far") {

    final<-gof_far(f_data, H, M, alpha=0.05, pplot, residual, suppress_print_output=TRUE)
    if (residual == TRUE){
      #title_print <- sprintf("Goodness-of-fit test for FAR(1)\n\n")
      #null_print <- sprintf("Null hypothesis: FAR(1) model is adequate for the series.\n")
      #samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
      #lag_print <- sprintf("maximum lag H = %d\n", H)
      #p_val_print <- sprintf("p-value = %f\n", final$p_value)
      #message(c(title_print, null_print, null_print, samp_print, p_val_print))

      return(list(resid = final$res))

      #return(list(statistic = final$statistic, p_value = as.numeric(final$p_value), resid = final$res))
    }

    title_print <- sprintf("Goodness-of-fit test for FAR(1)\n\n")
    null_print <- sprintf("Null hypothesis: FAR(1) model is adequate for the series.\n")
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    lag_print <- sprintf("maximum lag H = %d\n", H)
    p_val_print <- sprintf("p-value = %f\n", final$p_value)
    message(c(title_print, null_print, samp_print, lag_print, p_val_print))

    #return(list(statistic = final$statistic, p_value = as.numeric(final$p_value)))

  } else if (test == "arch") {
    final<-gof_fGARCH(f_data, M, model="arch", H, pplot, max_eval=10000)

    title_print <- sprintf("Goodness-of-fit test for fARCH(1)\n\n")
    null_print <- sprintf("Null hypothesis: fARCH(1) model is adequate for the series.\n")
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    lag_print <- sprintf("maximum lag H = %d\n", H)
    p_val_print <- sprintf("p-value = %f\n", final)
    message(c(title_print, null_print, samp_print, lag_print, p_val_print))

  } else if (test == "garch") {
    final<-gof_fGARCH(f_data, M, model="garch", H, pplot, max_eval=10000)

    title_print <- sprintf("Goodness-of-fit test for fGARCH(1,1)\n\n")
    null_print <- sprintf("Null hypothesis: fGARCH(1,1) model is adequate for the series.\n")
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    lag_print <- sprintf("maximum lag H = %d\n", H)
    p_val_print <- sprintf("p-value = %f\n", final)
    message(c(title_print, null_print, samp_print, lag_print, p_val_print))

  }

}
