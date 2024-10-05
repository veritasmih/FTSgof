#' Test based on fACF
#'
#' This function performs a hypothesis test using a test statistic computed from functional autocovariance kernels of a FTS.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param H A positive integer specifying the maximum lag for which test statistic is computed.
#' @param iid A Boolean value. If given TRUE, the hypothesis test will use the strong-white
#' noise (SWN) assumption instead of the weak white noise (WWN) assumption.
#' @param M A positive integer specifying the number of Monte Carlo simulations used to approximate the null distribution under the WWN assumption.
#' If \eqn{M = NULL, M = \text{floor}((\max(150 - N, 0) + \max(100 - J, 0) + (J / \sqrt{2})))},
#' ensuring that the number of Monte Carlo simulations is adequate based on the dataset size.
#' @param alpha A numeric value between 0 and 1 indicating the significance level for the test.
#' @param pplot A Boolean value. If TRUE, the function will produce a plot of p-values of the test
#' as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.
#' @param suppress_raw_output A Boolean value. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output A Boolean value. If TRUE, the function will not print any
#' output to the console.
#'
#'
#' @details The test statistic is the sum of the squared \eqn{L^2}-norm of the sample autocovariance kernels:
#' \deqn{
#' KRS_{N,H} = N \sum_{h=1}^H \|\hat{\gamma}_{N,h}\|^2,
#' }
#' where
#'    \eqn{
#'    \hat{\gamma}_{N,h}(t,s)=N^{-1}\sum_{i=1}^{N-h} (X_i(t)-\bar{X}_N(t))(X_{i+h}(s)-\bar{X}_N(s))},
#'    \eqn{\bar{X}_N(t) = N^{-1} \sum_{i=1}^N X_i(t)}.
#' This test assesses the cumulative significance of lagged autocovariance kernels, up to a
#' user-selected maximum lag \eqn{H}. A higher value of \eqn{KRS_{N,H}} suggests a potential
#' departure of the observed series from white noise process. The approximated null
#' distribution of this statistic is developed under both the strong and weak white noise assumptions.
#'
#'
#' @return If suppress_raw_output = FALSE, a list that includes the test statistic, the \eqn{(1-\alpha)} quantile of the
#' limiting distribution, and the p-value from the specified hypothesis test. Additionally, if suppress_print_output = FALSE,
#' a summary is printed with a brief explanation of the test, the p-value, and relevant details about the test procedure.
#'
#'
#' @references
#' [1] Kokoszka P., Rice G., Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @examples
#' \donttest{
#' data(sp500) # S&P500 index
#' fACF_test(OCIDR(sp500), H = 10, pplot=TRUE)
#' }
#' @import stats
#' @export
fACF_test <- function(f_data, H = 10, iid = FALSE, M = NULL, pplot = FALSE,
                           alpha = 0.05, suppress_raw_output = FALSE,
                           suppress_print_output = FALSE) {
  if ((H < 1) | (H %% 1 != 0)) {
    stop("The parameter 'H' must be a positive integer and greater than 1.")
  }
  if ((alpha > 1) | (alpha < 0)) {
    stop("The 'alpha' parameter must be a value between 0 and 1.")
  }

  K <- H

  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. At least one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }

  if (pplot == TRUE) {
    
    lags = 1:20
    
    if (iid == FALSE) {
      
      pvalues_wn = rep(0, 20)
      for (h in lags){
        pvalues_wn[h] <- V_WS_quantile(f_data, h, alpha=alpha, M=M)$p_value
      }
      
      #par(mar = c(4, 4, 2, 1))
      plot(lags, pvalues_wn, ylim=c(0,1.05), pch = 1, cex = 1, col='blue', xlab='H',
           ylab='p-value', main = 'p-values of autocovariance test')
      abline(0.05, 0 , col='red', lty='solid')
      
      legend('topleft',
             legend=c('p-values under WWN'),
             col=c('blue'), pch = 1, cex=1.1, bty = 'n')
      
    }else {
      
      pvalues_iid = rep(0, 20)
      for (h in lags){
        pvalues_iid[h] <- V_WS_quantile_iid(f_data, h, alpha=alpha)$p_value
      }
      
      #par(mar = c(4, 4, 2, 1))
      plot(lags, pvalues_iid, ylim=c(0,1.05), pch = 4, cex = 1, col='black', xlab='H',
           ylab='p-value', main = 'p-values of autocovariance test')
      abline(0.05, 0 , col='red', lty='solid')
  
      legend('topleft',
             legend=c('p-values under SWN'),
             col=c('black'), pch = 4, cex=1.1, bty = 'n')
  }
  }
  
  if (iid == FALSE) {
      results <- V_WS_quantile(f_data, K, alpha=alpha, M=M)
      if (suppress_print_output == FALSE) {
        title_print <- sprintf("Autocovariance Test (weak white noise assumption)\n\n")
        test_type <- 'the series is a weak white noise\n'
        null_print <- sprintf("null hypothesis: %s", test_type)
        p_val_print <- sprintf("p-value = %f\n", results$p_value)
        samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
        lag_print <- sprintf("maximum lag = %d\n", K)
        mc_print <- sprintf("number of monte-carlo simulations = %d\n\n\n", M)
        message(c(title_print, null_print, p_val_print, samp_print,
                  lag_print, mc_print))
      }
      if (suppress_raw_output == FALSE) {
        return(results)
      }

    }else {
      results <- V_WS_quantile_iid(f_data, K, alpha=alpha)
      if (suppress_print_output == FALSE) {
        title_print <- sprintf("Autocovariance Test (iid assumption)\n\n")
        test_type <- 'the series is a strong white noise\n'
        null_print <- sprintf("null hypothesis: %s", test_type)
        p_val_print <- sprintf("p-value = %f\n", results$p_value)
        samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
        lag_print <- sprintf("maximum lag = %d\n\n\n", K)
        message(c(title_print, null_print, p_val_print, samp_print,
                  lag_print))
      }

      if (suppress_raw_output == FALSE) {
        return(results)
      }

   }
}
