#' Test based on fSACF
#'
#' This function performs a hypothesis test using a test statistic computed from functional spherical autocorrelation coefficients of a FTS.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} represents the number of discrete points in a grid and \eqn{N} represents the sample size.
#' @param H A positive integer specifying the maximum lag for which test statistic is computed.
#' @param alpha A numeric value between 0 and 1 indicating the significance level for the test.
#' @param pplot A Boolean value. If TRUE, the function will produce a plot of p-values of the test
#' as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.
#' @param suppress_raw_output A Boolean value. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output A Boolean value. If TRUE, the function will not print any
#' output to the console.
#'
#' @details The test statistic is the sum of the squared \eqn{L^2}-norm of the sample spherical autocorrelation coefficients:
#' \deqn{
#' S_{N,H} = N \sum_{h=1}^H \|\tilde{\rho}_{h}\|^2,
#' }
#' where \eqn{\tilde\rho_h=\frac{1}{N}\sum_{i=1}^{N-h} \langle \frac{X_i - \tilde{\mu}}{\|X_i - \tilde{\mu}\|}, \frac{X_{i+h} - \tilde{\mu}}{\|X_{i+h} - \tilde{\mu}\|} \rangle},
#'    and \eqn{\tilde{\mu}} is the estimated spatial median of the series.
#' This test assesses the cumulative significance of lagged spherical autocorrelation coefficients, up to a
#' user-selected maximum lag \eqn{H}. A higher value of \eqn{S_{N,H}} suggests a potential
#' departure of the observed series from white noise process. The approximated null
#' distribution of this statistic is developed under the strong white noise assumptions.
#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the \eqn{(1-\alpha)} quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE.
#' @references
#' [1] Yeh CK, Rice G, Dubin JA (2023). “Functional spherical autocorrelation: A robust estimate of
#' the autocorrelation of a functional time series.” Electronic Journal of Statistics, 17, 650–687.
#' @examples
#' \donttest{
#' data(Spanish_elec)
#' fACF_test(Spanish_elec, H = 20)
#' }
#' @importFrom rainbow fts
#' @importFrom ftsa ftsm
#'
#' @export
fSACF_test <- function(f_data, H = 10, alpha = 0.05, pplot= FALSE,
                              suppress_raw_output=FALSE, suppress_print_output=FALSE) {
  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. Atleast one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }

  if ((H < 1) | (H %% 1 != 0)) {
    stop("The 'components lag must be a positive integer.")
  }

  res_raw <- my_new_receipt(t(f_data), max(H,20))

  J = NROW(f_data)
  N = NCOL(f_data)

  res_raw = as.matrix(res_raw)

  res_final<-calc_BP_test(rho = res_raw[,2], H = max(H,20),
               alpha = alpha, N = N, cp = res_raw[1,5])

  if (pplot == TRUE) {

    pvalues_iid = res_final[1:20,4]
    lags = 1:20

    #par(mar = c(4, 4, 2, 1))
    plot(lags, pvalues_iid, ylim=c(0,1.05), pch = 4, cex = 1, col='black', xlab='H',
         ylab='P-value', main = 'p-values of spherical test')
    abline(0.05, 0 , col='red', lty='solid')

    legend('topleft',
           legend=c('p-values under SWN'),
           col='black', pch = 4, cex=1.1, bty = 'n')
  }

  results <- list(statistic = res_final[H,2], quantile = res_final[H,3], p_value = res_final[H,4])


  if (suppress_print_output == FALSE) {
    title_print <- sprintf("Spherical Test\n\n")
    test_type <- 'the series is a strong white noise\n'
    null_print <- sprintf("null hypothesis: %s", test_type)
    p_val_print <- sprintf("p-value = %f\n", results$p_value)
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    lag_print <- sprintf("maximum lag = %d\n\n\n", H)
    message(c(title_print, null_print, p_val_print, samp_print, lag_print))
  }
  if (suppress_raw_output == FALSE) {
    return(results)
  }


}

