#' Goodness-of-fit test for FAR(1)
#'
#' @description It fits a FAR(1) model and then assesses the cumulative significance of lagged
#' autocovariance operators from the model residuals, up to a user-selected maximum lag \eqn{H}.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param H A positive integer specifying the maximum lag for which test statistics are computed.
#' @param alpha A numeric value between 0 and 1 specifying the significance level.
#' @param M A positive integer specifying the number of Monte Carlo simulations used to approximate the null distribution.
#' If \eqn{M = NULL, M = \text{floor}((\max(150 - N, 0) + \max(100 - J, 0) + (J / \sqrt{2})))},
#' ensuring that the number of Monte Carlo simulations is adequate based on the dataset size.
#' @param pplot A Boolean value. If TRUE, the function will produce a plot of p-values of the test
#' as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.
#' @param suppress_raw_output A Boolean value, FALSE by default. If TRUE, the function will not return the list
#' containing the p-value, quantile, and statistic.
#' @param suppress_print_output A Boolean value, FALSE by default. If TRUE, the function will not print any
#' output to the console.
#' @param residual A data frame. If TRUE, the function will provide the residuals obtained from fitting the FAR(1) model.

#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the \eqn{(1-\alpha)} quantile of the
#' limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#' containing a short description of the test, the p-value, and additional information about the test if
#' suppress_print_output = FALSE.
#'
#' @references
#' [1] Kim, M., Kokoszka, P., & Rice, G. (2023). White noise testing for functional time series. Statistic Surveys, 17, 119-168.
#'
#' @examples
#' \donttest{
#' yd_far <- dgp.far(J=50, N=100, S=0.7, p=2, kernel = "Gaussian", burn_in = 50)
#' gof_far(yd_far, H=5, pplot=TRUE)
#' }
#' @import stats
#' @export
gof_far<-function(f_data, H=10, M=NULL, alpha=0.05, pplot= FALSE, residual = FALSE, suppress_raw_output=FALSE, suppress_print_output=FALSE)
{
  if (!requireNamespace('fda')) {
    stop("Please install the 'fda' package to perform the GOf for FAR(1) test.")
  }

  if (suppress_raw_output == TRUE & suppress_print_output == TRUE) {
    stop("Current choice of parameters will produce no output. At least one of the parameters
         'suppress_raw_output' or 'suppress_print_output' must be FALSE.")
  }

  J <- dim(f_data)[1]
  N <- dim(f_data)[2]

  basis <- fda::create.bspline.basis(rangeval=c(0,1), nbasis=20, norder=4)
  fd.data <- fda::smooth.basis(0:(J-1)/(J-1), f_data, basis)$fd
  pca <- fda::pca.fd(fd.data[1:(N-1)], nharm=20, fda::fdPar(fd.data), centerfns=TRUE)
  kN <- which.max(cumsum(pca$values)/sum(pca$values)>0.90)

  score <- as.matrix(pca$scores[,1:kN])
  eigenval <- pca$values[1:kN]
  xi <- score%*%diag(1/eigenval,kN)%*%t(score)/N

  X <- fda::eval.fd(fd.data,0:(J-1)/(J-1))
  X <- X-rowMeans(X)
  eg <- X[,2:N]-X[,2:N]%*%xi
  eg <- cbind(X[,1],eg)

  f <- array(NA,c(N-1,J,max(H,20)))
  for(h in 1:max(H,20))
  {
    f[,,h] <- xi[,h:(N-1)]%*%t(eg[,1:(N-h)])
  }

  if (pplot == TRUE) {

    pvalues_iid = rep(0, 20)
    lags = 1:20

    for (h in lags){
      pvalues_iid[h] <- V_WS_quantile_far(eg, f, lag=h, alpha, M)$p_value
    }

    #par(mar = c(4, 4, 2, 1))
    plot(lags, pvalues_iid, ylim= c(0,1.05), pch = 4, cex = 1, col='black', xlab='H',
         ylab='p-value', main = 'p-values of goodness-of-fit for FAR(1)')
    abline(0.05, 0 , col='red', lty='solid')

    legend('topleft',
           legend=c('p-values under the FAR(1) null'),
           col=c('black'), pch = c(4), cex=1.1, bty = 'n')
  }

  results <-V_WS_quantile_far(eg, f, lag=H, alpha, M)

  if (suppress_print_output == FALSE) {
    title_print <- sprintf("Goodness-of-fit test for FAR(1)\n\n")
    null_print <- sprintf("null hypothesis: FAR(1) model is adequate for the series.\n")
    p_val_print <- sprintf("p-value = %f\n", results$p_value)
    samp_print <- sprintf("sample size = %d\n", NCOL(f_data))
    lag_print <- sprintf("lag = %d\n\n\n", H)
    message(c(title_print, null_print, p_val_print, samp_print,
              lag_print))
  }
  if (suppress_raw_output == FALSE) {

    if (residual == TRUE) {
      return(list(statistic = results$statistic, quantile = results$quantile, p_value = results$p_val, res= eg))
    }
    return(results)


  }


}
