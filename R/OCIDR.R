#' Convert Original Price Data to OCIDRs
#'
#' @description This function converts original price data into over-night cumulative intraday return curves (OCIDRs).
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @return A matrix of OCIDRs with dimensions \eqn{J \times (N-1)}, where \eqn{J} is the number of discrete grid points and \eqn{N-1} is the adjusted sample size.
#' @examples
#' \donttest{
#' data(sp500)
#' OCIDR(sp500)
#' }
#' @export
OCIDR<-function(f_data){
  temp=as.matrix(f_data)
  J=nrow(f_data)
  N=ncol(f_data)
  ocidr_re=matrix(0,J,N-1)
  ocidr_re<-log(temp[,2:N]) - matrix(log(temp[J,1:N-1]), J, N-1, byrow=TRUE)
  return(ocidr_re)
}















