#' Ornstein–Uhlenbeck Process Generator
#'
#'
#' @description It generates iid functional curve data following the Ornstein–Uhlenbeck process.
#'
#' @param J The number of grid points in each curve observation.
#' @param N The sample size.
#'
#' @return A (grid points) x (number of observations) matrix for iid sequences, where the finite realization of curves are stored in columns.
#' @export
#'
#' @importFrom MASS mvrnorm
#'
#' @details
#' The Ornstein–Uhlenbeck process is given by:
#' \eqn{x_i(t)=e^{-t/2}W_i(e^t)}, \eqn{t \in [0,1]},\cr
#' where \eqn{W_i(t)} is a standard Brownian Motion.
#'
#' @examples
#' 
#' \donttest{
#' # Generate discrete evaluations of 100 iid curves
#' # that each curve is realized on 50 grid points.
#' yd_ou = dgp.ou(J = 50, N = 100)
#' }
#'
dgp.ou <- function(J, N){
  times=1:J/J
  # control the covariance structure
  comat=matrix(NA,J,J)
  for (i in 1:J){
    comat[i,]=exp(1)^(-times[i]/2-times/2)*pmin(exp(1)^(times[i]),exp(1)^(times))
  }
  fiid=mvrnorm(n = N, mu = c(rep(0,J)), Sigma = comat, empirical = FALSE)

  return(t(fiid))
}
