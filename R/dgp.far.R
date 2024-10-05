#' FAR(p) Data Generator
#'
#' @description It generates functional data that follows a functional autoregressive process of order \eqn{p}, denoted as FAR(p). The generated data consists of curves evaluated at discrete grid points.
#'
#' @param J The number of grid points for each curve observation.
#' @param N The sample size, representing the number of curves to be generated.
#' @param S The serial dependence factor for the kernel used in the FAR(p) process. Default is 0.5.
#' @param p The order of the autoregressive process. Default is 1.
#' @param kernel The type of kernel function \eqn{\psi} used for the autoregressive process. Can be "Gaussian" or "Wiener". Default is "Gaussian".
#' @param burn_in The number of initial points discarded to eliminate transient effects. Default is 50.
#'
#' @return A \eqn{J \times N} matrix where each column contains a curve evaluated at \eqn{J} grid points, generated from the FAR(p) model.
#'
#' @details The functional autoregressive model of order \eqn{p} is given by:
#' \deqn{X_i(t) -\mu(t) = \sum_{j=1}^{p} \Psi(X_{i-j}-\mu)(t) + \epsilon_i(t),}
#' where \eqn{\Psi(X)(t) = \int \psi(t,s)X(s) dt} is the kernel operator, and \eqn{\epsilon_i(t)} are i.i.d. errors generated from a standard Brownian motion process.
#' The mean function \eqn{\mu} is assumed to be zero in the generating process.
#'
#' @importFrom sde BM
#' @examples
#' \donttest{
#' # Generate discrete evaluations of 200 curves, each observed at 50 grid points.
#' yd_far = dgp.far(J = 50, N = 200, S = 0.7, p = 2, kernel = "Gaussian", burn_in = 50)
#' }
#'
dgp.far <- function (J, N, S = 0.5, p = 1, kernel = "Gaussian", burn_in = 50)
{
  grid <- (1:J)/J

  brown_motion<-function (N, J)
  {
    motion <- matrix(nrow = J, ncol = N)
    for (i in 1:N) {
      motion[, i] <- as.vector(BM(N = J - 1))
    }
    as.array(motion)
  }

  error_mat <- brown_motion(N + burn_in, J)

  if (kernel == "Gaussian") {
    phi <- outer(grid,grid,function(t, s) {exp(-((t^2 + s^2)/2))})
  }
  else if (kernel == "Wiener") {
    phi <- outer(grid,grid, function(t, s) {pmin(t,s)})
  }

  phi_norm_w_out_c <- sqrt(sum(phi^2)/(J^2))

  c <- S/phi_norm_w_out_c

  ## generate FAR(p)
  far_mat <- matrix(0, nrow = J, ncol = N + burn_in)
  far_mat[, 1:p] <- error_mat[, 1:p]
  for (i in (p+1):(N + burn_in)) {

    for (j in 1:p) {
      far_mat[,i] <- far_mat[,i] + c*t(phi)%*%far_mat[,i - j]/J
    }
    far_mat[,i] <- far_mat[,i] + error_mat[,i]

  }

  far_mat[, (burn_in + 1):(burn_in + N)]
}
