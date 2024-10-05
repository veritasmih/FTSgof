#' Functional ARCH/GARCH Process Generator
#'
#' @description It generates functional curve data following the functional ARCH(1) or GARCH(1,1) process.
#'
#' @param J The number of grid point in each curve observation.
#' @param N The sample size.
#' @param type A string to switch data generating process between "arch" - functional ARCH and "garch" - functional GARCH.
#' @param alpha_par The ARCH kernel coefficient function in the conditional volatility equation. If it is missing, "\eqn{12 t (1-t) s (1-s)}" is used to generate FGARCH, and "\eqn{16 t (1-t) s (1-s)}" is used to generate FARCH, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#' @param beta_par The GARCH kernel coefficient function in the conditional volatility equation. If it is missing, "\eqn{12 t (1-t) s (1-s)}" is used to generate FGARCH, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#'
#' @return List of generated processes:
#' @return garch_mat: FARCH/GARCH sequences, where the finite realization of curves are stored in columns;
#' @return sigma_mat: Conditional volatility sequences,  where the finite realization of curves are stored in columns.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#'
#' @details
#' If \eqn{X_i(t)} follows an FARCH(1) process,\cr
#' \eqn{X_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) X^2_{i-1}(s)ds}.\cr
#'
#' If \eqn{X_i(t)} follows an FGARCH(1,1) process,\cr
#' \eqn{X_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) X^2_{i-1}(s)ds+\int \beta(t,s) \sigma^2_{i-1}(s)ds},\cr
#' where the innovation \eqn{\varepsilon_i(t)} follows an Ornstein–Uhlenbeck process \code{\link{dgp.ou}}, and the constant coefficient \eqn{\omega(t)=0.1t(1-t)}.
#'
#' @seealso \code{\link{dgp.ou}}
#' @examples
#' \donttest{
#' # Generate discrete evaluations of 100 fGARCH curves that
#' # each curve is realized on 50 grid points.
#' yd = dgp.fgarch(J = 50, N = 100, type = "garch")
#' yd_garch = yd$garch_mat
#' }
#' @references
#' [1] Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the ARCH model. Econometric Theory. 29(2), 267-288. <doi:10.1017/S0266466612000345>.\cr
#' [2] Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#'
dgp.fgarch <- function(J, N, type, alpha_par=NULL, beta_par=NULL){
  times = 1:J/J
  int_approx <- function(x){
    temp_n = nrow(x)
    return((1/(temp_n)) * sum(x))
  }

  switch(type,
         arch = {
           no_proc = 1000 + N
           sim_sigma2_matrix = sim_garch_matrix = matrix(NA, J, no_proc)
           error_matrix = dgp.ou(J = J, N = no_proc)

           if(is.null(alpha_par) == TRUE) {
             alpha_par = function(t,s){
               return(16 * t * (1-t) * s * (1-s))
             }
           }
           delta_par = function(t){
             return(0.1* t * (1-t))
           }
           #fill the initial values
           sim_sigma2_matrix[,1] = delta_par(times)
           sim_garch_matrix[,1] = delta_par(times) * error_matrix[,1]

           #fill in the rest of sigma2 and functional garch process matrices:
           for(j in 2:no_proc){
             #first fill in sigma2 column:
             for(i in 1:J){
               vector_op = alpha_par(times[i], times) * ((sim_garch_matrix[,(j-1)])^2)
               sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(as.matrix(vector_op))
             }
             #fill in garch process values for column:
             sim_garch_matrix[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
           }
         },
         garch = {
           no_proc = 1000 + N
           sim_sigma2_matrix = sim_garch_matrix = matrix(NA, J, no_proc)
           error_matrix = dgp.ou(J = J, N = no_proc)

           if( is.null(alpha_par) == TRUE) {
             alpha_par = function(t,s){
               return(12 * t * (1-t) * s * (1-s))
             }
           }
           if( is.null(beta_par) == TRUE) {
             beta_par = function(t,s){
               return(12 * t * (1-t) * s * (1-s))
             }
           }
           delta_par = function(t){
             return(0.5* t * (1-t))
           }
           #fill the initial values
           sim_sigma2_matrix[,1] = delta_par(times)
           sim_garch_matrix[,1] = delta_par(times) * error_matrix[,1]

           #fill in the rest of sigma2 and functional garch process matrices:
           for(j in 2:no_proc){
             #first fill in sigma2 column:
             for(i in 1:J){
               vector_for_alpha_op = alpha_par(times[i], times) * ((sim_garch_matrix[,(j-1)])^2)
               vector_for_beta_op = beta_par(times[i], times) * ((sim_sigma2_matrix[,j-1]))
               sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(as.matrix(vector_for_alpha_op)) + int_approx(as.matrix(vector_for_beta_op))
             }
             #fill in garch process values for column:
             sim_garch_matrix[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
           }
         },
         stop("Enter something to switch me!"))
  return(list(garch_mat = sim_garch_matrix[,1001:(1000 + N)],
              sigma_mat = sim_sigma2_matrix[,1001:(1000 + N)]))
}
