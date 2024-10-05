#' Test for Conditional Heteroscedasticity of Functional Time Series
#'
#' @description It tests the null hypothesis that the objective functional curve data is not conditionally heteroscedastic. If a small p-value rejects the null hypothesis,  the curves exhibit conditional heteroscedasticity.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param H A positive integer specifying the maximum lag for which test statistic is computed.
#' @param stat_Method A string specifying the test method to be used in the "ch" test. Options include:
#' \describe{
#'   \item{"norm"}{Uses \eqn{V_{N,H}}.}
#'   \item{"functional"}{Uses \eqn{M_{N,H}}.}
#' }
#' @param pplot A Boolean value. If TRUE, the function will produce a plot of p-values of the test
#' as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.
#'
#' @return
#'  A list that includes the test statistic and the p-value will be returned.
#'
#' @export
#' @importFrom graphics abline plot
#'
#' @import stats
#'
#' @details
#' Given the objective curve data \eqn{X_i(t)}, for \eqn{1\leq i \leq N}, \eqn{t\in[0,1]}, the test aims at distinguishing the hypotheses:\cr \cr
#' \eqn{H_0}: the sequence \eqn{X_i(t)} is IID; \cr
#' \eqn{H_1}: the sequence \eqn{X_i(t)} is conditionally heteroscedastic. \cr \cr
#' Two portmanteau type statistics are applied: \cr \cr
#' 1. the norm-based statistic: \eqn{V_{N,H}=N\sum_{h=1}^H\hat{\gamma}^2_{X^2}(h)}, where \eqn{\hat{\gamma}^2_{X^2}(h)} is the sample autocorrelation of the time series \eqn{||X_1||^2,\dots,||X_N||^2}, and \eqn{H} is a pre-set maximum lag length.\cr \cr
#' 2. the fully functional statistic \eqn{M_{N,H}=N\sum_{h=1}^H||\hat{\gamma}_{X^2,N,h}||^2}, where the autocovariance kernel \eqn{\hat{\gamma}_{X^2,N,h}(t,s)=N^{-1}\sum_{i=1}^{N-h}[X_i^2(t)-\bar{X}^2(t)][X^2_{i+h}(s)-\bar{X}(s)]}, for \eqn{||\cdot ||} is the \eqn{L^2} norm, and \eqn{\bar{X}^2(t)=N^{-1}\sum_{i=1}^N X^2_i(t)}.
#'
#' @examples
#' 
#' \donttest{
#' # generate discrete evaluations of the iid curves under the null hypothesis.
#' yd_ou = dgp.ou(50, 100)
#'
#' # test the conditional heteroscedasticity.
#' fCH_test(yd_ou, H=5, stat_Method="functional")
#' }
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Tests for conditional heteroscedasticity of functional data. Journal of Time Series Analysis. 41(6), 733-758. <doi:10.1111/jtsa.12532>.\cr
fCH_test <- function (f_data, H=10, stat_Method="functional", pplot=FALSE){

  yd = f_data
  K = H

  ## functions to compute the norm statistics
  normx2 <- function(x){ # compute the L2 norm
    N=ncol(x)
    grid_point=nrow(x)
    y1=matrix(NA,1,N)
    x=x-rowMeans(x)
    for(i in 1:N){
      y1[,i]=(1/grid_point)*sum(x[,i]^2)
    }
    y1=y1-mean(y1)
    return(y1)
  }

  gamma1 <- function(x,lag_val){ # compute the covariance of normx2 objective series
    N=ncol(x)
    x=x-mean(x)
    gamma_lag_sum=0
    for(j in 1:(N-lag_val)){
      gamma_lag_sum=gamma_lag_sum+x[j]*x[(j+lag_val)]
    }
    gamma1_h=gamma_lag_sum/N
    return(gamma1_h)
  }

  Port1.stat <- function(dy,Hf){ # the norm statistics
    vals=rep(0,Hf)
    N=ncol(dy)
    y2=normx2(dy)
    sigma2=gamma1(y2,0)
    for(h in 1:Hf){
      gam1hat=gamma1(y2,h)
      vals[h]=N*(gam1hat/sigma2)^2 # the square of autocorrelation function
    }
    stat=sum(vals)
    pv=(1-pchisq(stat,df=Hf))

    return(list(statistic = stat, p_value = pv))
  }

  ## functions to compute the fully functional statistics
  func1 <- function(x){
    N=ncol(x)
    grid_point=nrow(x)
    x=x-rowMeans(x)
    y2=matrix(NA,grid_point,N)
    for(i in 1:N){
      y2[,i]=x[,i]^2
    }
    y2=y2-rowMeans(y2)
    return(y2)
  }

  gamma2 <- function(x, lag_val){ # compute the covariance matrix
    N=ncol(x)
    gam2_sum=0
    for (j in 1:(N-lag_val)){
      gam2_sum=gam2_sum+x[,j]%*%t(x[,(j+lag_val)])
    }
    gam2=gam2_sum/N
    return(gam2)
  }

  Port2.stat <- function(dy,Hf){ # the fully functional statistics
    N=ncol(dy)
    grid_point=nrow(dy)
    vals=rep(0,Hf)
    y2=func1(dy)
    for(h in 1:Hf){
      gam2hat=gamma2(y2,h)
      vals[h]=N*((1/grid_point)^2)*sum(gam2hat^2)}
    stat=sum(vals)
    return(stat)
  }

  # functions to compute the p_value of the fully functional statistics

  etaM <- function(dy){
    N=dim(dy)[2]
    grid_point=dim(dy)[1]
    dy=func1(dy)
    sum1=0
    for(k in 1:N){
      sum1=sum1+(dy[,k])^2
    }
    return(sum1/N)
  }

  mean.W.2 <- function(dy, Hf){ # mean function
    sum1=0
    grid_point=dim(dy)[1]
    etal=etaM(dy)
    sum1=(sum(etal)/grid_point)^2
    sum1=sum1*Hf
    return(sum1)
  }

  etaM_cov <- function(dy){
    N=dim(dy)[2]
    grid_point=dim(dy)[1]
    dy=func1(dy)
    sum1=0
    for(k in 1:N){
      sum1=sum1+(dy[,k])%*%t(dy[,k])
    }
    return(sum1/N)
  }

  etapopNvarMC2 <- function(dy, Hf){ # covariance function
    N=dim(dy)[2]
    grid_point=dim(dy)[1]
    x=etaM_cov(dy)^2
    sum1=sum(x)/(grid_point^2)
    sum1=2*Hf*(sum1^2)
    return(sum1)
  }

  Port.test.2 <- function(dy, Hf){ # Welch–SaNerthwaite approximation
    stat=Port2.stat(dy, Hf)
    res2=mean.W.2(dy, Hf)
    res3=etapopNvarMC2(dy, Hf)
    beta=res3/(2*res2)
    nu=(2*res2^2)/res3
    pv = 1-pchisq(stat/beta, df=nu)
    return(list(statistic = stat, p_value = pv))
  }

  if(is.null(stat_Method) == TRUE) {
    stat_Method="functional"
  }

  switch(stat_Method,
         norm = {
           stat_p=Port1.stat(yd,K)
         },
         functional = {
           stat_p=Port.test.2(yd,K)
         },
         stop("Enter something to switch me!"))

  if(pplot == FALSE) {
    stat_p=stat_p
  }
  else if(pplot==TRUE){
    kseq=1:20
    pmat=c(rep(1,20))
    for (i in 1:20){
      switch(stat_Method,
             norm = {
               pmat[i]=Port1.stat(yd,kseq[i])[2]
             },
             functional = {
               pmat[i]=Port.test.2(yd,kseq[i])[2]
             },
             stop("Enter something to switch me!"))
    }
    x<-1:20

    #par(mar = c(4, 4, 2, 1))
    plot(x,pmat, col='black',ylim=c(0,1.05), pch = 4, cex = 1,
         xlab="H",ylab="p-values", main="p-values of conditional heteroscedasticity test", )
    abline(0.05, 0 , col='red', lty='solid')
    legend('topleft',legend=c("p-values under SWN"),col='black', pch = 4, cex=1.1, bty = 'n')


  }
  return(stat_p)
}



