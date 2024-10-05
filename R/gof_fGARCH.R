#' Goodness-of-fit Test for Functional ARCH/GARCH Model
#'
#' @description It tests the goodness-of-fit of functional ARCH/GARCH models by accounting for the effect of functional GARCH parameter estimation.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.
#' @param M A positive integer specifying the number of basis functions.
#' @param model A string to indicate which model will be estimated: "arch" - FARCH(1); "garch" - FGARCH(1,1).
#' @param H A positive integer specifying the maximum lag for which test statistics are computed.
#' @param pplot A Boolean value. If TRUE, the function will produce a plot of p-values of the test
#' as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.
#' @param max_eval The maximum number of evaluations of the optimization function, used in the "arch" and "garch" tests.
#'
#'
#' @return p-value.
#'
#' @importFrom graphics abline plot
#' @importFrom nloptr cobylac
#' @importFrom fda create.bspline.basis smooth.basis eval.fd pca.fd fdPar
#'
#' @details
#' It tests the goodness-of-fit of the fARCH(1) or fGARCH(1,1) models.
#' It fits the model to the input data and applies the test \eqn{M_{N,H}} in \code{\link{fport_wn}} to the model residuals.
#' The asymptotic distribution is adjusted to account for the estimation effect,
#' because the model residual depends on the joint asymptotics of the innovation process and
#' the estimated parameters. We assume that the kernel parameters are consistently estimated
#' by the Least Squares method proposed in Aue et al. (2017).
#' Then, the asymptotic distribution of the statistic \eqn{M_{N,H}} is given in Theorem 3.1
#' in Rice et al. (2020).
#'
#'
#' @examples
#' \donttest{
#' # generate discrete evaluations of the FGARCH process.
#' set.seed(42)
#' yd = dgp.fgarch(J=50, N=200, type = "garch")$garch_mat
#'
#' # test the adequacy of the FARCH(1) model.
#' gof_fGARCH(yd, M=2, model = "arch", H=10, pplot=TRUE)
#' }
#' @references
#' [1] Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#'
#' [2] Rice, G., Wirjanto, T., Zhao, Y. (2020). Tests for conditional heteroscedasticity of functional data. Journal of Time Series Analysis. 41(6), 733-758. <doi:10.1111/jtsa.12532>.\cr
#'
#' @export
gof_fGARCH <- function (f_data, M, model, H=10, pplot=NULL, max_eval=10000){

  yd = f_data
  K=H
  sample_size = ncol(yd)
  grid_point = nrow(yd)

  yd=yd-rowMeans(yd)
  times=rep(0,grid_point)
  for(i in 1:grid_point){
    times[i]=i/grid_point
  }
  y2m=yd*yd
  basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=32,norder=4)
  y_sq_fdsmooth=smooth.basis(argvals=times,y=y2m,fdParobj=basis_obj)
  y_sq_fd=y_sq_fdsmooth$fd
  pca_obj=pca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd), centerfns=TRUE)  #centerfns=FALSE.
  eigen_functs=pca_obj$harmonics
  eigen_v=pca_obj$values
  tve=eigen_v/sum(eigen_v)
  ortho_basis_matrix=matrix(0,nrow=grid_point,ncol=10)
  for(j in 1:10){
    indi_basis=eval.fd(times,eigen_functs[j])
    indi_basis[indi_basis<0]<-0
    ortho_basis_matrix[,j]=indi_basis
  }
  basis = as.matrix(ortho_basis_matrix[,1:M])

  int_approx=function(x){
    temp_n=NROW(x)
    return((1/temp_n)*sum(x))}

  y_vec=matrix(0,nrow=M,ncol=sample_size)
  for(i in 1:sample_size){
    for(j in 1:M){y_vec[j,i]=int_approx(y2m[,i]*basis[,j])}}

  switch(model,
         arch = {
           get_theta=function(theta,M){
             #first M entries are d
             d=theta[1:M]
             #Next p*M^2 entries are A, A_inds are indices of A matrix params
             num_mat_params=(M^2-M*(M-1)/2)
             A_inds=M+(1:(num_mat_params))
             theta_A=theta[A_inds]

             curr_A_vals=theta_A[1:num_mat_params]
             A=matrix(0,ncol=M,nrow=M)
             diag(A)=curr_A_vals[1:M]
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
             A=t(A)
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
             return(list(ds=d,As=A))
           }

           function_to_minimize2=function(x,data=y_vec){
             sample_size=ncol(data)
             M=nrow(data)
             pam = get_theta(x,M)
             d=pam$d
             A=pam$As
             sigma_sq_proj_coefs=matrix(0,M,sample_size)
             sigma_sq_proj_coefs[,1]=d

             for(i in 2:(sample_size)){
               first_part=d
               second_part=A%*%data[,i-1]
               sigma_sq_proj_coefs[,i]=first_part+second_part
             }

             s_theta=0
             for (i in 2:sample_size) {s_theta=s_theta+t(data[,i]-sigma_sq_proj_coefs[,i])%*%(data[,i]-sigma_sq_proj_coefs[,i])}
             r=sum(s_theta)
             return(r)
           }

           #nonlinear inequality constraints on parameters
           eval_g0<-function(x,Mn=M){
             conpam=get_theta(x,Mn)
             A = conpam$As
             upb = 1-(10^-20)
             h <- -(sum(A)-upb)
             return(h)
           }
           num_params=M^2-M*(M-1)/2+M
           stav=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))
           ress=nloptr::cobyla(x0 = stav,fn = function_to_minimize2, lower=c(rep(10^-20,num_params)),
                               upper=c(rep(1,num_params)),
                               hin = eval_g0, nl.info = FALSE,
                               control = list(maxeval=max_eval))

           para=as.numeric(ress$par)
           pam_hat = get_theta(para,M)
           d_hat=pam_hat$d
           A_hat=pam_hat$As

           d_M=0
           for (i in 1:M){
             d_M=d_M+d_hat[i]*basis[,i]}

           alpha_M=matrix(0,grid_point,grid_point)
           for(i in 1:grid_point){
             for(j in 1:grid_point){
               temp_alpha_m=0
               for(t in 1:M){
                 for(s in 1:M){temp_alpha_m=temp_alpha_m+A_hat[t,s]*basis[i,t]*basis[j,s]                  }
               }
               alpha_M[i,j]=temp_alpha_m
             }}

           sigma_fit=matrix(1,grid_point,(sample_size+1))
           sigma_fit[,1]=d_M
           for(j in 2:(sample_size+1)){
             #first fill in sigma2 column:
             for(i in 1:grid_point){
               fit_alpha_op = alpha_M[i,] * ((yd[,(j-1)])^2)
               sigma_fit[i,j] = d_M[i] + int_approx(fit_alpha_op)
             }      }

           error_fit=yd/sqrt(abs(sigma_fit[,1:sample_size]))

           sigma_2_proj_coefs=matrix(0,M,sample_size)
           sigma_2_proj_coefs[,1]=d_hat
           for(i in 2:(sample_size)){
             first_part=d_hat
             second_part=A_hat%*%y_vec[,i-1]
             sigma_2_proj_coefs[,i]=first_part+second_part
           }
         },
         garch = {
           get_theta=function(theta,M){
             #first M entries are d
             d=theta[1:M]

             #Next p*M^2 entries are A, A_inds are indices of A matrix params
             num_mat_params=(M^2-M*(M-1)/2)
             A_inds=M+(1:(num_mat_params))
             B_inds=sfsmisc::last(A_inds)+(1:(num_mat_params))

             theta_A=theta[A_inds]
             theta_B=theta[B_inds]

             curr_A_vals=theta_A[1:num_mat_params]
             A=matrix(0,ncol=M,nrow=M)
             diag(A)=curr_A_vals[1:M]
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
             A=t(A)
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]

             #see above comment
             curr_B_vals=theta_B[1:num_mat_params]
             B=matrix(0,ncol=M,nrow=M)
             diag(B)=curr_B_vals[1:M]
             B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]
             B=t(B)
             B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]

             return(list(ds=d,As=A,Bs=B))
           }

           function_to_minimize2=function(x,data=y_vec){
             sample_size=ncol(data)
             M=nrow(data)
             pam = get_theta(x,M)
             d=pam$d
             A=pam$As
             B=pam$Bs
             sigma_sq_proj_coefs=matrix(0,M,sample_size)
             sigma_sq_proj_coefs[,1]=d

             for(i in 2:(sample_size)){
               first_part=d
               second_part=0
               for (l in 1:(i-1)){
                 first_part=first_part+B^(l-1)%*%d
                 second_part=second_part+B^(l-1)%*%(A%*%data[,i-l])}
               sigma_sq_proj_coefs[,i]=first_part+second_part
             }
             s_theta=0
             for (i in 2:sample_size) {s_theta=s_theta+t(data[,i]-sigma_sq_proj_coefs[,i])%*%(data[,i]-sigma_sq_proj_coefs[,i])}
             r=sum(s_theta)
             return(r)
           }

           #nonlinear inequality constraints on parameters
           eval_g0<-function(x,Mn=M){
             h <- numeric(2)
             upb = 1-(10^-20)
             conpam=get_theta(x,M)
             A=conpam$As
             B=conpam$Bs
             h[1] <- -(sum(A)-upb)
             h[2] <- -(sum(B)-upb)
             return(h)
           }
           num_params=2*(M^2-M*(M-1)/2)+M
           stav=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))
           ress=nloptr::cobyla(x0 = stav,fn = function_to_minimize2, lower=c(rep(10^-20,num_params)),
                               upper=c(rep(1,num_params)),
                               hin = eval_g0, nl.info = FALSE,
                               control = list(maxeval=max_eval))
           
           para=as.numeric(ress$par)
           pam_hat = get_theta(para,M)
           d_hat=pam_hat$d
           A_hat=pam_hat$As
           B_hat=pam_hat$Bs

           d_M=0
           for (i in 1:M){
             d_M=d_M+d_hat[i]*basis[,i]}

           alpha_M=matrix(0,grid_point,grid_point)
           beta_M=matrix(0,grid_point,grid_point)
           for(i in 1:grid_point){
             for(j in 1:grid_point){
               temp_alpha_m=0
               temp_beta_m=0
               for(t in 1:M){
                 for(s in 1:M){
                   temp_alpha_m=temp_alpha_m+A_hat[t,s]*basis[i,t]*basis[j,s]
                   temp_beta_m=temp_beta_m+B_hat[t,s]*basis[i,t]*basis[j,s]
                 }
               }
               alpha_M[i,j]=temp_alpha_m
               beta_M[i,j]=temp_beta_m
             }}

           sigma_fit=matrix(1,grid_point,(sample_size+1))
           sigma_fit[,1]=d_M
           for(j in 2:(sample_size+1)){
             #first fill in sigma2 column:
             for(i in 1:grid_point){
               fit_alpha_op = alpha_M[i,] * ((yd[,(j-1)])^2)
               fit_beta_op = beta_M[i,] * ((sigma_fit[,j-1]))
               sigma_fit[i,j] = d_M[i] + int_approx(fit_alpha_op) + int_approx(fit_beta_op)
             }     }

           error_fit=yd/sqrt(abs(sigma_fit[,1:sample_size]))

           sigma_2_proj_coefs=matrix(0,M,sample_size)
           sigma_2_proj_coefs[,1]=d_hat
           for(i in 2:(sample_size)){
             first_part=d_hat
             second_part=A_hat%*%y_vec[,i-1]+B_hat%*%sigma_2_proj_coefs[,i-1]
             sigma_2_proj_coefs[,i]=first_part+second_part
           }
         },
         stop("Enter something to switch me!"))

  #### goodness-of-fit test

  gof_pvalue<-function(error_fit,sigma_fit,ortho_basis_matrix,y_2_proj_coefs,sigma_2_proj_coefs,para,K){

    ## functions to calculate the fully functional statistics.
    int_approx_tensor<-function(x){# x is a 4-dimensional tensor
      dt=length(dim(x))
      temp_n=nrow(x)
      return(sum(x)/(temp_n^dt))}

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

    ## functions to calculate covariance operator

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

    H0_f<-function(y_2_proj_coefs1){
      # H0 is a M x (M+2M^2) matrix
      M=nrow(y_2_proj_coefs1)
      n=ncol(y_2_proj_coefs1)

      ded=diag(M)
      dey=matrix(0, M, M^2)
      h0_sum=0
      #des=matrix(0, M, M^2)########  for garch
      for (i in 2:n){
        if (M>1){
          for (j in 1:M){
            dey[j,(((j-1)*M+1):(j*M))]=y_2_proj_coefs1[,(i-1)]}}
        else {dey[1,1]=y_2_proj_coefs1[,(i-1)]}
        h0_c=cbind(ded,dey)
        h0_sum=h0_sum+h0_c}

      return(h0_sum/(n-1))
      #return(rbind(ded,dey,des))  #  for garch
    }

    Q0_f<-function(y_2_proj_coefs1){
      # Q0 is a (M+2M^2) x (M+2M^2) matrix
      M=nrow(y_2_proj_coefs1)
      n=ncol(y_2_proj_coefs1)

      dey=matrix(0, M, M^2)
      #des=matrix(0, M, M^2)########  for garch
      q0_sum=0
      for (i in 2:n){
        if (M>1){
          for (j in 1:M){
            dey[j,(((j-1)*M+1):(j*M))]=y_2_proj_coefs1[,(i-1)]}}
        else {dey[1,1]=y_2_proj_coefs1[,(i-1)]}
        ded=diag(M)
        q0_c=cbind(ded,dey)
        q0_e=t(q0_c)%*%q0_c
        q0_sum=q0_sum+q0_e
      }
      return(q0_sum/(n-1))
    }

    J0_f<-function(y_2_proj_coefs1,sigma_2_proj_coefs){
      # J0 is a M x M matrix
      sumdys=0
      for (i in 1: ncol(y_2_proj_coefs1)){
        dys=y_2_proj_coefs1[,i]-sigma_2_proj_coefs[,i]
        sumdys=sumdys+dys%*%t(dys)
      }
      return(sumdys/ncol(y_2_proj_coefs1))
    }

    covariance_op_new<-function(error_fit,sigma_fit,ortho_basis_matrix,y_2_proj_coefs1,sigma_2_proj_coefs,para,h,g){
      ortho_basis_matrix=ortho_basis_matrix/(sum(ortho_basis_matrix)/nrow(error_fit))
      TT = dim(error_fit)[2]
      p = dim(error_fit)[1]

      M = dim(y_2_proj_coefs1)[1]
      sum1 = 0
      sum1_m = 0
      error_fit = func1(error_fit)
      for (k in 1:TT){
        sum1 = sum1 + (error_fit[,k]) %*% t(error_fit[,k])
        sum1_m = sum1_m + (error_fit[,k])^2
      }

      A = (sum1/TT)%o%(sum1/TT)
      cov_meanA = (sum1_m/TT)

      H0_est = H0_f(y_2_proj_coefs1)
      Q0_est = Q0_f(y_2_proj_coefs1)
      J0_est = J0_f(y_2_proj_coefs1,sigma_2_proj_coefs)
      innermatrix = solve(Q0_est)%*%t(H0_est)%*%J0_est%*%H0_est%*%solve(Q0_est)

      # kronecker %x%
      #sum2 = array(0, c(p, p, p, p))
      rep_permu<-function(M){
        permu=matrix(NA,M^2,2)
        permu[,1]=c(rep(c(1:M),each=M))
        permu[,2]=c(rep(c(1:M),M))
        return(permu)
      }
      permu = rep_permu(M)

      d_hat = para[1:M]
      A_hat = matrix(para[(M+1):(M+M^2)],nrow=M,ncol=M)

      G_g = array(0, c((M+M^2),p,p))

      for (k in 2:(TT-g)){
        p2_g=matrix(0,(M+M^2),p)
        parti_g=matrix(0,(M+M^2),1)
        for (m in 1:(M+M^2)){
          if ( m <= M ){
            vecm1=c(replicate(M,0))
            vecm1[m]=1
            partd = sum(vecm1)
            parti_g[m,1]=partd
          }
          else{
            vecm2=matrix(0,M,M)
            vecm2[permu[m-M,1],permu[m-M,2]]=1
            parta = sum(vecm2%*%y_2_proj_coefs1[,k+g-1])
            parti_g[m,1]= parta
          }
        }
        for (j in 1:(M+M^2)){
          if (j<=M){p2_g[j,]=ortho_basis_matrix[,j]*parti_g[j,1]}#(j-(ceiling(j/M)-1)*M)
          else{
            p2_g[j,]=(rowSums(ortho_basis_matrix[,(permu[j-M,1])]%o%ortho_basis_matrix[,(permu[j-M,2])])/p)*parti_g[j,1]}
        }
        G_g = G_g + ((t(replicate((M+M^2),(1/sigma_fit[,k+g])))*p2_g)%o%(error_fit[,k]))
      }
      G_g = G_g/(TT-g-1)

      G_h = array(0, c((M+M^2),p,p))

      for (k in 2:(TT-h)){
        p2_h=matrix(0,(M+M^2),p)
        parti_h=matrix(0,(M+M^2),1)
        for (m in 1:(M+M^2)){
          if ( m <= M ){
            vecm1=c(replicate(M,0))
            vecm1[m]=1
            partd = sum(as.matrix(vecm1))
            parti_h[m,1]=partd
          }
          else{
            vecm2=matrix(0,M,M)
            vecm2[permu[m-M,1],permu[m-M,2]]=1
            parta = sum(vecm2%*%y_2_proj_coefs1[,k+h-1])
            parti_h[m,1]= parta
          }

        }
        for (j in 1:(M+M^2)){
          if (j<=M){p2_h[j,]=ortho_basis_matrix[,j]*parti_h[j,1]}#(j-(ceiling(j/M)-1)*M)
          else{
            p2_h[j,]=(rowSums(ortho_basis_matrix[,(permu[j-M,1])]%o%ortho_basis_matrix[,(permu[j-M,2])])/p)*parti_h[j,1]}
        }
        G_h=G_h + ((t(replicate((M+M^2),(1/sigma_fit[,k+h])))*p2_h)%o%(error_fit[,k]))
      }
      G_h = G_h/(TT-h-1)

      arr_3_prod<-function(mat1,Lmat,mat2){
        p=dim(mat1)[2]
        L=dim(Lmat)[1]
        mat1t=aperm(mat1)
        x=array(0,c(p,p,p,p))
        for (i in 1:p){
          for (j in 1:p){
            x[i,j,,] = x[i,j,,] + mat1t[i,,]%*%Lmat%*%mat2[,,j]
          }
        }
        return(x)
      }
      B = arr_3_prod(G_h,innermatrix,G_g)
      rm(parti_h)
      rm(parti_g)
      rm(p2_h)
      rm(p2_g)

      sum3 = array(0, c(p, p, p, p))

      arr_prod<-function(mat1,Lmat){
        p=dim(mat1)[2]
        L=dim(Lmat)[1]
        mat1t=aperm(mat1)
        x=matrix(0,p,p)
        for (i in 1:p){
          x[i,] = x[i,] + mat1t[i,,]%*%Lmat
        }
        return(x)
      }

      for (k in 2:(TT-g)){
        smtx=matrix(0,M,M^2)
        for (m in 1:M^2){
          vecm2=matrix(0,M,M)
          vecm2[permu[m,1],permu[m,2]]=1
          parta = sum(vecm2%*%y_2_proj_coefs1[,k-1])
          smtx[permu[m,1],m] = parta}#t(y_2_proj_coefs1[,(k-1)])}
        partii=cbind(diag(M),smtx)
        innercov2 = as.matrix((solve(Q0_est) %*% t(partii) %*% (y_2_proj_coefs1[,k]-sigma_2_proj_coefs[,k])),(M+M^2),1)
        EG_g=arr_prod(G_g,innercov2)
        sum3 = sum3 + (error_fit[,k]) %o% (error_fit[,k+h]) %o% EG_g
      }
      C=sum3/(TT-g-1)
      D=array(0, c(p, p, p, p))
      for (i in 1:p){
        for (j in 1:p){
          D[i,j,,]=C[,,i,j]
        }
      }

      cov_est=B+C+D
      cov_mean=matrix(NA,p,p)
      for (i in 1:p){
        for (j in 1:p){
          cov_mean[i,j]=cov_est[i,j,i,j]
        }
      }
      return(list(A,cov_est,cov_mean,cov_meanA))
    }

    ## compute the fully functional statistics and p_values

    stat = Port2.stat(error_fit,K)
    mu_cov=0
    sigma2_cov=0

    len=20 # parameters in Monte Carlo integration
    len2=15
    grid_point=dim(error_fit)[1]

    rrefind <- floor(grid_point * sort(runif(len, 0, 1)))
    rrefind[which(rrefind == 0)] = 1
    r_samp <- c(0,rrefind) # c(0, rrefind) we define v_0 = 0 for computation of the product in rTrap
    xd <- diff(r_samp / grid_point)
    gmat_A <- xd %o% xd
    gmat <- xd %o% xd %o% xd %o% xd

    rrefind2 <- floor(grid_point * sort(runif(len2, 0, 1)))
    rrefind2[which(rrefind2 == 0)] = 1
    r_samp2 <- c(0,rrefind2) # we define v_0 = 0 for computation of the product in rTrap
    xd2 <- diff(r_samp2 / grid_point)
    gmat2_A <- xd2 %o% xd2
    gmat2 <- xd2 %o% xd2 %o% xd2 %o% xd2

    mu_cov=0
    sigma2_cov=0

    if (K>1){ # for K is greater than 1
      for (hh in 1:K){
        mean_A=array(0, c(len))
        mean_mat=array(0, c(len, len))
        cov_mat_A=array(0, c(len, len, len, len))
        cov_mat=array(0, c(len, len, len, len))

        covop_trace=covariance_op_new(error_fit[rrefind,],sigma_fit[rrefind,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind,],length(rrefind),1),y_2_proj_coefs,sigma_2_proj_coefs,para,hh,hh)
        mean_A[2:len] = covop_trace[[4]][-len]
        mean_mat[2:len,2:len] = covop_trace[[3]][-len,-len]
        mu_cov=mu_cov+sum((1/2) * (covop_trace[[4]]+mean_A)*xd)^2
        mu_cov=mu_cov+sum((1/2) * (covop_trace[[3]]+mean_mat)*gmat_A)^2

        cov_mat_A[2:len,2:len,2:len,2:len] = covop_trace[[1]][-len,-len,-len,-len]
        cov_mat[2:len,2:len,2:len,2:len] = covop_trace[[2]][-len,-len,-len,-len]
        sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_trace[[1]]^2+cov_mat_A^2)*gmat)
        sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_trace[[2]]^2+cov_mat^2)*gmat)
      }

      for (h1 in 1:K){
        for (h2 in (h1+1):K){
          if((abs(h2-h1)) < 3){
            cov_mat=array(0, c(len, len, len, len))

            covop_est=covariance_op_new(error_fit[rrefind,],sigma_fit[rrefind,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind,],length(rrefind),1),y_2_proj_coefs,sigma_2_proj_coefs,para,h1,h2)
            cov_mat[2:len,2:len,2:len,2:len] = covop_est[[2]][-len,-len,-len,-len]
            sigma2_cov=sigma2_cov+2*2*sum((1/2)*(covop_est[[2]]^2+cov_mat^2)*gmat)
          }

          if((abs(h2-h1)) >= 3){
            cov_mat=array(0, c(len2, len2, len2, len2))
            covop_est=covariance_op_new(error_fit[rrefind2,],sigma_fit[rrefind2,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind2,],length(rrefind2),1),y_2_proj_coefs,sigma_2_proj_coefs,para,h1,h2)
            cov_mat[2:len2,2:len2,2:len2,2:len2] = covop_est[[2]][-len2,-len2,-len2,-len2]
            sigma2_cov=sigma2_cov+2*2*sum((1/2)*(covop_est[[2]]^2+cov_mat^2)*gmat2)
          }
        }
      }
    }
    else{ # for K=1
      mean_A=array(0, c(len))
      mean_mat=array(0, c(len, len))
      cov_mat_A=array(0, c(len, len, len, len))
      cov_mat=array(0, c(len, len, len, len))

      covop_est=covariance_op_new(error_fit[rrefind,],sigma_fit[rrefind,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind,],length(rrefind),1),y_2_proj_coefs,sigma_2_proj_coefs,para,1,1)

      mean_A[2:len] = covop_est[[4]][-len]
      mean_mat[2:len,2:len] = covop_est[[3]][-len,-len]

      mu_cov=mu_cov+sum((1/2) * (covop_est[[4]]+mean_A)*xd)^2
      mu_cov=mu_cov+sum((1/2) * (covop_est[[3]]+mean_mat)*gmat_A)^2

      cov_mat_A[2:len,2:len,2:len,2:len] = covop_est[[1]][-len,-len,-len,-len]
      cov_mat[2:len,2:len,2:len,2:len] = covop_est[[2]][-len,-len,-len,-len]

      sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_est[[1]]^2+cov_mat_A^2)*gmat)
      sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_est[[2]]^2+cov_mat^2)*gmat)
    }
    sigma2_cov=2*sigma2_cov

    beta = sigma2_cov/(2 * mu_cov)
    nu = (2 * mu_cov^2)/sigma2_cov

    return(1 - pchisq(stat/beta, df = nu))
  }

  if( is.null(K) == TRUE) {
    K=20}
  stat_pvalue=gof_pvalue(error_fit,sigma_fit,basis,y_vec,sigma_2_proj_coefs,para,K)

  if( is.null(pplot) == TRUE) {
    stat_pvalue=stat_pvalue
  }
  else if(pplot==1){
    kseq=1:20
    pmat=c(rep(1,20))
    for (i in 1:20){
      pmat[i]=gof_pvalue(error_fit,sigma_fit,basis,y_vec,sigma_2_proj_coefs,para,kseq[i])
    }
    #par(mar = c(4, 4, 2, 1))
    x<-1:20


    if(model=="arch"){
    plot(x,pmat, col="black", ylim=c(0,1.05), pch = 4, cex = .8,
         xlab="H",ylab="p-values",main="p-values of goodness-of-fit for fGARCH(1,0)")
    abline(h=0.05, lty = 'solid', col="red")
    legend('topleft',legend=c("p-values under the fGARCH(1,0) null"), col=c("black"), pch = 4, cex=1.1, bty = 'n')

    }else if(model=="garch"){
      plot(x,pmat, col="black", ylim=c(0,1.05), pch = 4, cex = .8,
           xlab="H",ylab="p-values",main="p-values of goodness-of-fit for fGARCH(1,1)")
      abline(h=0.05, lty = 'solid', col="red")
      legend('topleft',legend=c("p-values under the fGARCH(1,1) null"), col=c("black"), pch = 4, cex=1.1, bty = 'n')

    }



  }




  return(stat_pvalue)
}



