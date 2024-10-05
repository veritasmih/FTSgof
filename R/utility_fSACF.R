# Store all utility functions used in the 'fSACF' function

gsj <- function(x, y){
  n <- nrow(x)
  m <- ncol(x)
  if (n == 1){
    x <-  t(x)
    y <- t(y)
    n <- nrow(x)
    m <- ncol(x)
  }

  if (n < 3)
    return(NULL)
  ind <- apply(x, 2, order)
  x <- apply(x, 2, sort)

  # Also sort the y accordingly using the indicator above
  for (j in 1:m) {
    y[ ,j] <- y[ind[ ,j], j]
  }

  a <- (x[3:n, ] - x[2:(n-1), ]) / (x[3:n, ] - x[1:(n-2), ])
  b <- (x[2:(n-1), ] - x[1:(n-2), ])/(x[3:n, ]-x[1:(n-2), ])
  c <- a^2 + b^2 + 1
  e <- a * y[1:(n-2), ] + b * y[3:n, ] - y[2:(n-1),]

  estvar <- colMeans(e^2/c)
  return(estvar)
}

calc_norm_matrix <- function(obs){
  Nt <- ncol(obs)
  N <- nrow(obs)
  sapply(1:N, function(i){
    sqrt(sum(obs[i, ]^2))/sqrt(Nt-1)
  })
}

SpMed <- function(tt, x, dtyp = 's'){
  tt <- as.matrix(tt, ncol = length(tt), nrow = 1)
  eps <- 2^(-52)
  # Input check
  if (nargs() < 2)
    stop("Not enough input variables")

  n <- nrow(x)
  m <- ncol(x)

  if (nrow(tt) > 1)
    tt <- t(tt)

  if (ncol(tt) != m)
    stop('Dimensions of T and X not compatible')

  if (!(dtyp != 'n' | dtyp != 's'))
    stop('Wrong input value for DTYP')

  # Inner-product matrix
  A <- 0.5 * (x[, 1:(m-1)] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 1:(m-1)]) +
                x[, 2:m] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 2:m]))

  if( tolower(dtyp) == "n"){
    se2 <- gsj(kronecker(matrix(1, 1, n), tt), t(x))
    A1 <- A - diag(se2) * (tt[m] - tt[1]);
    if (min(eigen(A1)>-1e-6))
      A <- A1
    else{
      message('Corrected inner-product matrix is not nonnegative definite \n Using uncorrected matrix instead ')
    }
  }

  ## Iterative minimization from sample mean
  w <- matrix(1, nrow = n, ncol = 1)/n # a very naive way
  norms <- sqrt(diag(A) + as.vector(t(w) %*% A %*% w) - 2 * A %*% w) # pay attention to the R syntax that it does not like to do arithmetic on a number to vector
  f <- sum(norms)
  err <- 1
  iter <- 0

  while( err > 1E-5 && iter < 50){
    iter <- iter + 1
    f0 <- f
    if (min(norms < eps)){
      i0 = find(norms < eps)
      w <- matrix(0, nrow = n, ncol = 1)
      w[i0] <- 1/length(i0)
    } else{
      w <- 1/norms
      w <- w/sum(w)
    }

    norms <- sqrt(diag(A) + as.vector(t(w) %*% A %*% w) - 2 * A %*% w) # same here
    f <- sum(norms)
    err <- abs(f/f0 - 1)
  }

  med <- t(w) %*% x
  return(list(med = med, w = w, norms = norms))
}

calc_var <- function(xbr){
  Nt <- ncol(xbr)
  N <- nrow(xbr)
  res <- matrix(NA, Nt, Nt)
  for(s in 1:Nt){
    for(t in 1:Nt){
      val <- 0
      for(i in 1:N){
        val <- val + xbr[i, s] * xbr[i, t]
      }
      res[s, t] <- (val/N)^2
    }
  }
  sum(res)/(Nt-1)^2
}


calc_var_raw <- function(x_raw){
  N <- nrow(x_raw)
  Nt <- ncol(x_raw)

  Spatial_med <- as.numeric(SpMed(tt = seq(0, 1, length.out = Nt), x_raw)$med)
  x_cen <- sapply(1:Nt, function(k){
    x_raw[, k] - Spatial_med[k]
  })

  my_norm <- calc_norm_matrix(x_cen)
  xbr <- x_cen/my_norm

  res <- matrix(NA, Nt, Nt)
  for(s in 1:Nt){
    for(t in 1:Nt){
      val <- 0
      for(i in 1:N){
        val <- val + xbr[i, s] * xbr[i, t]
      }
      res[s, t] <- (val/N)^2
    }
  }
  sum(res)/(Nt-1)^2
}



my_new_receipt <- function(obs, H){
  N <- nrow(obs)
  Nt <- ncol(obs)

  norm_raw <- calc_norm_matrix(obs)

  SpMedian <- as.numeric(SpMed(seq(0, 1, length.out = Nt), obs)$med)

  obs_center <- sapply(1:N, function(i){
    obs[i, ] - SpMedian
  }) |>  t()

  norm_center <- calc_norm_matrix(obs_center)

  rho_raw <- sapply(1:H, FUN = function(h){
    res <- sapply(1:(N-h), function(i){
      sum(obs[i, ]/norm_raw[i] * obs[i+h, ]/norm_raw[i+h])/(Nt-1)
    })
    sum(res)/N
  })

  rho_cen <- sapply(1:H, FUN = function(h){
    res <- sapply(1:(N-h), function(i){
      sum(obs_center[i, ]/norm_center[i] * obs_center[i+h, ]/norm_center[i+h])/(Nt-1)
    })
    sum(res)/N
  })
  my_var_raw <- calc_var(obs/norm_raw)
  my_var_cen <- calc_var_raw(obs)
  cbind(1:H, rho_raw,  rho_cen, std_0_raw = rep(sqrt(my_var_raw),H), std_0_cen = rep(sqrt(my_var_cen),H))

}



calc_BP_test <- function(rho, H, alpha, N, cp) {
  results <- lapply(1:H, function(h) {
    val <- N * sum(rho[1:h]^2) / cp^2
    qval <- qchisq(alpha, df = h, lower.tail = FALSE)
    pval <- pchisq(val, df = h, lower.tail = FALSE)
    data.frame(lag = h, val_bp = val, qval_bp = qval, pval_bp = pval)
  })

  do.call(rbind, results)
}

