#' 3D Rainbow Plot for Functional Time Series
#'
#' @description This function creates a 3D rainbow plot of a functional time series (FTS) using the \code{rgl} package.
#' The data is color-coded using a rainbow color scheme, with earlier curves in red and later curves in violet.
#' A 360-degree interactive view is available through an embedded widget.
#'
#' @param f_data A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in the grid and \eqn{N} is the sample size.
#' @return A 3D rainbow plot of the functional time series data.
#'
#' @examples
#' \donttest{
#' data(Spanish_elec) # Daily Spanish electricity price profiles
#' rainbow3D(Spanish_elec)
#'
#' data(sp500) # S&P500 index data
#' rainbow3D(OCIDR(sp500))
#' }
#' @importFrom rgl open3d plot3d lines3d axes3d title3d rglwidget
#' @export
#'
rainbow3D <- function(f_data) {
  .onLoad <- function(libname, pkgname) {
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("The 'rgl' package is required but not installed. Please install it using `install.packages(\"rgl\")`. Additionally, 'rgl' requires XQuartz to be installed on your system. You can download and install XQuartz from https://www.xquartz.org/.")
    }
  }


  f_data <- t(f_data)
  N <- dim(f_data)[1]
  J <- dim(f_data)[2]

  x <- 1:N
  y <- as.matrix(f_data)


  open3d()

  plot3d(x, 1:ncol(y), y, type = "n", xlab = "", ylab = "", zlab = "", box = FALSE)

  cols <- rainbow(1000)[(x - min(x)) / (max(x) - min(x)) * 1000 * 7/9 + 1]
  for (i in seq_along(x)) {
    lines3d(x[i], 1:ncol(y), y[i, ], col = cols[i], lwd = 2.5)
  }

  axes3d()
  title3d(xlab = "1:N", ylab = "1:J", zlab = "")

  rglwidget()
  
}

