% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OCIDR.R
\name{OCIDR}
\alias{OCIDR}
\title{Convert Original Price Data to OCIDRs}
\usage{
OCIDR(f_data)
}
\arguments{
\item{f_data}{A \eqn{J \times N} matrix of functional time series data, where \eqn{J} is the number of discrete points in a grid and \eqn{N} is the sample size.}
}
\value{
A matrix of OCIDRs with dimensions \eqn{J \times (N-1)}, where \eqn{J} is the number of discrete grid points and \eqn{N-1} is the adjusted sample size.
}
\description{
This function converts original price data into over-night cumulative intraday return curves (OCIDRs).
}
\examples{
\donttest{
data(sp500)
OCIDR(sp500)
}
}
