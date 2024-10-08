% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fSACF_test.R
\name{fSACF_test}
\alias{fSACF_test}
\title{Test based on fSACF}
\usage{
fSACF_test(
  f_data,
  H = 10,
  alpha = 0.05,
  pplot = FALSE,
  suppress_raw_output = FALSE,
  suppress_print_output = FALSE
)
}
\arguments{
\item{f_data}{A \eqn{J \times N} matrix of functional time series data, where \eqn{J} represents the number of discrete points in a grid and \eqn{N} represents the sample size.}

\item{H}{A positive integer specifying the maximum lag for which test statistic is computed.}

\item{alpha}{A numeric value between 0 and 1 indicating the significance level for the test.}

\item{pplot}{A Boolean value. If TRUE, the function will produce a plot of p-values of the test
as a function of maximum lag \eqn{H}, ranging from \eqn{H=1} to \eqn{H=20}, which may increase the computation time.}

\item{suppress_raw_output}{A Boolean value. If TRUE, the function will not return the list
containing the p-value, quantile, and statistic.}

\item{suppress_print_output}{A Boolean value. If TRUE, the function will not print any
output to the console.}
}
\value{
If suppress_raw_output = FALSE, a list containing the test statistic, the \eqn{(1-\alpha)} quantile of the
limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
containing a short description of the test, the p-value, and additional information about the test if
suppress_print_output = FALSE.
}
\description{
This function performs a hypothesis test using a test statistic computed from functional spherical autocorrelation coefficients of a FTS.
}
\details{
The test statistic is the sum of the squared \eqn{L^2}-norm of the sample spherical autocorrelation coefficients:
\deqn{
S_{N,H} = N \sum_{h=1}^H \|\tilde{\rho}_{h}\|^2,
}
where \eqn{\tilde\rho_h=\frac{1}{N}\sum_{i=1}^{N-h} \langle \frac{X_i - \tilde{\mu}}{\|X_i - \tilde{\mu}\|}, \frac{X_{i+h} - \tilde{\mu}}{\|X_{i+h} - \tilde{\mu}\|} \rangle},
   and \eqn{\tilde{\mu}} is the estimated spatial median of the series.
This test assesses the cumulative significance of lagged spherical autocorrelation coefficients, up to a
user-selected maximum lag \eqn{H}. A higher value of \eqn{S_{N,H}} suggests a potential
departure of the observed series from white noise process. The approximated null
distribution of this statistic is developed under the strong white noise assumptions.
}
\examples{
\donttest{
data(Spanish_elec)
fACF_test(Spanish_elec, H = 20)
}
}
\references{
[1] Yeh CK, Rice G, Dubin JA (2023). “Functional spherical autocorrelation: A robust estimate of
the autocorrelation of a functional time series.” Electronic Journal of Statistics, 17, 650–687.
}
