#' Daily Eurodollar Futures Curves
#'
#' @description This dataset contains daily Eurodollar futures curves from February 9, 1998 to June 5, 1998 (\eqn{N=82}).
#' A Eurodollar futures contract represents an obligation to deliver 1,000,000 USD to a bank
#' outside the United States at a specified time.
#' The Eurodollar futures curves consist of daily settlement prices for these contracts,
#' available at monthly delivery dates for the first six months and quarterly delivery dates
#' up to 10 years into the future. These curves are preprocessed using cubic splines,
#' following Kargin and Onatski (2008), to transform the raw data into smooth curves on
#' a grid of 114 equally spaced points (\eqn{J=114}).
#'
#' @format A matrix with columns representing
#' the daily settlement prices as observed functions.
#'
#' @references
#' Kargin V, Onatski A (2008). Curve forecasting by functional autoregression. Journal of
#' Multivariate Analysis, 99, 2508–2526.
#'
#' @usage
#' data(EF)
#'
'EF'


