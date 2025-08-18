#' Nonparametric Limits
#' @name nonp_limits
#'
#' @description
#' Calculate shortest nonparametric reference limits, at given confidence level.
#'
#' @usage
#' nonp_limits(X, RR=TRUE, perc=0.95, cover=0.9)
#'
#' @param X  the data set whose reference limit is sought.
#' @param RR  *optional* (default of TRUE) - the way percentiles are defined:
#'    * FALSE -- use the Weibull limit i/(n+1),
#'    * TRUE  -- use the Hazen limit (i-0.5)/n.
#' @param perc	*optional* (default of 0.95) - the two-sided probability.
#' @param cover	*optional* (default of 0.9) - the confidence level of the CI for the reference limit.
#'
#' @details
#' The reference limits are estimated as lower and upper percentiles of the sample.
#' There are many ways of defining sample percentiles of which the code has
#' two options – the Hazen limit (preferred) and the Weibull limit
#' (which has historical precedents).
#'
#' The confidence limits on these percentiles come from standard binomial methodology.
#'    * Each interval is a pair of order statistics from the sample.
#'    * The codes find the order statistics with the required 90% coverage
#'    but use order statistics with indices as close together as possible.
#'    * An implication of this is that the 90% confidence interval on each
#'    reference limit has a total tail area less than 10% but does not
#'    necessarily have less than 5% in each tail.
#'
#' If the sample size is too small (<91 for the default settings), confidence
#' limits can not be computed and will be returned as NA.
#' If it is even smaller (<21 for the default settings) the estimated
#' reference limits will also be returned as NA.
#'
#' @returns A list containing the following components:
#'
#'   \item{lower}{the lower limit and confidence interval}
#'   \item{upper}{the upper limit and confidence interval}
#'   \item{a}{the index of the order statistic defining the lower limit of the CI}
#'   \item{b}{the index of the order statistic defining the upper limit of the CI}
#'   \item{coverage}{the actual confidence level the interval achieves}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @example /Example/nonp_limits_example.R
#'
#' @references Horn PS, Peske AJ (2005). Reference intervals: a user’s guide. Washington (DC): AACC Press.
#'
#' @importFrom stats pbinom
#' @importFrom stats dbinom
#'
#' @export

nonp_limits = function(X, RR=TRUE, perc=0.95, cover=0.9) {

  X         <- X[!is.na(X)]
  N         <- length(X)
  oneperc   <- (1 + perc) / 2
  targorder <- 0.5 + N * oneperc        # Hazen definition
  if (!RR) targorder <- oneperc * (N+1) # Weibull definition
  sorter    <- sort(X)
  upper     <- NA
  lower     <- NA
  if (targorder < N & targorder > 1) {
    q         <- floor(targorder)
    q1        <- q+1
    leftover  <- targorder - q
    notleft   <- 1 - leftover
    upper     <- sorter[q]     * notleft + leftover * sorter[q1]
    lower     <- sorter[N+1-q] * notleft + leftover * sorter[N+1-q1]
  }

  # Search for CI
  a      <- 1
  b      <- N-1
  sofar  <- dbinom(0, N, oneperc) + dbinom(N, N, oneperc)
  down   <- dbinom(a , N, oneperc)
  up     <- dbinom(b, N, oneperc)

  for (looper in 1:N) {
    if (down < up) {
      if (sofar + down > 1 - cover) break
      sofar <- sofar + down
      a     <- a+1
      down <- dbinom(a, N, oneperc)
    } else {
      if (sofar + up > 1 - cover) break
      sofar <- sofar + up
      b <- b - 1
      up <- dbinom(b, N, oneperc)
    }
  }
  a <- a
  b <- b + 1
  below  <- pbinom(a-1, N, oneperc)
  atend  <- pbinom(b-1, N, oneperc)
  calc   <- atend - below
  uplow  <- sorter[a]
  uphi   <- sorter[b]
  lowlow <- sorter[N+1-b]
  lowhi  <- sorter[N+1-a]
  lower=c(lower,lowlow,lowhi)
  upper=c(upper,uplow,uphi)
  if (sofar > 1-cover) {
    lower[2:3] <- rep(NA,2)
    upper[2:3] <- rep(NA,2)
    a     <- NA
    b     <- NA
  }
  return(list(lower=lower, upper=upper, a=a, b=b, coverage=1-sofar))
}
