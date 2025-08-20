#' Parametric Limits, based on Normal
#' @name para_limits
#'
#' @description
#' Parametric confidence limit for a normal or normal-in-the-middle sample.
#'
#' @usage
#' para_limits(mean, sd, N, censor=0, winsor=0, perc=0.95, cover=0.9)
#'
#' @param mean  the mean of a complete normal sample, or
#'   the intercept of the QQ regression of a censored or winsorized sample.
#' @param sd  the sd of a complete normal sample, or
#'   the slope of the QQ regresion of a censorised or winsorized sample.
#' @param N	 the sample size
#' @param censor	*optional* (default of 0) - the number of left-censored readings
#' @param winsor	*optional* (default of 0) - the number winsorized in each tail
#' @param perc	*optional* (default of 0.95) - the two-sided coverage sought
#' @param cover	*optional* (default of 0.9) - the confidence level of the CI
#'   on the reference limit
#'
#' @details This function computes two-sided reference limits and their
#' confidence intervals for data that are normal;
#' normal across the reference interval; or censored normal.
#' The reference limits are conventional mean + z*se, and their
#' confidence intervals come from the delta method.
#'
#' @returns A list containing the following components:
#'
#'   \item{lower}{the lower reference limit and its CI}
#'   \item{upper}{the upper reference limit and its CI}
#'   \item{effn}{the effective sample size with censoring or winsorization}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @examples
#' # parameters
#' mu    <- 40
#' sigma <- 10
#' n     <- 120
#' # identifying winsoring
#' wins  <- trunc(n/40)
#' # replicable randomization
#' set.seed(1069)
#' X     <- mu + sigma*rnorm(n)
#' # replicable randomization with heavy tails
#' set.seed(1069)
#' HT    <- mu + sigma * rt(n, 5)
#'
#' # visual settings
#' par(mfrow=c(2,2))
#' par(oma=c(0,0,2,0))
#' # plot to compare
#' base <- QQnorm(X, main="Base normal", showsum=TRUE)
#' title("Illustrating QQnorm with para_limits", outer=TRUE)
#' basew <- QQnorm(X, main="Winsorized", winsor=wins, showsum=TRUE)
#' ht   <- QQnorm(HT, main="Heavy tail", showsum=TRUE)
#' htw  <- QQnorm(HT, main="Winsorized", winsor=wins, showsum=TRUE)
#'
#' # evaluate and review
#' norm_results <- para_limits(mean(X), sd(X), n)
#' norm_results
#' # evaluate and review with tails
#' tailed_results <- para_limits(htw$intercept, htw$slope, n, winsor=wins)
#' tailed_results
#'
#' @references Horn PS, Peske AJ (2005). Reference intervals: a user’s guide. Washington (DC): AACC Press.
#'
#' @importFrom stats qnorm
#'
#' @export

para_limits <- function(mean, sd, N, censor=0, winsor=0, perc=0.95, cover=0.9) {
  oneper   <- (1 + perc) / 2
  effn     <- N
  cfrac    <- 1 - censor / N
  if (winsor > 0) effn <- N - 3.5 * winsor
  if (censor > 0) effn <- N / (1.38 - 0.37*cfrac)^2
  zeecov   <- abs(qnorm(0.5*(1-cover)))
  zeelim   <- abs(qnorm(oneper))
  MoE      <- sd * zeecov * sqrt(1 / effn + zeelim^2 / (2 * effn))
  lower    <- mean-zeelim*sd + c(0, -MoE, MoE)
  upper    <- mean+zeelim*sd + c(0, -MoE, MoE)
  return(list(lower=lower, upper=upper, effn=effn))
}
