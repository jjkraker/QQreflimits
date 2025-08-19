#' P value testing for normality of (transformed) data
#' @name BCr_Pval
#'
#' @description
#' Computes the P value of the quantile-quantile correlation coefficient (QQCC).
#'
#' @usage
#' BCr_Pval(correl, n, censor=0, winsor=0, isBC=FALSE, is2pBC=FALSE)
#'
#' @param correl the QQ correlation coefficient
#' @param n      the sample size
#' @param censor *optional* (default of 0) - the number of readings censored
#'   on the left
#' @param winsor *optional* (default of 0) - the number of readings winsorized
#'   in each tail
#' @param isBC   *optional* (default of FALSE) - if TRUE,
#'   the QQCC is after Box-Cox transformation
#' @param is2pBC *optional* (default of FALSE) - if TRUE,
#'   the QQCC is after a shifted Box-Cox transformation.
#'
#' @details
#' Lower-level function, called by other functions in package.
#' It takes information from a quantile-quantile regression, along with the
#' circumstances leading up to it, to produce a P value testing for normality.
#'
#' @returns
#'
#'   \item{Pval}{the P value}
#'
#' @author Douglas M. Hawkins, Jessica J. Kraker <krakerjj@uwec.edu>
#'
#' @example /Example/BCr_Pval_example.R
#'
#' @references Hawkins DM, Esquivel RN (2024). A Quantile-Quantile Toolbox for
#' Reference Intervals.  *The Journal of Applied Laboratory Medicine*, **9:2**, 357-370.
#'
#' @importFrom stats pnorm
#'
#' @export

BCr_Pval = function(correl, n, censor=0, winsor=0, isBC=FALSE, is2pBC=FALSE) {

  lam     <- -0.1
  cencons <- 30
  inx     <- 1

  if (isBC) inx <- 3
  if (is2pBC) inx <- 5
  if (winsor > 0) inx <- inx + 1
  if (censor > 0) {
    inx <- 1
    if (isBC) inx <- 2
  }

  # Uncensored constants

  A <- c( 1.99182,  3.1201 ,  1.40524,  2.80948,  0.92454,  2.55957)[inx]
  B <- c(-1.80221, -2.1146 , -1.78156, -2.16401, -1.7704 , -2.18845)[inx]
  D <- c( 0.67168,  0.4413 ,  0.59406,  0.42878,  0.53297,  0.3827 )[inx]
  E <- c( 0.02561,  0.08462,  0.03245,  0.07453,  0.03563,  0.07581)[inx]

  # Chebyshev censoring constants

  cmn <- rbind(c( 0.2152433,  0.6254158, -1.0697954, -3.0768692),
               c(-1.6455875,  0.6829029, -1.7539082, -2.8199981))

  csd <- rbind(c(1.25926732, -0.01723923, -12.69680564, -3.04292825),
               c(0.94621791, -0.03559444,  -7.52895893, -1.64387074))

  Pval <- NA
  fail <- censor & (winsor > 0 | is2pBC)
  if (!fail) { # Don't get P value with censoring and (winsor or 2pBC)
    Y   <- ((1-correl)^lam - 1) / lam
    lnn <- log(n + cencons)
    if (censor == 0) {
      fitmean <- A + B * lnn
      fitsd   <- D + E * lnn
    } else {
      cens    <- censor / n
      fitmean <- cmn[inx,1] + cmn[inx,2]*(cmn[inx,3]+lnn)*(cmn[inx,4]+cens)
      fitsd   <- csd[inx,1] + csd[inx,2]*(csd[inx,3]+lnn)*(csd[inx,4]+cens)
    }
    Zee     <- (Y - fitmean) / fitsd
    Pval <- pnorm(-Zee)
  }
  return(Pval)
}
