# Parametric confidence limit for a normal or normal-in-the-middle sample.
# Arguments: 
# mean:   the mean of a complete normal sample, or
#         the intercept of the QQ regression of a censored or winsorized sample.
# sd:     the sd of a complete normal sample, or
#         the slope of the QQ regresion of a censorised or winsorized sample
# N:      the sample size
# winsor:  the number winsorised in each tail
# censor:  the number of left-censored readings
# perc: the two-sided coverage sought
# cover:   the confidence level of the CI on the reference limit.
#
# Return:
# lower: the lower reference limit and its CI
# upper: the upper reference limit and its CI 
# effn:  the effective sample size with censoring or winsorization
  
para_limits <- function(mean, sd, N, winsor=0, censor=0, perc=0.95, cover=0.9) {
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
  
  