# Getting Box Cox transformation to normality and finding reference value
# using golden section search
# Arguments:
# X:          the data set to be transformed.
# perc:      the two-sided cover of the reference range computed.
# cover:    the confidence level of the CI computed for the reference limits.
# censor:     the number of left-censored data.
# winsor:     the number of data winsorized on each side.
# bottom:     the smallest Box-Cox power to be considered.
# top:        the largest Box-Cox power to be considered.
# epsilon:    a tolerance limit for convergence
# neff:       effective sample size, computed by the code but can be overridden.
# CI_corrfac: correction factor for CIs, computed by code but can be overridden.
#
# Returns:
# bestr:      the maximized QQ correlation coefficient.
# bestpow:    the fitted Box Cox power.
# bestxform:  the fitted Box Cox transform of the data.
# lower:      the lower reference limit and CI on the original scale.
# upper:      the upper reference limit and CI on the original scale.
# BClower:    the lower reference limits and CI on the transformed scale.
# BCupper:    the upper reference limits and CI on the transformed scale.
# meanof:     the mean of the Box-Cox transform.
# sdof:       the sd of the Box-Cox transform.
# intercept:  the intercept of the fitted QQ regression.
# slope:      the slope of the fitted QQ regression.
# Pval:       the P value of the QQ correlation. 
 
BC_limits = function(X, perc=0.95, cover=0.9, censor = 0, winsor=0, 
  bottom=-3, top=3, epsilon=0.0001, neff=NA, CI_corrfac=NA) {
  zalf   <- qnorm((1+perc)/2)	 
  cilim  <- qnorm((1+cover)/2) 
  if (sum(is.na(X)) > 0) print(paste("BC limits cases, missing", length(X), sum(is.na(X))))
  X      <- X[!is.na(X)]              # Just in case some NA slip through  
  n      <- length(X)

  if (is.na(neff)) neff <- n 

  if (is.na(CI_corrfac)) CI_corrfac <- 0.68 - 5.09/n    # n adjustment for CIs

# Optional effective sample size
  lower  <- 1 + max(c(censor, winsor))
  upper  <- n - winsor
  effn   <- neff - 3.5 * winsor  # Use this if there is no censoring
  if (censor > 0) {              # Use this for censoring, no winsorizing
    keepfrac <- 1 - censor / n
	effn  <- neff / (1.38 - 0.37*keepfrac)^2
    }
  keep   <- lower:upper
  ns     <- qnorm(((1:n)-0.5)/n)
  sorted <- sort(X)
  gr  <-  1.618
  alltry <- NULL
  allfun <- NULL
  bestr  <- -1
  flag   <- 0
  for (looper in 1:20) {
    h <- top - bottom
    if (h < epsilon) break
	c <- top    - h / gr
	d <- bottom + h / gr
    bcpowc  <- (sorted^c - 1) / c
	bcpowd  <- (sorted^d - 1) / d
	fc      <- cor(bcpowc[keep], ns[keep])  
	fd      <- cor(bcpowd[keep], ns[keep]) 
    if (is.na(fc+fd)) {
      cat(sprintf("Error in BC_Limits. Trial powers %5.5f %5.5f %6.5f %6.5f\n", c, d, fc, fd))
      flag <- flag + 1	
      if (flag == 1) {	  
	    print("Input data first fail" )
        print(sorted[keep])
	    print("c xform")
	    print(bcpowc[keep])
		print("scores")
		print(ns[keep])
        }
      if (is.na(fc)) fc <- -1     # Get away from these values
	  if (is.na(fd)) fd <- -1
	  }
	alltry  <- c(alltry, c, d)    # Log all powers tried
	allfun  <- c(allfun, fc, fd)

    if (fc > fd) {
	  top <- d
	  if (fc > bestr) {                    # Store best seen so far
	    bestr <- fc
	    bestpow <- c
	    bestxform <- (X^c-1)/c
		}
	  } else {
	  bottom <- c
	  if (fd > bestr) {                    # Store best seen so far
	    bestr <- fd
	    bestpow <- d
	    bestxform <- (X^d-1)/d
		}
	  }  
	}   
  fitlin  <- lm(sort(bestxform)[keep] ~ ns[keep])
  cofs    <- coef(fitlin)
  cut     <- cofs[1]
  slope   <- cofs[2]
  c4      <- 4*(n-1) / (4*n-3)  
  meanof  <- mean(bestxform)     # Use mean and sd of data if complete
  sdof    <- sd  (bestxform) / c4
  useeffn <- effn * CI_corrfac    
  usemean <- meanof
  usesd   <- sdof  
 
  if (censor+winsor > 0) {   # Use the QQ plot regression instead
    usemean <- cut
	usesd   <- slope
	}
  estu     <- usemean + zalf*usesd
  estl     <- usemean - zalf*usesd
  se       <- usesd * sqrt(1 / useeffn + zalf^2/(2*(useeffn-1)))   
  MoE      <- cilim * se
  offset   <- c(0, -1, 1) * MoE
 
  ci       <- c(estl+offset, estu+offset)
  BClower  <- estl + offset
  BCupper  <- estu + offset
  lower    <- (1 + bestpow * BClower) ^ (1 / bestpow)
  upper    <- (1 + bestpow * BCupper) ^ (1 / bestpow)
  
  Pval     <- BCr_Pval(bestr, n, censor=censor, winsor=winsor, isBC = TRUE)
  return(list(bestr = bestr, bestpow=bestpow, bestxform=bestxform, 
    lower=lower, upper=upper, BClower=BClower, BCupper=BCupper, Pval=Pval))
  }
