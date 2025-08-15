#  Calculate shortest nonparametric reference limit or limit of blank with CI 
# Arguments:
# X:     the data set whose reference limit is sought.
# RR:    the way percentiles are defined:
#        FALSE -- use the Weibull limit i/(n+1), 
#        TRUE  -- use the Hazen limit (i-0.5)/n.
# perc:  the two-sided probability
# cover: the confidence level of the CI for the reference limit.
#  
# Returns:
#
# lower:    the lower limit and confidence interval. 
# upper:    the upper limit and confidence interval.
# a:        the index of the order statistic defining the lower limit of the CI.
# b:        the index of the order statistic defining the upper limit of the CI.
# coverage: the actual confidence level the interval achieves.

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
