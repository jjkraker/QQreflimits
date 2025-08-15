# formerly "qqdoug" makes QQ plot of complete or censored data,
# fits regression line to complete uncensored portion of data,
# calculates the QQ correlation coefficient of the fitted line,
# reports the P value of the QQCC as calculated by BCR_pval.
# Arguments:
# X:    the numeric data vector to be plotted.  Censored values should 
#          be reported as the censoring value.
# main:    the header text.
# ylab:    the name of the variable.
# censor:  the number of data censored on the left.
# winsor:  the number of data winsored in each tail.
# joinem:  if TRUE the plot is drawn as a segmented line, if FALSE the  
#          individual points are plotted as x if winsorized, else *.
# ylim:    optional y limits on the plot.
# isBC:    if true, the data set is a Box-Cox transform of the original data.
# is2pBC:  if true, the data set is a shifted Box-Cox tranform.
#          These two parameters are relevant to the P value calculation.
# doplot:  if true, the QQ plot is drawn.
# showP:   if true, the QQCC P value is shown on the plot. 
# fitline: if true, the QQ regression line is plotted 
# showsum: if true, the intercept and slope of the QQ regression are shown.

# Returns:
# correl:    the QQ correlation coefficient.
# Pval:      the P value of the QQCC.
# mean:      the mean of the data.
# sd:        the sd of the data.
# intercept: the intercept of the QQ regression line, used in place of
#            the mean when there is censoring or winsorization. 
# slope:     the slope of the QQ regression line, used in place of
#            the sd when there is censorin or winsorization.

QQnorm = function(X, main="", ylab="", censor=0, winsor=0,
  joinem=FALSE, ylim=c(NA,NA), isBC=FALSE, is2pBC=FALSE, doplot=TRUE, 
  showP=TRUE, fitline=TRUE, showsum=FALSE) {
  mean    <- NA
  sd      <- NA
  if (censor + winsor == 0) {  # Don't get mean or sd if censored or winsorized
    mean    <- mean(X)
    sd      <- sd  (X)
	}
  winsor  <- trunc(winsor)	  
  n       <- length(X)
  pch     <- rep(1, n)
  if (winsor > 0) {
    pch[1:winsor] <- 4
	pch[(n-winsor+1):n] <- 4
	}
  ns      <- qnorm(((1:n)-0.5)/n)
  X    <- sort(X)
  range   <- (censor+1) : n
  QQrange <- (max(c(censor,winsor))+1) : (n-winsor)
  ltype   <- "p"
  if (joinem) ltype = "l"
  if (is.na(sum(ylim))) ylim <- c(min(X), max(X))
  correl  <- cor(ns[QQrange], X[QQrange])
  Pval    <- BCr_Pval(correl, n, censor, winsor, isBC, is2pBC) 
  leg     <- paste("r=", round(correl,4))
  leg2    <- sprintf("P=%6.5f", Pval)
  linef   <- lm(X[QQrange] ~ ns[QQrange])
  if (doplot) {
	plot(ns[range], X[range], main=main, xlab="Normal scores", 
      ylab=ylab, type=ltype, ylim=ylim, pch=pch)
	if (showP) leg <- paste(leg, leg2, sep="\n")  
	legend("topleft", leg, bty="n")

	if (fitline) abline(linef)
	bl   <- ""
	summ <- sprintf("cut %4.4g\nslope  %4.4g",
	  coef(linef)[1], coef(linef)[2])
	if (showsum) bl <- paste(bl, summ, sep="")
	if (winsor > 0) {
	  vline <- 0.5*(ns[winsor]+ns[winsor+1])
	  abline(v= vline, lty=3)
	  abline(v=-vline, lty=3)
	  }
	legend("right", legend=bl, bty="n")  
	}	
  return(list(correl=correl, Pval=Pval, mean=mean, sd=sd, 
    intercept=coef(linef)[1], slope=coef(linef)[2]))
  }
 