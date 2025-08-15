
setwd("h:/research/qq_toolbox/r")
source("QQnorm.r")
source("BCR_Pval.r")
source("para_limits.r")
mu    <- 40
sigma <- 10
n     <- 120
wins  <- trunc(n/40)
dev.new()
par(mfrow=c(2,2))
par(oma=c(0,0,2,0))
set.seed(1069)
X     <- mu + sigma*rnorm(n)
base <- QQnorm(X, main="Base normal", showsum=TRUE)	
title("Illustrating QQnorm with para_limits", outer=TRUE)
para_limits(mean(X), sd(X), n)
basew <- QQnorm(X, main="Winsorized", winsor=wins, showsum=TRUE)

set.seed(1069)  
HT    <- mu + sigma * rt(n, 5)
ht   <- QQnorm(HT, main="Heavy tail", showsum=TRUE)
htw  <- QQnorm(HT, main="Winsorized", winsor=wins, showsum=TRUE)
para_limits(htw$intercept, htw$slope, n, winsor=wins)



  
  
