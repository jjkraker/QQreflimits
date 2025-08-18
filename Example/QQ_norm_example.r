# parameters
mu    <- 40
sigma <- 10
n     <- 120
# identifying winsoring
wins  <- trunc(n/40)
# replicable randomization
set.seed(1069)
X     <- mu + sigma*rnorm(n)
# replicable randomization with heavy tails
set.seed(1069)
HT    <- mu + sigma * rt(n, 5)

# visual settings
par(mfrow=c(2,2))
par(oma=c(0,0,2,0))
# plot to compare
base <- QQnorm(X, main="Base normal", showsum=TRUE)
title("Illustrating QQnorm with para_limits", outer=TRUE)
basew <- QQnorm(X, main="Winsorized", winsor=wins, showsum=TRUE)
ht   <- QQnorm(HT, main="Heavy tail", showsum=TRUE)
htw  <- QQnorm(HT, main="Winsorized", winsor=wins, showsum=TRUE)
