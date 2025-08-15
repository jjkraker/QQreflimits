
setwd("h:/research/qq_toolbox/r")
source("QQnorm.r")
source("BCR_Pval.r")
source("para_limits.r")
source("nonp_limits.r")
mu    <- 40
sigma <- 10
n     <- 120
set.seed(1069)
X     <- mu + sigma*rnorm(n)
nonp_limits(X)


  
  
