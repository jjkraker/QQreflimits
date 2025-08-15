setwd("h:/research/toolbox/r")
source("BC_limits.r")

set.seed(1069)
mul    <- 3.6
sigmal <- 0.75
X      <- exp(mul + sigmal*rnorm(120))
BC_limits(X)
HT     <- X
HT[c(1,2,3,4)] <- HT[c(1,2,3,4)] * c(0.5, 0.5, 2, 2)
BC_limits(HT)
BC_limits(HT, winsor=3)

