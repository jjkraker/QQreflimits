# parameters
mul    <- 3.6
sigmal <- 0.75

# replicable randomization
set.seed(1069)
X      <- exp(mul + sigmal*rnorm(120))

# evaluate and review
BC_results <- BC_limits(X)
BC_results$bestpow
BC_results$bestr
# original-scale limits
BC_results$lower[1]; BC_results$upper[1]
cat("\nWith 90% [default] confidence, the lower limit is between",
    signif(BC_results$lower[2],5), "and", signif(BC_results$lower[3],5),
    ";\n while the upper limit is between",
    signif(BC_results$upper[2],5),"and",signif(BC_results$upper[3],5),".\n\n")

# adjust to have heavy tails
HT     <- X
HT[c(1,2,3,4)] <- HT[c(1,2,3,4)] * c(0.5, 0.5, 2, 2)

# evaluate and review
BC_HT_results <- BC_limits(HT)
BC_HT_results$lower; BC_HT_results$upper
# winsorized
BC_HT_wins_results <- BC_limits(HT, winsor=3)
BC_HT_wins_results$lower; BC_HT_wins_results$upper

