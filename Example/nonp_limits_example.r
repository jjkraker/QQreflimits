# parameters
mu    <- 40
sigma <- 10
n     <- 120

# replicable randomization
set.seed(1069)
X     <- mu + sigma*rnorm(n)

# evaluate and review
perc_used = 0.95
nonp_results <- nonp_limits(X, perc=perc_used)
cat("\nThe reference limits are the", 100*(1-perc_used)/2, "th percentile",
    signif(nonp_results$lower[1],5), "and the", 100*(1+perc_used)/2, "th percentile",
    signif(nonp_results$upper[1],5), ".\n")

cat("\nAnd with ", round(100*nonp_results$coverage,2),"% confidence, the lower limit is between",
    signif(nonp_results$lower[2],5), "and", signif(nonp_results$lower[3],5),
    ";\n while the upper limit is between",
    signif(nonp_results$upper[2],5),"and",signif(nonp_results$upper[3],5),".\n\n")
