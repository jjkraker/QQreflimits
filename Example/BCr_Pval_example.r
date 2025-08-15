setwd("h:/research/qq_toolbox/r")
source("BCr_Pval.r")

BCr_Pval(c(0.993, 0.99), 120)
BCr_Pval(c(0.993, 0.99), 120, winsor=3)
BCr_Pval(c(0.993, 0.99), 120, winsor=3)
BCr_Pval(c(0.993, 0.99), 120, isBC=TRUE)
BCr_Pval(c(0.993, 0.99), 120, isBC=TRUE, winsor=3)
