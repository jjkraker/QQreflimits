# compute the Pvalue for two QQCC's
BCr_Pval(c(0.993, 0.99), 120)
# if censored
BCr_Pval(c(0.993, 0.99), 120, censor=3)
# if winsorized
BCr_Pval(c(0.993, 0.99), 120, winsor=3)
# on Box-Cox transformed data
BCr_Pval(c(0.993, 0.99), 120, isBC=TRUE)
# on Box-Cox transformed data, and winsorized
BCr_Pval(c(0.993, 0.99), 120, isBC=TRUE, winsor=3)
