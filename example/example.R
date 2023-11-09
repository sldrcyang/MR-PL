rm(list = ls())
setwd("your/path")

library(pls)
library(glmnet)
library(hdi)
source('mr_pl.R')


#load g_matrix, exposure_matrix, outcome, and association results
g_matrix0 = read.table("g_matrix.txt", header = T, stringsAsFactors = F)
exposure_matrix0 = read.table("exposure_matrix.txt", header = T, stringsAsFactors = F)
outcome = read.table("outcome.txt", header = T, stringsAsFactors = F)
gwas_assoc0 = read.table("gwas_assoc.txt", header = T, stringsAsFactors = F)

#scale data
g_matrix0 = scale(g_matrix0)
exposure_matrix0 = scale(exposure_matrix0)
outcome = scale(outcome)


#cutoff: the p-value threshold used to select instrumental snps from gwas_assoc.txt
#wcc: whether to perform winner's curse correction
#c: the tuning parameter in the winner's curse correction, and we recommended the selection for c as follows: 
    #Step 1. Determine whether the winner's curse correction is required; if the gwas_assoc.txt used to selected IVs is based on another independent dataset, then c=0;
    #Step 2. Find the grid search range for the parameter c. Here, we strongly recommend to select c as 20 based on our results;
    #Step 3. Examine the robustness of MR results at neighboring values of the selected c (e.g. c = 15 or 25).
#pleiotropy_test: whether to perform pleiotropy test


#Example —— if gwas_assoc.txt used to selected IVs is based on the same dataset
mr_result = mr_pl(g_matrix0, 
                  exposure_matrix0, 
                  outcome, 
                  gwas_assoc, 
                  cutoff = 5e-8, 
                  wcc = TRUE, 
                  c = 20,   #c=15/25
                  pleiotropy_test = TRUE)
#Example —— if gwas_assoc.txt used to selected IVs is based on another independent dataset
mr_result = mr_pl(g_matrix0, 
                  exposure_matrix0, 
                  outcome, 
                  gwas_assoc, 
                  cutoff = 5e-8, 
                  wcc = FALSE, 
                  c = 0, 
                  pleiotropy_test = TRUE)


#exposures with non-zero causal estimate and its corresponding p-value from lasso projection method
mr_result$main_results
#p value of pleiotropy test
mr_result$pleiotropy_test.p   #if p < 0.05, there exists horizontal pleiotropy
#the instrumental snps used for MR analysis (after winner's curse correction)
mr_result$iv_include
#the exposures used for MR analysis (after winner's curse correction)
mr_result$exposure_include
