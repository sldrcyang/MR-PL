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
#c: the tuning parameter in the winner's curse correction, which is recommended to set as 20; if wcc = FALSE, then c = 0
#pleiotropy_testwhether to perform pleiotropy test

#if gwas_assoc.txt used to selected IVs is based on the same dataset
mr_result = mr_pl(g_matrix0, 
                  exposure_matrix0, 
                  outcome, 
                  gwas_assoc, 
                  cutoff = 5e-8, 
                  wcc = TRUE, 
                  c = 20, 
                  pleiotropy_test = TRUE)
#if gwas_assoc.txt used to selected IVs is based on another independent dataset
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
