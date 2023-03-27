rm(list = ls())
setwd("your/path")
library(pls)
library(parallel)
library(glmnet)
library(lars)
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


#if gwas summary statistics used to selected IVs is based on another independent dataset, then wcc = F
mr_result = mr_pl(g_matrix0, exposure_matrix0, outcome, gwas_assoc, cutoff=5e-8, wcc = T, c=20, pleiotropy_test = T)

#exposures with non-zero causal estimate and its corresponding p-value from lasso projection method
mr_result$main_results
#p value of pleiotropy test
mr_result$pleiotropy_test.p   #if p < 0.05, there exists horizontal pleiotropy
#the instrumental snps used for MR analysis (after winner's curse correction)
mr_result$iv_include
#the exposures used for MR analysis (after winner's curse correction)
mr_result$exposure_include

