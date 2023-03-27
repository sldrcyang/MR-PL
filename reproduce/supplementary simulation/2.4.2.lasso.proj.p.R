rm(list = ls())
setwd("/home1/yanganyi/Desktop/MR")
library(pls)
library(parallel)
library(glmnet)
library(lars)
library(hdi)


cutoff = 5e-8
total_setting = 48
total_replication = 100

g_matrix0 = load("stimulate/genotype/g_matrix.ukb_10000sample_chr22_r20.6_5995snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


pred_pls = function(plsdata1_Z, plsdata2_Z){
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = 20, validation="CV", jackknife=TRUE)
  
  RMSE = RMSEP(pls.fit)
  best_ncomp = which.min(RMSE$val [2,1,-1])
  cat(best_ncomp, '\n')
  
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = best_ncomp, validation="none")
  coef = as.matrix(as.data.frame(coef(pls.fit)))
  
  pls.pred = as.data.frame(predict(pls.fit, plsdata1_Z, ncomp = best_ncomp))
  return(pls.pred)
}



run_each_setting = function(setting){
  library(pls)
  library(parallel)
  library(hdi)
  
  lasso.proj.p = data.frame(stringsAsFactors = F)
  for (idx in 1:total_replication){
    load(paste("stimulate/exposure-outcome/stimulation_4/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    
    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break}
    iv_snp = unique(result_snp$snp)
    g_matrix = g_matrix0[, iv_snp]
    
    pls.pred = as.matrix(pred_pls(g_matrix, exposure_matrix))
    
    set.seed(1)
    
    ################################ lasso.proj.p######################################
    outlasso = lasso.proj(x = pls.pred, y = outcome_matrix)
    p = as.vector(outlasso$pval)
    
    ###lasso.proj.p
    tmp.p = data.frame(setting=setting, replicate = rep(idx, 20), image_name=1:20, lasso_proj_p=p)
    lasso.proj.p = rbind(lasso.proj.p, tmp.p)
  }
    
  write.table(lasso.proj.p, paste('stimulate/results/stimulation_4/pls+lasso/pls+lasso.proj_p.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  cat('setting_', setting, 'is finished.', '\n')
  return(0)
}


library(doParallel)
library(foreach)

cl = makeCluster(14) 
registerDoParallel(cl)
result = foreach(setting = 1:total_setting, .combine = 'rbind') %dopar% run_each_setting(setting)
stopCluster(cl)






