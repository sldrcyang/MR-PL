rm(list = ls())
setwd("/public/home/yanganyi/Desktop/MR")
library(pls)
library(parallel)
library(glmnet)
library(lars)
library(hdi)


cutoff = 5e-8
total_setting = 48
total_replication = 100
u_tao = 0.2; sd_tao = 0.05

g_matrix0 = load("stimulate/genotype/g_matrix.10000sample_5000snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


pred_pls = function(plsdata1_Z, plsdata2_Z){
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = 20, validation="CV", jackknife=TRUE)
  
  RMSE = RMSEP(pls.fit)
  best_ncomp = which.min(RMSE$val [2,1,-1])
  cat(best_ncomp, '\n')
  
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = best_ncomp, validation="none")
  coef = as.matrix(as.data.frame(coef(pls.fit))) 
  
  pls.pred = as.data.frame((predict(pls.fit, plsdata1_Z, ncomp = best_ncomp)))
  return(pls.pred)
}




run_each_setting = function(setting){
  library(pls)
  library(parallel)
  library(hdi)
  
  lasso.proj.p = data.frame(stringsAsFactors = F)
  
  for (idx in 1:total_replication){
    load(paste("/public/mig_old_storage/home1/yanganyi/Desktop/MR/stimulate/exposure-outcome/stimulation_2/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    rm(outcome_matrix)
    load(paste("stimulate/exposure-outcome/stimulation_8/u_tao",u_tao,"sd_tao",sd_tao,"/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    
    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break} 
    iv_snp = unique(result_snp$snp)
    g_matrix = g_matrix0[, iv_snp]
    
    pls.pred = as.matrix(pred_pls(g_matrix, exposure_matrix)) 
    
    set.seed(1)
    grid = 10 ^ seq(10, -2, length=100)
    
    ################################ lasso.proj.p######################################
    outlasso = lasso.proj(x = pls.pred, y = outcome_matrix)  
    p = as.vector(outlasso$pval)
    
    tmp.p = data.frame(setting=setting, replicate = rep(idx, 20), image_name=1:20, lasso_proj_p=p)
    lasso.proj.p = rbind(lasso.proj.p, tmp.p)
  }
    
  write.table(lasso.proj.p, paste('stimulate/results/stimulation_8/pls+lasso/u_tao',u_tao,"sd_tao",sd_tao,'/pls+lasso.proj_p.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  cat('setting_', setting, 'is finished.', '\n')
  return(0)
}


library(doParallel)
library(foreach)

cl = makeCluster(16) 
registerDoParallel(cl)
result = foreach(setting = 1:total_setting, .combine = 'rbind') %dopar% run_each_setting(setting)
stopCluster(cl)

