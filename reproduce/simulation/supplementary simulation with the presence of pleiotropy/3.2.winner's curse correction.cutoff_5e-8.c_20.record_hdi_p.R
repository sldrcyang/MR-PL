rm(list = ls())
setwd("/public/home/yanganyi/Desktop/MR")
library(pls)
library(parallel)
library(glmnet)
library(lars)
library(hdi)


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



run_each = function(setting, cutoff, c){
  library(pls)
  library(glmnet)
  library(hdi)
  
  lasso.proj.p = data.frame(stringsAsFactors = F)
  
  for (idx in 1:total_replication){
    load(paste("/public/mig_old_storage/home1/yanganyi/Desktop/MR/stimulate/exposure-outcome/stimulation_2/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    rm(outcome_matrix)
    load(paste("stimulate/exposure-outcome/stimulation_8/u_tao",u_tao,"sd_tao",sd_tao,"/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    
    result_snp1 = result_snp[result_snp$Pr...t.. < cutoff, ]
    
    #winner's curse correction------------------
    result_snp1$abs.t = abs(result_snp1$t.value)
    if (c != 0){
      result_snp1$abs.lasso.t = result_snp1$abs.t - c / (-log10(cutoff)+log10(nrow(g_matrix0)))
      result_snp1[result_snp1$abs.lasso.t < 0, 'abs.lasso.t'] = 0
      result_snp1$corrected.p = 2 * (1 - pnorm(result_snp1$abs.lasso.t))
      result_snp1 = result_snp1[result_snp1$corrected.p < cutoff, ]
    }
    
    exposure_include = sort(unique(result_snp1$exposure))
    exposure_matrix = exposure_matrix[, exposure_include]

    iv_snp = unique(result_snp1$snp)
    if (length(iv_snp) <=1) {next}
    g_matrix = g_matrix0[, iv_snp]
    
    pls.pred = as.matrix(pred_pls(g_matrix, exposure_matrix))
    
    
    set.seed(1)
    grid = 10 ^ seq(10, -2, length=100)
    
    ################################ lasso.proj.p######################################
    outlasso = lasso.proj(x = pls.pred, y = outcome_matrix)
    p = as.vector(outlasso$pval)

    full_p = rep(1, 20)
    full_p[exposure_include] = p
    if_include = rep(FALSE, 20); if_include[exposure_include] = TRUE

    #causal effect
    tmp.p = data.frame(setting=setting, replicate=idx, image_name=1:20, lasso_proj_p=full_p, if_include=if_include)
    lasso.proj.p = rbind(lasso.proj.p, tmp.p)
  } 

  write.table(lasso.proj.p, paste('stimulate/results/winners curse/stimulation_8/u_tao',u_tao,"sd_tao",sd_tao,'/proj_p.cutoff_',cutoff,'.c_',c,'.setting_', 
                           setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  cat('setting_', setting, 'is finished.', '\n')
  return(0)
}




total_setting = 48
total_replication = 100
cutoff_set = c(5e-8)  
c_set = c(20)
u_tao = 0.2; sd_tao = 0.05


library(doParallel)
library(foreach)

cl = makeCluster(8) 
registerDoParallel(cl)

for (cutoff in cutoff_set){
  for (c in c_set){
    result = foreach(setting = 1:total_setting, .combine = 'rbind') %dopar% run_each(setting, cutoff, c)
  }
}
stopCluster(cl)
