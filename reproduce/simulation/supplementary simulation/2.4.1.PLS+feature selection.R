rm(list = ls())
setwd("/home1/yanganyi/Desktop/MR")
library(pls)
library(parallel)
library(glmnet)
library(lars)


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
  
  pls.pred = as.data.frame((predict(pls.fit, plsdata1_Z, ncomp = best_ncomp)))
  return(pls.pred)
}



run_each_setting = function(setting){
  library(pls)
  library(parallel)
  library(glmnet)
  library(lars)
  
  mse.1 = data.frame(stringsAsFactors = F); coef.1 = data.frame(stringsAsFactors = F)
  mse.2 = data.frame(stringsAsFactors = F); coef.2 = data.frame(stringsAsFactors = F)
  mse.3 = data.frame(stringsAsFactors = F); coef.3 = data.frame(stringsAsFactors = F)
  mse.4 = data.frame(stringsAsFactors = F); coef.4 = data.frame(stringsAsFactors = F) 
  
  
  for (idx in 1:total_replication){
    load(paste("stimulate/exposure-outcome/stimulation_4/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    
    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break}
    iv_snp = unique(result_snp$snp)
    g_matrix = g_matrix0[, iv_snp]
    
    pls.pred = as.matrix(pred_pls(g_matrix, exposure_matrix))
    
    
    set.seed(1)
    grid = 10 ^ seq(10, -2, length=100)
    
    ################################ 1.ridge ######################################
    reg.mod = glmnet(pls.pred, outcome_matrix, alpha=0, lambda=grid)
    cv.out = cv.glmnet(pls.pred, outcome_matrix, alpha=0)
    bestlam = cv.out$lambda.min
    reg.coef = predict(reg.mod, type = 'coefficients', s = bestlam) 
    
    coef_name = reg.coef@Dimnames[[1]]
    coef_val = reg.coef@x
    
    ###1).mse
    mse_beta = mean((coef_val[-1] -  beta_xToy)^2)
    tmp.mse = data.frame(setting=setting, replicate = idx, mse_beta = mse_beta, rmse_beta = sqrt(mse_beta))
    mse.1 = rbind(mse.1, tmp.mse)
    ###2).causal effect
    tmp.coef = data.frame(setting=setting, replicate = rep(idx, 20), image_name=coef_name[-1], coef=coef_val[-1])
    tmp.coef$true_beta = beta_xToy
    coef.1 = rbind(coef.1, tmp.coef)
    
    
    ################################ 2.lasso ######################################
    reg.mod = glmnet(pls.pred, outcome_matrix, alpha=1, lambda=grid)
    cv.out = cv.glmnet(pls.pred, outcome_matrix, alpha=1)
    bestlam = cv.out$lambda.min
    
    reg.coef = predict(reg.mod, type = 'coefficients', s=bestlam)
    coef_name = reg.coef@Dimnames[[1]] [reg.coef@i +1]
    coef_val = reg.coef@x
    if ("(Intercept)" %in% coef_name){
      select_exposure_idx = match(coef_name, reg.coef@Dimnames[[1]]) - 1
      select_exposure_idx = select_exposure_idx[-1]
    } else{
      select_exposure_idx = match(coef_name, reg.coef@Dimnames[[1]])
    }
    
    ###1).mse
    full_coef_val = rep(0, 20); full_coef_val[select_exposure_idx] = coef_val[-1]
    mse_beta = mean((full_coef_val - beta_xToy) ^ 2)
    tmp.mse = data.frame(setting=setting, replicate = idx, mse_beta = mse_beta, rmse_beta = sqrt(mse_beta))
    mse.2 = rbind(mse.2, tmp.mse)
    ###2).causal effect
    tmp.coef = data.frame(setting=setting, replicate = rep(idx, 20), image_name=1:20, coef=full_coef_val)
    tmp.coef$true_beta = beta_xToy
    coef.2 = rbind(coef.2, tmp.coef)
    
    
    ################################ 3.Elastic net ######################################
    alphalist = seq(0.05, 0.95, by=0.05)
    elasticnet.cv.out = lapply(alphalist, function(a){
      cv.out = cv.glmnet(pls.pred, outcome_matrix, alpha=a)
      bestlam = cv.out$lambda.min
      minmse = cv.out[["cvm"]] [match(bestlam, cv.out[["lambda"]])]
      c(a, bestlam, minmse)
    })
    elasticnet.cv.out = data.frame(t(data.frame(elasticnet.cv.out)))
    bestalpha = elasticnet.cv.out[which.min(elasticnet.cv.out[,3]),1]
    bestlam = elasticnet.cv.out[which.min(elasticnet.cv.out[,3]),2]
    
    reg.mod = glmnet(pls.pred, outcome_matrix, alpha=bestalpha, lambda=grid)
    reg.coef = predict(reg.mod, type = 'coefficients', s=bestlam)
    coef_name = reg.coef@Dimnames[[1]] [reg.coef@i +1]
    coef_val = reg.coef@x  
    
    if ("(Intercept)" %in% coef_name){
      select_exposure_idx = match(coef_name, reg.coef@Dimnames[[1]]) - 1
      select_exposure_idx = select_exposure_idx[-1]
    } else{
      select_exposure_idx = match(coef_name, reg.coef@Dimnames[[1]])
    }
    
    ###1).mse
    full_coef_val = rep(0, 20); full_coef_val[select_exposure_idx] = coef_val[-1]
    mse_beta = mean((full_coef_val - beta_xToy) ^ 2)
    tmp.mse = data.frame(setting=setting, replicate = idx, mse_beta = mse_beta, rmse_beta = sqrt(mse_beta))
    mse.3 = rbind(mse.3, tmp.mse)
    ###2).causal effect
    tmp.coef = data.frame(setting=setting, replicate = rep(idx, 20), image_name=1:20, coef=full_coef_val)
    tmp.coef$true_beta = beta_xToy
    coef.3 = rbind(coef.3, tmp.coef)
    
    
    ################################ 4.LARS ######################################
    cv.out = try(cv.lars(pls.pred, outcome_matrix,  K=10, plot.it = F, se = T, type = "lar", mode="step"))
    if("try-error" %in% class(cv.out))  {next}
    beststep = cv.out$index[which.min(cv.out$cv)]
    
    reg.mod = lars(pls.pred, outcome_matrix, type = "lar")
    reg.coef = reg.mod$beta [beststep, ]
    
    ###1).mse
    mse_beta = mean((reg.coef -  beta_xToy)^2)
    tmp.mse = data.frame(setting=setting, replicate = idx, mse_beta = mse_beta, rmse_beta = sqrt(mse_beta))
    mse.4 = rbind(mse.4, tmp.mse)
    ###2).causal effect
    tmp.coef = data.frame(setting=setting, replicate = rep(idx, 20), image_name=1:20, coef=reg.coef)
    tmp.coef$true_beta = beta_xToy
    coef.4 = rbind(coef.4, tmp.coef)
    
    
    cat(idx, '\n')
  }
  
  
  write.table(mse.1, paste('stimulate/results/stimulation_4/pls+ridge/pls+ridge.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef.1, paste('stimulate/results/stimulation_4/pls+ridge/pls+ridge.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  write.table(mse.2, paste('stimulate/results/stimulation_4/pls+lasso/pls+lasso.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef.2, paste('stimulate/results/stimulation_4/pls+lasso/pls+lasso.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  write.table(mse.3, paste('stimulate/results/stimulation_4/pls+elasticnet/pls+elasticnet.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef.3, paste('stimulate/results/stimulation_4/pls+elasticnet/pls+elasticnet.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  write.table(mse.4, paste('stimulate/results/stimulation_4/pls+lars/pls+lars.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef.4, paste('stimulate/results/stimulation_4/pls+lars/pls+lars.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  cat('setting_', setting, 'is finished.', '\n')
  
  return(0)
}



library(doParallel)
library(foreach)

cl = makeCluster(15) 
registerDoParallel(cl)
result = foreach(setting = 1:total_setting, .combine = 'rbind') %dopar% run_each_setting(setting)
stopCluster(cl)


