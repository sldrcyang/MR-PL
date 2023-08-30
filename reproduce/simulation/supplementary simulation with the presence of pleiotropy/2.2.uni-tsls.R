rm(list = ls())
setwd("/public/home/yanganyi/Desktop/MR/")
library(sem)


cutoff = 5e-8
total_setting = 48
total_replication = 100
u_tao = 0.2; sd_tao = 0.05
dir.create(paste('stimulate/results/stimulation_8/uni-tsls/u_tao',u_tao,"sd_tao",sd_tao, sep = ''),  showWarnings = F)


g_matrix0 = load("stimulate/genotype/g_matrix.10000sample_5000snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


run_each_setting = function(setting){
  mse = data.frame(stringsAsFactors = F)  #y/estimate
  coef = data.frame(stringsAsFactors = F) #predicted/observed
  sargan = data.frame(stringsAsFactors = F)  #results of sargan test
  
  for (idx in 1:total_replication){
    load(paste("stimulate/exposure-outcome/stimulation_2/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    rm(outcome_matrix)
    load(paste("stimulate/exposure-outcome/stimulation_8/u_tao",u_tao,"sd_tao",sd_tao,"/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    

    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break}
    
    
    beta_2sls = data.frame(stringsAsFactors = F)
    tmp1.sargan = data.frame(stringsAsFactors = F)
    for (n_col in 1:ncol(exposure_matrix)){
      result_snp1 = subset(result_snp, exposure == n_col) 
      iv_snp = unique(result_snp1$snp)
      g_matrix = g_matrix0[, iv_snp]
      
      lm.fit1 = lm(exposure_matrix[,n_col] ~ g_matrix)
      exposure_matrix_predicted = predict(lm.fit1, newdata = data.frame(g_matrix))
      lm.fit2 = lm(outcome_matrix ~ exposure_matrix_predicted)

      beta_2sls.indiv = data.frame(t(summary(lm.fit2)$coefficients[2,]))
      beta_2sls.indiv$image_name = n_col
      beta_2sls = rbind(beta_2sls, beta_2sls.indiv)
      
      ############################ Sargan test ########################
      outcome_predicted = predict(lm.fit2, newdata = data.frame(exposure_matrix_predicted))
      residual = outcome_matrix - outcome_predicted
      
      reg.sargan = lm(residual ~ g_matrix)  #regress the residuals on the full set of instruments
      
      R2 = summary(reg.sargan)$r.squared
      sargan.stat = nrow(g_matrix) * R2
      df = ncol(g_matrix) - 1
      sargan.p = 1 - pchisq(sargan.stat, df = df) 
      
      tmp2.sargan = data.frame(setting=setting, replicate=idx, image_name = n_col,sargan_p = sargan.p)
      tmp1.sargan = rbind(tmp1.sargan, tmp2.sargan)
    }
    
    
    ###1.record mse
    mse_beta = mean((beta_2sls$Estimate -  beta_xToy)^2)
    tmp.mse = data.frame(setting=setting, replicate=idx, mse_beta=mse_beta, rmse_beta=sqrt(mse_beta))
    mse = rbind(mse, tmp.mse)
    
    ###2.record causal effect
    tmp.coef = data.frame(beta_2sls[,1:4])
    tmp.coef$true_beta = beta_xToy
    tmp.coef$setting = setting
    tmp.coef$replicate = idx
    tmp.coef$image_name = 1:20
    tmp.coef = tmp.coef[,c((ncol(tmp.coef)-2):ncol(tmp.coef), 1:(ncol(tmp.coef)-3))]
    coef = rbind(coef, tmp.coef)  
    
    ###3.record the p-value of sargan test
    sargan = rbind(sargan, tmp1.sargan)
    
    cat(idx, '\n')
  } 
  
  write.table(mse, paste('stimulate/results/stimulation_8/uni-tsls/u_tao',u_tao,"sd_tao",sd_tao,'/uni-tsls.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef, paste('stimulate/results/stimulation_8/uni-tsls/u_tao',u_tao,"sd_tao",sd_tao,'/uni-tsls.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(sargan, paste('stimulate/results/stimulation_8/uni-tsls/u_tao',u_tao,"sd_tao",sd_tao,'/uni-tsls.sargan.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  cat('setting_', setting, 'is finished.', '\n')
  return(0)
}


library(doParallel)
library(foreach)

cl = makeCluster(8) 
registerDoParallel(cl)
result = foreach(setting = 1:total_setting, .combine = 'rbind') %dopar% run_each_setting(setting)
stopCluster(cl)
