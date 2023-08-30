rm(list = ls())
setwd("/home1/yanganyi/Desktop/MR/")
library(sem)


cutoff = 5e-8
total_setting = 56
total_replication = 100

g_matrix0 = load("stimulate/genotype/g_matrix.10000sample_5000snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


###multi-tsls-iwas
for (setting in 1:total_setting){
  mse = data.frame(stringsAsFactors = F)  #y/estimate
  coef = data.frame(stringsAsFactors = F) #predicted/observed
  
  for (idx in 1:total_replication){
    load(paste("/public/mig_old_storage/home1/yanganyi/Desktop/MR/stimulate/exposure-outcome/stimulation_2/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break}  
    
    
    for (n_col in 1:ncol(exposure_matrix)){
      result_snp1 = subset(result_snp, exposure == n_col) 
      iv_snp = unique(result_snp1$snp)
      g_matrix = g_matrix0[, iv_snp]
      
      lm.fit1 = lm(exposure_matrix[,n_col] ~ g_matrix)
      exposure_matrix_predicted_ncol = predict(lm.fit1, newdata = data.frame(g_matrix))
      if (n_col == 1) {
        exposure_matrix_predicted = data.frame(exposure_matrix_predicted_ncol, stringsAsFactors = F)
        names(exposure_matrix_predicted)[n_col] = n_col
      }
      if (n_col > 1){
        exposure_matrix_predicted = cbind(exposure_matrix_predicted, exposure_matrix_predicted_ncol)
        names(exposure_matrix_predicted)[n_col] = n_col
      }
    }
    exposure_matrix_predicted = as.matrix(exposure_matrix_predicted)
    lm.fit2 = lm(outcome_matrix ~ exposure_matrix_predicted)
    beta_2sls = summary(lm.fit2)$coefficients
    
    ###1. record mse
    mse_beta = mean((beta_2sls[2:nrow(beta_2sls), 1] -  beta_xToy)^2)
    tmp.mse = data.frame(setting=setting, replicate=idx, mse_beta=mse_beta, rmse_beta=sqrt(mse_beta))
    mse = rbind(mse, tmp.mse)
    ###2.record fitted coef
    tmp.coef = data.frame(beta_2sls[2:nrow(beta_2sls), ])
    tmp.coef$true_beta = beta_xToy
    tmp.coef$setting = setting
    tmp.coef$replicate = idx
    tmp.coef$image_name = 1:20
    tmp.coef = tmp.coef[,c((ncol(tmp.coef)-2):ncol(tmp.coef), 1:(ncol(tmp.coef)-3))]
    coef = rbind(coef, tmp.coef)
    
    cat(idx, '\n')
  } 
  
  write.table(mse, paste('stimulate/results/stimulation_2/multi-tsls-iwas/multi-tsls-iwas.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef, paste('stimulate/results/stimulation_2/multi-tsls-iwas/multi-tsls-iwas.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  
  cat('setting_', setting, 'is finished.', '\n')
}