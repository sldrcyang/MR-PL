rm(list = ls())
setwd("/home1/yanganyi/Desktop/MR/")
library(sem)


cutoff = 5e-8
total_setting = 48
total_replication = 100

g_matrix0 = load("stimulate/genotype/g_matrix.ukb_10000sample_chr22_r20.6_5995snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


###multi-tsls
for (setting in 1:total_setting){
  mse = data.frame(stringsAsFactors = F)  #y/estimate
  coef = data.frame(stringsAsFactors = F) #predicted/observed
  
  for (idx in 1:total_replication){
    load(paste("stimulate/exposure-outcome/stimulation_4/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break}
    iv_snp = unique(result_snp$snp)
    g_matrix = g_matrix0[, iv_snp]
    
    
    tsls = tsls(outcome_matrix ~ exposure_matrix, ~ g_matrix)
    beta_2sls = summary(tsls)$coefficients
    
    ###1.mse
    mse_beta = mean((beta_2sls[2:nrow(beta_2sls), 1] -  beta_xToy)^2)
    tmp.mse = data.frame(setting=setting, replicate=idx, mse_beta=mse_beta, rmse_beta=sqrt(mse_beta))
    mse = rbind(mse, tmp.mse)
    ###2.coef
    tmp.coef = data.frame(beta_2sls[2:nrow(beta_2sls), ])
    tmp.coef$true_beta = beta_xToy
    tmp.coef$setting = setting
    tmp.coef$replicate = idx
    tmp.coef$image_name = 1:20
    tmp.coef = tmp.coef[,c((ncol(tmp.coef)-2):ncol(tmp.coef), 1:(ncol(tmp.coef)-3))]
    coef = rbind(coef, tmp.coef)
    
    cat(idx, '\n')
  } 
  
  write.table(mse, paste('stimulate/results/stimulation_4/multi-tsls/multi-tsls.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef, paste('stimulate/results/stimulation_4/multi-tsls/multi-tsls.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  
  cat('setting_', setting, 'is finished.', '\n')
}
