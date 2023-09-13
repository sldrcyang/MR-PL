rm(list = ls())
setwd("/public/home/yanganyi/Desktop/MR")
library(pls)
library(parallel)
library(glmnet)
library(lars)


g_matrix0 = load("stimulate/genotype/g_matrix.10000sample_5000snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


pred_pls = function(plsdata1_Z, plsdata2_Z){
  if (ncol(plsdata1_Z) > 20) {
    pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = 20, validation="CV", jackknife=TRUE)
    RMSE = RMSEP(pls.fit)
    best_ncomp = which.min(apply(RMSE$val [2,1:20,-1], 2, mean))
  }
  else {best_ncomp = ncol(plsdata1_Z)}
  
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = best_ncomp, validation="none")
  coef = as.matrix(as.data.frame(coef(pls.fit)))
  
  pls.pred = as.data.frame((predict(pls.fit, plsdata1_Z, ncomp = best_ncomp)))
  names(pls.pred) = gsub(paste('.',best_ncomp, ' comps',sep=''), '', names(pls.pred))
  return(pls.pred)
}


run_each = function(setting, cutoff, c){
  library(pls)
  #library(parallel)
  library(glmnet)
  #library(lars)
  
  mse.2 = data.frame(stringsAsFactors = F); coef.2 = data.frame(stringsAsFactors = F)
  sargan.2 = data.frame(stringsAsFactors = F)
  
  for (idx in 1:total_replication){
    load(paste("stimulate/exposure-outcome/stimulation_2/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
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
    
    ################################ 2.lasso ######################################
    reg.mod = glmnet(pls.pred, outcome_matrix, alpha=1, lambda=grid)  
    cv.out = cv.glmnet(pls.pred, outcome_matrix, alpha=1)
    bestlam = cv.out$lambda.min
    
    reg.coef = predict(reg.mod, type = 'coefficients', s = bestlam)
    coef_name = reg.coef@Dimnames[[1]] [reg.coef@i +1]
    coef_val = reg.coef@x  
    if ("(Intercept)" %in% coef_name){
      select_exposure_idx = coef_name[-1]
      select_exposure_idx = as.numeric(gsub('Y', '', select_exposure_idx))
    } else{
      select_exposure_idx = as.numeric(gsub('Y', '', select_exposure_idx))
    }
    
    #lasso.proj.p---
    beta = rep(0, length(exposure_include))
    beta[ select_exposure_idx ] = coef_val[-1]
    
    outcome_predict = predict(reg.mod, newx = pls.pred, s = bestlam)  #predict the outcome
    residuals = outcome_matrix[,1] - outcome_predict[,1]
    sd = sd(residuals)
    
    outlasso = lasso.proj(x = pls.pred, y = outcome_matrix, betainit= beta, sigma = sd)
    p = as.vector(outlasso$pval)
    
    
    ##Sargan test ##----------------------------------
    outcome_predict = predict(reg.mod, s = bestlam, newx = pls.pred)
    residuals = outcome_matrix - outcome_predict
    
    reg.sargan = lm(residuals ~ g_matrix)  #regress the residuals on the full set of instruments
    
    R2 = summary(reg.sargan)$r.squared
    sargan.stat = nrow(g_matrix) * R2
    df = ncol(g_matrix) - ncol(pls.pred)
    sargan.p = 1 - pchisq(sargan.stat, df = df)
    #----------------------------------------------------------------------
    
    
    ###1).mse
    full_coef_val = rep(0, 20)
    full_coef_val[exposure_include] = 0
    full_coef_val[ exposure_include[select_exposure_idx] ] = coef_val[-1]
    if_include = rep(FALSE, 20); if_include[exposure_include] = TRUE
    
    mse_beta = mean((full_coef_val - beta_xToy) ^ 2)
    mse_beta_analyse = mean((full_coef_val[exposure_include] - beta_xToy[exposure_include]) ^ 2)
    tmp.mse = data.frame(setting=setting, replicate=idx, cutoff=cutoff, c=c, 
                         mse_beta=mse_beta, mse_beta_analyse=mse_beta_analyse)
    mse.2 = rbind(mse.2, tmp.mse)
    
    ###2).causal effect
    tmp.coef = data.frame(setting=setting, replicate=idx, image_name=1:20, coef=full_coef_val, if_include=if_include)
    tmp.coef$true_beta = beta_xToy
    coef.2 = rbind(coef.2, tmp.coef)
    
    ###3).sargan
    tmp.sargan.2 = data.frame(setting=setting, replicate=idx, sargan_p = sargan.p)
    sargan.2 = rbind(sargan.2, tmp.sargan.2)
    
    ###4)p
    full_p = rep(NA, 20)
    full_p[exposure_include] = p
    if_include = rep(FALSE, 20); if_include[exposure_include] = TRUE
    tmp.p = data.frame(setting=setting, replicate = rep(idx, 20), image_name=1:20, lasso_proj_p=full_p, if_include=if_include)
    lasso.proj.p = rbind(lasso.proj.p, tmp.p)
  } 
  
  write.table(mse.2, paste('stimulate/results/winners curse/stimulation_8/u_tao',u_tao,"sd_tao",sd_tao,'/mse.cutoff_',cutoff,'.c_',c,'.setting_', 
                           setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef.2, paste('stimulate/results/winners curse/stimulation_8/u_tao',u_tao,"sd_tao",sd_tao,'/coef.cutoff_',cutoff,'.c_',c,'.setting_', 
                           setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(sargan.2, paste('stimulate/results/winners curse/stimulation_8/u_tao',u_tao,"sd_tao",sd_tao,'/sargan.cutoff_',cutoff,'.c_',c,'.setting_', 
                            setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(lasso.proj.p, paste('stimulate/results/winners curse/stimulation_8/proj_p/proj_p.cutoff_',cutoff,'.c_',c,'.setting_', 
                                  setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  cat('setting_', setting, 'is finished.', '\n')
  return(0)
}



total_replication = 100
total_setting = 48
cutoff_set = c(5e-8) 
c = 20
u_tao = 0.2; sd_tao = 0.05
dir.create(paste('/public/home/yanganyi/Desktop/MR/stimulate/results/winners curse/stimulation_8/u_tao',u_tao,'sd_tao',sd_tao, sep = ''),  showWarnings = F)


for (setting in 1:total_setting){
  for (cutoff in cutoff_set){
    result = run_each(setting, cutoff, c)
    cat('setting_',setting, ', cutoff=', cutoff, 'is OK.\n')
  }
}
