mr_pl = function(g_matrix0, exposure_matrix0, outcome, gwas_assoc, cutoff, wcc, c, pleiotropy_test){
  #Extract snps with gwas p-value < cutoff as instruments (significantly associated with at least one exposure)
  gwas_assoc = gwas_assoc0[gwas_assoc0$p < cutoff, ]
  
  #winner's curse correction###########################
  if (wcc == T){
    gwas_assoc$abs.z = abs(gwas_assoc$z)
    gwas_assoc$corrected.abs.z = gwas_assoc$abs.z - c / (-log10(cutoff)+log10(nrow(g_matrix0)))
    gwas_assoc[gwas_assoc$corrected.abs.z < 0, 'corrected.abs.z'] = 0
    gwas_assoc$corrected.p = 2 * (1 - pnorm(gwas_assoc$corrected.abs.z))
    gwas_assoc = gwas_assoc[gwas_assoc$corrected.p < cutoff, ]
  }

  exposure_include = sort(unique(gwas_assoc$exposure_id))
  exposure_matrix = exposure_matrix0[, exposure_include]
  
  iv_include = unique(gwas_assoc$snp_id) #instrument SNPs
  if (length(iv_include) <=1) {result = NULL}
  g_matrix = g_matrix0[, iv_include]
  

  pls.pred = as.matrix(pred_pls(g_matrix, exposure_matrix)) 
  
  
  #lasso ---
  grid = 10 ^ seq(10, -2, length=100)
  reg.mod = glmnet(pls.pred, as.matrix(outcome), alpha=1, lambda=grid)
  cv.out = cv.glmnet(pls.pred, outcome, alpha=1)
  bestlam = cv.out$lambda.min
  
  reg.coef = predict(reg.mod, type = 'coefficients', s = bestlam) 
  coef_name = reg.coef@Dimnames[[1]] [reg.coef@i +1]
  coef_val = reg.coef@x   #non-zero effect estimate
  if ("(Intercept)" %in% coef_name){
    coef_val = coef_val[-which(coef_name == '(Intercept)')]
    coef_name = coef_name[-which(coef_name == '(Intercept)')]
  }
  
  outcome_predict = predict(reg.mod, newx = pls.pred, type = 'response', s=bestlam)  #predict the outcome
  R2_outcome = cor(outcome[,1], outcome_predict[,1])   # the prediction R2 of the outcome
  
  #lasso.proj.p ---
  out.lasso.proj = lasso.proj(x = pls.pred, y = outcome)
  p = as.vector(out.lasso.proj$pval)
  p = p[match(coef_name,exposure_include)]
  
  #pleiotropy test ---
  if (pleiotropy_test == T){
    outcome.pred = predict(reg.mod, newx = pls.pred, s = bestlam)
    residual = data.frame(outcome - outcome.pred)
    
    reg.sargan = lm(residual[,1] ~ g_matrix)  #regress the residuals on the full set of instruments
    R2 = summary(reg.sargan)$adj.r.squared
    sargan.stat = nrow(g_matrix) * R2
    df = ncol(g_matrix) - ncol(pls.pred)
    sargan.p = 1 - pchisq(sargan.stat, df = df)   #if p < 0.05, there exists horizontal pleiotropy
  } else{
    sargan.stat = NULL
    sargan.p = NULL
  }
  
  #record results
  result = list()
  result[['main_results']] = data.frame(exposure_name=coef_name, causal_estimate = coef_val, lasso_proj_p = p)
  result[['wcc']] = wcc
  result[['c']] = c
  result[['cutoff']] = cutoff
  result[['pleiotropy_test.stat']] = sargan.stat
  result[['pleiotropy_test.p']] = sargan.p
  result[['iv_include']] = iv_include
  result[['exposure_include']] = exposure_include
  result[['R2']] = R2_outcome
  
  return(result)
}


pred_pls = function(plsdata1_Z, plsdata2_Z){
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = ncol(plsdata2_Z), validation="CV", jackknife=TRUE)
  RMSE = RMSEP(pls.fit) 
  best_ncomp = which.min(RMSE$val [2,1,-1]) 
  
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = best_ncomp, validation="none") 
  coef = as.matrix(as.data.frame(coef(pls.fit)))
  
  pls.pred = as.data.frame((predict(pls.fit, plsdata1_Z, ncomp = best_ncomp)))
  names(pls.pred) = gsub(paste('.',best_ncomp, ' comps',sep=''), '', names(pls.pred)) 
  return(pls.pred)
}

