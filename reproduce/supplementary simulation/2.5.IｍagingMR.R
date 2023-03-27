rm(list = ls())
setwd("F:/类脑/1_实操文件/MR/")
library(data.table)
library(matlabr)
library(R.matlab) 
library(MendelianRandomization) 
library(psych) 
library(psychTools)
library(factoextra)


cutoff = 5e-8
total_setting = 48
total_replication = 100

g_matrix0 = load("stimulate/genotype/g_matrix.ukb_10000sample_chr22_r20.6_5995snp.RData")
g_matrix0 = eval(parse(text = g_matrix0))
rm(g_matrix)


#### linear regression of a factor and single SNP
factor_snp.reg = function(factor, dataX){
  factor_snp.lm = matrix(NA, nrow = ncol(dataX), ncol = 4) 
  colnames(factor_snp.lm) = c("Beta_exp", "Beta_SE_exp", "Tvalue_exp", "Pvalue_exp")
  for (i in 1:ncol(dataX)) {
    fm = summary(lm(factor ~ dataX[,i] ))
    factor_snp.lm[i,1] = fm$coefficients[2,1]
    factor_snp.lm[i,2] = fm$coefficients[2,2]
    factor_snp.lm[i,3] = fm$coefficients[2,3]
    factor_snp.lm[i,4] = fm$coefficients[2,4]
  }
  return(factor_snp.lm)
}


#FUNCTION: Run Factor Analysis
fa.mr <- function(dataM, dataX, las_sel_col, ulmres.out) {
  if (length(las_sel_col) == 1){
    factor = dataM
    factor_snp.lm = factor_snp.reg(factor, dataX)
    #IVW
    mmr_gs_sub = mr_allmethods(mr_input(bx = factor_snp.lm[,1], bxse = factor_snp.lm[,2],
                                        by = ulmres.out[,1], byse = ulmres.out[,2]),  method = "all")
    beta_exposure = rep(0, 20)
    beta_exposure[las_sel_col] = mmr_gs_sub@Values[4,2]
  } 
  
  
  else{
    searchN = fa.parallel(dataM, main="Parallel Analysis", plot=F)
    factorN = searchN[["nfact"]]
    
    if (is.na(factorN)){factorN = 1}
    
    #Factor Analysis
    fa_results = fa(dataM, nfactors = factorN, scores="tenBerge", rotate = 'none') 
    factors = fa_results[["scores"]]
    factors_loadings = fa_results[["loadings"]][,]
    
    #causal effect
    factors_beta = c()
    for (n in 1:factorN){
      factor_snp.lm = factor_snp.reg(factors[,n], dataX)
      #IVW
      mmr_gs_sub = mr_allmethods(mr_input(bx = factor_snp.lm[,1], bxse = factor_snp.lm[,2],
                                          by = ulmres.out[,1], byse = ulmres.out[,2]),  method = "all")
      factors_beta = c(factors_beta, mmr_gs_sub@Values[4,2])
    }
    
    beta_exposure = rep(0, 20)
    if (factorN == 1){
      beta_exposure[las_sel_col] = factors_loadings * factors_beta
    } else{
      beta_exposure[las_sel_col] = factors_loadings %*% factors_beta} 
    
  }
  
  return(beta_exposure)
}



save_imagingMR_LAS_results <- function(setting, idx, dataX, dataY, dataM, p, q, alpha){
  ######################### 2. MR-single IV #########################
  ## 2a. mediator ~ SNP regression
  ulmres.med <- vector("list", q)
  for (j in 1:q) {
    med = dataM[,j]
    temp.res = matrix(NA, p, 4) 
    colnames(temp.res) = c("Beta_exp", "Beta_SE_exp","Tvalue_exp", "Pvalue_exp")
    for (i in 1:p) {
      fm = summary(lm(med ~ dataX[,i]))
      temp.res[i,1] = fm$coefficients[2,1]
      temp.res[i,2] = fm$coefficients[2,2]
      temp.res[i,3] = fm$coefficients[2,3]
      temp.res[i,4] = fm$coefficients[2,4]
    }
    
    ulmres.med[[j]] = temp.res
    names(ulmres.med)[[j]] = colnames(dataM)[j]
  }
  
  snp_med_pval <- sapply(ulmres.med, function(x){x[,4]})
  snp_med_neglogP <- apply(snp_med_pval, 2, function(x){-log10(x)})
  
  ## 2b. outcome ~ SNP regression
  ulmres.out = matrix(NA, p, 4) 
  colnames(ulmres.out) = c("Beta_out", "Beta_SE_out","Tvalue_out", "Pvalue_out")
  for (i in 1:ncol(dataX)) {
    fm = summary(lm(dataY ~ dataX[,i]))
    ulmres.out[i,1] = fm$coefficients[2,1]
    ulmres.out[i,2] = fm$coefficients[2,2]
    ulmres.out[i,3] = fm$coefficients[2,3]
    ulmres.out[i,4] = fm$coefficients[2,4]
  }
  
  
  ## 2c. univariate MR for each SNP-mediator-outcome pair 
  uniMR.1snp = c() 
  uniMR.1snp.pvalue <- matrix(NA, p, q)
  uniMR.1snp.fullresults <- vector("list", p) ## MR fitting results for all SNP vs. all fa
  for (h in 1:p) {  ## MR using single SNP, h:snp的index
    uniMR = data.frame(matrix(NA, q, 5))
    colnames(uniMR) = c("Estimate","StdError","CILower","CIUpper","Pvalue")
    rownames(uniMR) = names(ulmres.med)
    uniMR.fullresults <- vector("list", q) ## MR fitting results for one SNP vs. all fa
    for (k in 1:q) {
      MRobj = mr_ivw(mr_input(bx = ulmres.med[[k]][h,1], bxse=ulmres.med[[k]][h,2],
                              by = ulmres.out[h,1], byse = ulmres.out[h,2]))  #IVW方法
      uniMR[k,1] = MRobj$Estimate
      uniMR[k,2] = MRobj$StdError
      uniMR[k,3] = MRobj$CILower
      uniMR[k,4] = MRobj$CIUpper
      uniMR[k,5] = MRobj$Pvalue
      
      uniMR.fullresults[[k]] = MRobj
      names(uniMR.fullresults)[k] = names(ulmres.med)[k]
    }
    uniMR.1snp[[h]] = uniMR
    uniMR.1snp.pvalue[h,] <- uniMR$Pvalue
    uniMR.1snp.fullresults[[h]] = uniMR.fullresults
  }
  
  uniMR.1snp.neglogP <- apply(uniMR.1snp.pvalue, 2, function(x){-log10(x)}) 
  colnames(uniMR.1snp.neglogP) <- colnames(dataM) 
  
  ## B-H correction 
  uniMR_bh <- matrix(p.adjust(uniMR.1snp.pvalue, method="BH"), p, q)
  uniMR_bh_cut <- (uniMR_bh < alpha)
  uniMR_sel_col <- which(colSums(uniMR_bh_cut)>0)
  uniMR_sel_row <- which(rowSums(uniMR_bh_cut)>0)
  
  comb2_neglogP <- uniMR.1snp.neglogP ## combine two p-values
  comb2_neglogP[snp_med_pval>alpha] <- 0 

  comb2_neglogP <- comb2_neglogP + matrix(runif(p*q,0,1e-10), p, q)

  
  
  ######################### 3. LAS submatrix #########################
  ## save the data for passing to Matlab 
  write.table(comb2_neglogP, row.names=F, col.names=F, sep=",",
              file=paste("MR_software/ImagingMR-LAS/trial_results/las_input.setting_",setting,".", idx,".csv", sep="")) #不区分的话并行会打架
  
  ## use LAS matlab function in the mtba package 
  code <- c("cd 'F:/类脑/1_实操文件/MR/MR_software/ImagingMR-LAS/trial_results/';", 
            paste("W=readtable('las_input.setting_", setting,".", idx, ".csv');", sep=""),
            "W1=table2array(W);",
            "addpath('F:/类脑/2_数据报告资料/软件包/bi-cluster_matlab_toolbox.win_nix/mtba_win');",
            "res=LAS(W1, 1);", 
            "rowSel=res.RowxNum;",
            "colSel=res.NumxCol;",
            paste("save('sub_matrix.setting_",setting,".mat','rowSel','colSel','res')", sep=""))
  
  run_matlab <- run_matlab_code(code, verbose=FALSE)
  
  ## get the results back from Matlab
  fa.cluster <- readMat(paste("MR_software/ImagingMR-LAS/trial_results/sub_matrix.setting_",setting,".mat", sep=""))
  
  subdataM <- as.matrix(dataM[,fa.cluster[["colSel"]][[1]]==1])
  subdataX <- as.matrix(dataX[,fa.cluster[["rowSel"]][[1]]==1])
  
  las_sel_col <- which(fa.cluster[["colSel"]][[1]]==1)
  las_sel_row <- which(fa.cluster[["rowSel"]][[1]]==1)
  print(las_sel_col)
  print(las_sel_row)
  
  return(list(subdataM, las_sel_col, ulmres.out))
}





##################################### main #############################################
run_each_setting = function(setting){
  library(matlabr)
  library(R.matlab) 
  library(MendelianRandomization) 
  library(psych) 
  library(psychTools)
  library(factoextra)
  
  mse = data.frame(stringsAsFactors = F); coef = data.frame(stringsAsFactors = F)
  
  for (idx in 1:total_replication){
    load(paste("stimulate/exposure-outcome/stimulation_4/setting_",setting,'/',idx,".exposure-outcome-iv.RData",sep = ''))
    
    result_snp = result_snp[result_snp$Pr...t.. < cutoff, ]
    if (length(unique(result_snp$exposure)) != 20) {break} 
    iv_snp = unique(result_snp$snp)
    g_matrix = g_matrix0[, iv_snp]
    
    
    g_matrix = as.data.frame(g_matrix); exposure_matrix = as.data.frame(exposure_matrix)
    names(g_matrix) = paste('snp_', 1:ncol(g_matrix), sep='')
    names(exposure_matrix) = paste('exposure_', 1:ncol(exposure_matrix), sep='')
    
    dataX = g_matrix; dataM = exposure_matrix; dataY = outcome_matrix
    
    LAS_results = save_imagingMR_LAS_results(setting, idx, dataX, dataY, dataM, p=ncol(dataX), q=ncol(dataM), cutoff)
    subdataM = LAS_results[[1]]
    las_sel_col = LAS_results[[2]]
    ulmres.out = LAS_results[[3]]
    
    beta_exposure = fa.mr(subdataM, dataX, las_sel_col, ulmres.out)
    
    ###1).mse
    mse_beta = mean((beta_exposure -  beta_xToy)^2)
    tmp.mse = data.frame(setting = setting, replicate = idx, mse_beta = mse_beta, rmse_beta = sqrt(mse_beta))
    mse = rbind(mse, tmp.mse)
    ###2).causal effect
    tmp.coef = data.frame(setting = setting, replicate = rep(idx, 20), image_name = 1:20, coef = beta_exposure)
    tmp.coef$true_beta = beta_xToy
    coef = rbind(coef, tmp.coef)
    
    
    cat('replication_',idx,'is finished','\n')
  }
  write.table(mse, paste('stimulate/results/stimulation_4/las/las.mse.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  write.table(coef, paste('stimulate/results/stimulation_4/las/las.coef.setting_',setting,'.txt', sep=''), sep="\t", row.names = F, quote = F)
  
  return(0)
}



library(doParallel)
library(foreach)

cl = makeCluster(4)
registerDoParallel(cl)
result = foreach(setting = 1:total_setting, .combine = 'rbind') %dopar% run_each_setting(setting)
stopCluster(cl)



