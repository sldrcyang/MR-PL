rm(list = ls())
setwd("/home1/yanganyi/Desktop/MR/")
set.seed(2021)



#########################################simulate genotype matrix####################################
snp_coeff0 = read.table("true_data/genotype/UKBB/UKBB_extract_ld_snp_chr22_r20.6_sampled10000.xmat.gz", header=F, skip=2)
g_matrix = snp_coeff0[, 3:ncol(snp_coeff0)]
for (i in 1:ncol(g_matrix)){
  g_matrix[is.na(g_matrix[,i]),i] = mean(g_matrix[,i], na.rm = T)
}
g_matrix = as.matrix(g_matrix)
for (i in 1:ncol(g_matrix)){
  g_matrix[,i] = scale(g_matrix[,i], center = T, scale = T)
  if (i %% 1000 == 0 ) cat(i,'\n')
}

save(g_matrix, file=paste("stimulate/genotype/g_matrix.ukb_10000sample_chr22_r20.6_5995snp", ".RData", sep=""))

load("stimulate/genotype/g_matrix.ukb_10000sample_chr22_r20.6_5995snp.RData")




######################################simulate exposure and outcome matrix####################################
seed_vec = sample(1:1e4, 1000, replace=F) 

n_sample = nrow(g_matrix)
n_exposure = 20
n_snp = ncol(g_matrix)
n_causal_snp = 300 #5995 * 0.05
sd_error_x = 0.1
sd_error_y = 0.1


stimulation = function(setting, hg2.fomul, beta_ux_vector.fomul, beta_xToy.fomul, sd_error_x, sd_error_y){
  for (idx in 1:100) {
    print(paste(idx, "-th stimulation starts", sep=""))
    set.seed(seed_vec[idx])
    
    
    ###########################simulate effects of snps on exposures##############################
    hg2 = eval(parse(text = hg2.fomul))
    beta_g_matrix = matrix(NA, n_snp, n_exposure)
    causal_snp_idx_total = c()
    for (i in 1:n_exposure) {
      causal_snp_idx = sample(1:n_snp, size = n_causal_snp)
      causal_snp_idx_total = c(causal_snp_idx_total, causal_snp_idx)
      
      beta_causal_snp = rnorm(n_causal_snp, 0, sd = sqrt(hg2[i]/n_causal_snp))       
      beta_g_matrix[causal_snp_idx, i] = beta_causal_snp
      beta_g_matrix[-causal_snp_idx, i] = 0
    }
    causal_snp_idx_total = unique(causal_snp_idx_total)
    
    
    
    #######################simulate error matrix for each exposure 10000 * 20###############################
    error_x_matrix = matrix(NA, n_sample, n_exposure)
    for (i in 1:n_exposure) {
      error_x_matrix[,i] = rnorm(n_sample, mean=0, sd=sd_error_x)
    }
    
    
    beta_ux_vector = eval(parse(text = beta_ux_vector.fomul))

    
    u_matrix = matrix(NA, n_sample, n_exposure)
    for (i in 1:n_exposure){
      se_u = sqrt((1 - hg2[i] - sd_error_x^2)/sum(beta_ux_vector[-i]^2))
      u_matrix[,i] = rnorm(n_sample, 0, se_u)
      cat(se_u, '\n')
    }
    
    
    
    ###############################simulate exposure matrix 10000*20####################################
    exposure_matrix = matrix(NA, n_sample, n_exposure)
    for (i in 1:n_exposure){
      exposure_matrix[,i] = g_matrix %*% beta_g_matrix[,i] + u_matrix[,-i] %*% beta_ux_vector[-i] + error_x_matrix[,i]
    }
    
    
    beta_uy_vector = runif(n_exposure, 0,1)
    
    
    beta_xToy = eval(parse(text = beta_xToy.fomul))

    error_y_matrix = matrix(NA, n_sample, 1)
    error_y_matrix[,1] = rnorm(n_sample, mean=0, sd=sd_error_y)
    
    
    ##################################simulate outcome 10000 *1 ######################################
    outcome_matrix = matrix(NA, n_sample, 1)

    outcome_matrix[,1] = exposure_matrix %*% beta_xToy + u_matrix %*% beta_uy_vector + error_y_matrix[,1] 
    outcome_matrix = scale(outcome_matrix)  
    
    
    ################################generate gwas sumstat#######################
    result = data.frame(stringsAsFactors = F)
    for (i in 1:ncol(exposure_matrix)){
      for (j in 1:ncol(g_matrix)){
        lm.reg = summary(lm(exposure_matrix[,i] ~ g_matrix[,j]))
        tmp.result = data.frame(t(lm.reg$coefficients[2,]))
        tmp.result$exposure = i
        tmp.result$snp = j
        result = rbind(result, tmp.result)
      }
      cat(i, '\n')
    }
    
    result_snp = result[order(result$Pr...t..), ]
    
    
    cat(paste(idx, "-th replication finish", sep = ''), '\n')
    save(exposure_matrix, outcome_matrix, result_snp, beta_xToy, file=paste("stimulate/exposure-outcome/stimulation_4/", setting,'/', idx,".exposure-outcome-iv", ".RData", sep=""))
  }
}





##################################### generate setting summary ##################################
hg2_list = list('rep(0.2, n_exposure)', 'rep(0.3, n_exposure)', 
                'rep(0.4, n_exposure)', 'rep(0.5, n_exposure)') 
beta_ux_vector_list = list( 'runif(n_exposure, 0, 5)', 'runif(n_exposure, 5, 10)' ) 
beta_xToy_list = list('sample(c(-0.3,0,0.3), size=n_exposure, replace=T)',
                      'sample(c(-0.2,0,0.2), size=n_exposure, replace=T)', 
                      'sample(c(-0.1,0,0.1), size=n_exposure, replace=T)', 
                      'runif(n_exposure, -0.3, 0.3)',
                      'runif(n_exposure, -0.2, 0.2)', 
                      'runif(n_exposure, -0.1, 0.1)')

setting.summary = data.frame(stringsAsFactors = F)
count = 0
for (hg2 in hg2_list){
  for (beta_ux_vector in beta_ux_vector_list){
    for (beta_xToy in beta_xToy_list){
      count = count + 1
      setting = paste('setting_', count, sep='')
      
      setting.tmp = data.frame(setting=setting, hg2=hg2, beta_ux_vector=beta_ux_vector, beta_xToy=beta_xToy, stringsAsFactors=F)
      setting.summary = rbind(setting.summary, setting.tmp)
    }
  }
}



############################################# Main ###########################################
run_each_setting = function(id_row){
  setting = setting.summary[id_row, 'setting']
  hg2.fomul = setting.summary[id_row, 'hg2']
  beta_ux_vector.fomul = setting.summary[id_row, 'beta_ux_vector']
  beta_xToy.fomul = setting.summary[id_row, 'beta_xToy']
  
  stimulation(setting, hg2.fomul, beta_ux_vector.fomul, beta_xToy.fomul, sd_error_x=0.1, sd_error_y=0.1)
  return(0)
}


library(doParallel)
library(foreach)

cl = makeCluster(20) 
registerDoParallel(cl)
result = foreach(id_row = 1:nrow(setting.summary), .combine = 'rbind') %dopar% run_each_setting(id_row)
stopCluster(cl)

