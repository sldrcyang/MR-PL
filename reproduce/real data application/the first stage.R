rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)


##read image matrix---
image_matrix = fread('true_data/phenotype/UKBB/ukb.image_matrix.wm_mean_FA.txt', data.table = F, header = T, sep = "\t",stringsAsFactors = F)
#image excluded
image_exclude = paste(c(25056,25057,25064:25071,25102,25103), '-2.0', sep = '')
image_matrix = image_matrix[,-which(names(image_matrix) %in% image_exclude)]
names(image_matrix)[1] = 'IID'
for (i in 2:ncol(image_matrix)) {
  image_matrix[,i] = as.numeric(image_matrix[,i])
  image_matrix[is.na(image_matrix[,i]),i] = mean(image_matrix[,i], na.rm=T)
}


############################## regress out the influence of covariates ###############################
##read covariates
cov_processed = fread("true_data/phenotype/UKBB/ukb.cov_matrix.txt", data.table = F, header = T, stringsAsFactors = F)
names(cov_processed)[1] = 'IID'

sample_id = Reduce(intersect, list(image_matrix$IID, cov_processed$IID))
image_matrix = image_matrix[match(sample_id, image_matrix$IID),]
cov_processed = cov_processed[match(sample_id, cov_processed$IID),]

image_matrix_resid = as.data.frame(matrix(nrow = nrow(image_matrix), ncol = ncol(image_matrix)))
names(image_matrix_resid) = names(image_matrix)
image_matrix_resid$IID = image_matrix$IID

for (i in 2:ncol(image_matrix)) {
  data = data.frame(image_matrix[,i], cov_processed[,-1])
  names(data)[1] = c('phenotype')
  
  formul = paste("phenotype ~ ",names(cov_processed)[2],sep='')
  for (k in 3:ncol(data)){
    formul = paste(formul, "+", names(data)[k], sep="")
  }
  reg = lm(formul, data = data)
  #y = predict(reg, data)
  residual = residuals(reg)
  image_matrix_resid[,i] = residual
  
  cat(i, '\n')
}

for (i in 2:ncol(image_matrix_resid)) {
  image_matrix_resid[,i] = as.numeric(image_matrix_resid[,i])
  image_matrix_resid[,i] = scale(image_matrix_resid[,i], center = T, scale = T) 
}



##3)read genotype matrix ---
snp_matrix = read.table("true_data/genotype/UKBB/UKBB_extract_wm_meanFA.r2_0.1.xmat.gz", header=T, sep = " ", stringsAsFactors = F)
snp_matrix = snp_matrix[-1,]
snp_matrix = snp_matrix[,-1]
snp_matrix = snp_matrix[,-ncol(snp_matrix)]

####extract snps needed###
snp_info = read.table('true_data/selected_snps/zhuhongtu.signif_snps_FA.r2_0.1.txt', header=T, sep = "\t", stringsAsFactors = F)
snp_info = snp_info[snp_info$ID != 'IFO', ]
length(unique(snp_info$rsID))

image_ID_before = unique(snp_info$ID)

#winner's curse correction---------------------------------------------------------------------------------------
snp_info$abs.t = qnorm(1 - snp_info$p/2)
c = 20
cutoff = 5e-08
snp_info$abs.lasso.t = snp_info$abs.t - c / (-log10(cutoff) + log10(nrow(snp_matrix)))
snp_info[snp_info$abs.lasso.t < 0, 'abs.lasso.t'] = 0
snp_info$corrected.p = 2 * (1 - pnorm(snp_info$abs.lasso.t))
snp_info = snp_info[snp_info$corrected.p < cutoff, ]

snp_matrix = snp_matrix[, which(names(snp_matrix) %in% c('IID', snp_info$rsID))]
for (i in 2:ncol(snp_matrix)){
  snp_matrix[,i] = as.numeric(snp_matrix[,i])
  snp_matrix[is.na(snp_matrix[,i]),i] = mean(snp_matrix[,i], na.rm=T) 
}
ncol(snp_matrix)  

image_ID_after = unique(snp_info$ID)
setdiff(image_ID_before, image_ID_after)
image_ID_delete = c('25062-2.0', '25063-2.0')
image_matrix_resid = image_matrix_resid[, -which(names(image_matrix_resid) %in% image_ID_delete)]



################################## PLS prediction ###################################
library(pls)
library(parallel)

pred_pls = function(plsdata1_Z, plsdata2_Z, ncomp){
  pls.options(parallel = makeCluster(6, type = "PSOCK"))
  
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = ncomp, validation="CV", jackknife=TRUE)
  
  RMSE = RMSEP(pls.fit) 
  best_ncomp = which.min(apply(RMSE$val [2,1:ncomp,-1], 2, mean))
  cat(best_ncomp, '\n')
  
  pls.fit = plsr(plsdata2_Z ~ plsdata1_Z, ncomp = best_ncomp, validation = "none") 
  coef = as.matrix(as.data.frame(coef(pls.fit)))
  #jack_result = jack.test(pls.fit)
  
  stopCluster(pls.options()$parallel)  # The cluster should be stopped manually
  
  pls.pred = as.data.frame((predict(pls.fit, plsdata1_Z, ncomp = best_ncomp)))
  return(pls.pred)
}

sample_id = Reduce(intersect, list(snp_matrix$IID, image_matrix_resid$IID))
snp_matrix = snp_matrix[match(sample_id, snp_matrix$IID),]
image_matrix_resid = image_matrix_resid[match(sample_id, image_matrix_resid$IID),]


set.seed(1)
pls.pred = pred_pls(plsdata1_Z=scale(snp_matrix[,-1]), plsdata2_Z=scale(image_matrix_resid[,-1]), ncomp=ncol(image_matrix_resid[,-1]))   
names(pls.pred) = gsub(paste('-', unlist(strsplit(names(pls.pred)[1], "-"))[2], sep=''), '', names(pls.pred))
pls.pred$IID = snp_matrix$IID
pls.pred = pls.pred[, c(ncol(pls.pred), 1:(ncol(pls.pred)-1))]

save(snp_matrix, pls.pred, file = paste('true_data/results/ukb.c=',c,'.snp_matrix_r2_0.1.and.wm_mean_FA_pls_pred.RData', sep=''))



