#-------------------- outcome ~ pls.predict---------------------------
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)
library(glmnet)
library(hdi)


load("true_data/results/ukb.c=20.snp_matrix_r2_0.1.and.wm_mean_FA_pls_pred.RData")
outcome_matrix = fread('true_data/phenotype/UKBB/ukb.394_pheno_matrix.processed.allsamples.txt', data.table = F, header = T, sep = "\t", stringsAsFactors = F)
pheno_info = fread('true_data/phenotype/UKBB/ukb.394_pheno_info_nsample.select.txt', data.table = F, header = T, stringsAsFactors = F)

names(outcome_matrix)[1] = 'IID'
sample_id = Reduce(intersect, list(pls.pred$IID, outcome_matrix$IID))
pls.pred = pls.pred[match(sample_id, pls.pred$IID), ]
outcome_matrix = outcome_matrix[match(sample_id, outcome_matrix$IID), ]


set.seed(1)
grid = 10 ^ seq(10, -2, length=100)
result = data.frame(stringsAsFactors = F)
for (i in 2:ncol(outcome_matrix)) {
  fieldID = names(outcome_matrix)[i]
  valuetype = pheno_info[which(pheno_info$FieldID == fieldID), 'ValueType.after_mapping']
  
  data = cbind(outcome_matrix[,i], pls.pred[,-1])
  names(data)[1] = 'outcome'
  row.names(data) = pls.pred[,1]
  data = subset(data, data$outcome != 'NA')
  if (nrow(data) < 2000) next
  
  
  ###### LASSO ######
  if (valuetype == "Integer" | valuetype == "Continuous"){
    data$outcome = scale(data$outcome)
    
    family = "gaussian"
    type = "link"
  } else if(length(unique(data$outcome)) == 2){
    if (!( all(unique(data$outcome) %in% c(0,1)) )){
      data$outcome = gsub(unique(data$outcome)[1], 0, data$outcome)
      data$outcome = gsub(unique(data$outcome)[2], 1, data$outcome)
      data$outcome = as.numeric(data$outcome)
    }

    if(sum(data$outcome) / nrow(data) <= 0.001 | sum(data$outcome) / nrow(data) >= 0.999)  next

    valuetype = 'Binary Category'  
    
    family = "binomial"
    type = "response"
  } else if(length(unique(data$outcome)) > 2){
    next
  }
  
  data = as.matrix(data)
  
  reg.mod = glmnet(data[,-1], data[,1], alpha = 1, lambda = grid, family = family)
  cv.out = cv.glmnet(data[,-1], data[,1], alpha=  1, family = family)
  bestlam = cv.out$lambda.min
  cat(bestlam, '\n')
  
  reg.coef = predict(reg.mod, type = 'coefficients', s = bestlam)
  coef_name = reg.coef@Dimnames[[1]] [reg.coef@i +1]
  coef_val = reg.coef@x
  if ("(Intercept)" %in% coef_name){
    id = which(coef_name == "(Intercept)")
    coef_name = coef_name[-id]
    coef_val = coef_val[-id]
  }
  
  if (length(coef_val) == 0) {next}
  
  #desparsified lasso---
  predictions = predict(reg.mod, newx = data[,-1], s = bestlam)
  residuals = data[,1] - predictions
  sd = sd(residuals)
  beta = rep(0, ncol(data)-1)
  sorting = match(coef_name, colnames(data)[-1])
  beta[sorting] = coef_val
  
  outlasso = lasso.proj(x = data[,-1], y = data[,1], betainit=beta, sigma = sd) 
  p = outlasso$pval[coef_name]
  #----------------------------------------------------------------
  
  
  ############################ Sargan test ########################
  outcome.pred = predict(reg.mod, newx = data[,-1], s = bestlam)
  residual = data.frame(data[,'outcome'] - outcome.pred)
  residual$IID = row.names(data)
  sample_id = Reduce(intersect, list(residual$IID, snp_matrix$IID))
  residual = residual[match(sample_id, residual$IID),]
  snp_matrix1 =  snp_matrix[match(sample_id,  snp_matrix$IID),]
  
  reg.sargan = lm(residual[,1] ~ as.matrix(snp_matrix1[,-1]))
  R2 = summary(reg.sargan)$r.squared
  sargan.stat = nrow(snp_matrix1) * R2
  df = ncol(snp_matrix1) - ncol(pls.pred)
  sargan.p = 1 - pchisq(sargan.stat, df = df) 
  if (sargan.stat < qchisq(p = 0.95, df = df)) {
    sargan.test = 'pass'
  } else {sargan.test = 'reject'}
  
  cat(sargan.test, '\n')
  #################################################################
  
  
  tmp.result = data.frame(phenotype = fieldID, category = pheno_info[which(pheno_info$FieldID == fieldID), 'Category'], 
                          field = pheno_info[which(pheno_info$FieldID == fieldID), 'Field'], valuetype = valuetype,
                          image_name = coef_name, coef = coef_val, p_lasso = p, sargan.test = sargan.test, 
                          sargan.p = sargan.p, nsample = nrow(data))
  result = rbind(result, tmp.result)
  cat('The', i, 'th outcome is OK.', '\n')
}

result_output = subset(result, result$coef != 0)

#transform brain region id to brain region name
field_info = fread('true_data/phenotype/UKBB/Documentation/Data_Dictionary_Showcase.csv',
                   data.table = F, header = T, stringsAsFactors = F)
field_info = field_info[,c('FieldID', 'Field')]
result_output = merge(result_output, field_info, by.x = 'image_name', by.y = 'FieldID', all.x = T)
result_output = result_output[,c("phenotype","category","field","image_name","Field",
                                 "valuetype","coef","p_lasso","sargan.test","sargan.p","nsample")]
result_output$fdr = p.adjust(result_output$p_lasso, method = 'fdr', n = nrow(result_output))


write.table(result_output,'true_data/results/ukb.c=20.snp_r2_0.1.wm_mean_FA.370+.pls+lasso.coef_p.txt',sep="\t",row.names=F,quote=F)
