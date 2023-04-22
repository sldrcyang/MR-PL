############################ make matched ukb_FieldID-file_name pairs#################################
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")

path = "phenotype/UKBB/field_csv/"
file = list.files(path) 

dir = paste(path, file, sep = "")
n = length(dir) 

result = data.frame(stringsAsFactors = F)
for (i in 1:n){
  field_info = read.table(file = dir[i], sep = ',', nrows = 1, stringsAsFactors = F)
  field_info = field_info[1, -which(field_info[1,] == 'eid')]
  
  tmp.result = as.data.frame(t(field_info[1,]))
  tmp.result$csv_file = file[i]
  
  result = rbind(result, tmp.result)
}
names(result) = c('field_id', 'csv_file')
write.table(result, 'phenotype/UKBB_phenotype/ukb.field_id-csv_file.txt', sep="\t", row.names = F, quote = F)



##################################### extract sample ID and sample QC ###########################################
# british (exclude heterogeneity, missing rate, relatedness exclusions, samples with inconsistent genetic sex and sex)
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)

file = fread('/share/inspurStorage/home2/UKB_Preprocessed_Data/Behavior/ukb26427.csv', 
             data.table = F, sep = ',', header = T, stringsAsFactors = F)

#not British id 21000
file_subset = file[,c('eid', '21000-0.0')]
ex_id1 = subset(file_subset, file_subset$`21000-0.0` != '1001') $ `eid`  #not british
length(ex_id1)  #59028

# Recommended genomic analysis exclusions (poor heterozygosity/missingness) 22010
file_subset = file[,c('eid', '22010-0.0')]
ex_id2 = subset(file_subset, file_subset$`22010-0.0` == '1') $ `eid`
length(ex_id2)  #480

# retain the first pair in a Genetic relatedness pair
file_subset = file[,c('eid', '22011-0.0')]
file_subset = subset(file_subset, file_subset$`22011-0.0` != 'NA')
file_subset1 = file_subset[duplicated(file_subset$`22011-0.0`),]
ex_id3 = file_subset1$eid
length(ex_id3)  #8141

#  exclude samples with inconsistent genetic sex and sex
file_subset = file[,c('eid', '31-0.0', '22001-0.0')]
ex_id4 = c()
for (i in 1:nrow(file_subset)){
  if (!is.na(file_subset[i,2]) & !is.na(file_subset[i,3])){
    if (file_subset[i,2] != file_subset[i,3]){
      ex_id4 = c(ex_id4, file_subset$eid[i])
    }
  }
}
length(ex_id4)  #378

ex_id = Reduce(union, list(ex_id1, ex_id2, ex_id3, ex_id4))
length(ex_id)  #67565
write.table(ex_id, 'phenotype/UKBB/ukb.sample_exclude.txt', sep="\t", row.names = F, col.names = F, quote = F)



############################################ extract covariates ####################################################
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)

file1 = fread('/share/inspurStorage/home2/UKB_Preprocessed_Data/Behavior/ukb44306.csv', data.table = F, sep = ',', header = T, stringsAsFactors = F)
file2 = fread('/share/inspurStorage/home2/UKB_Preprocessed_Data/Behavior/ukb47737.csv', data.table = F, sep = ',', header = T, stringsAsFactors = F)

col_name1 =  c('eid', '31-0.0','54-0.0','21001-0.0','21003-0.0')
col_name2 = c('eid', '22000-0.0', paste('22009-0.', 1:20, sep=''))  #含有时间点的列
cov1 = file1[, which(names(file1) %in% col_name1)]
cov2 = file2[, which(names(file2) %in% col_name2)]

cov = merge(cov1, cov2, by='eid')

names(cov)
names(cov) = c('eid', 'sex', 'site', 'bmi', 'age', 'batch', paste('pca',1:20,sep = ''))

ex_id = read.table('phenotype/UKBB/ukb.sample_exclude.txt', sep="\t", header = F, stringsAsFactors = F)
cov = subset(cov, !(cov$eid %in% ex_id$V1)) 


cov$sex = gsub('0', 'F', cov$sex)
cov$sex = gsub('1', 'M', cov$sex)

unique_id = unique(cov$site)

unique_id = unique_id[!is.na(unique_id)]
for (i in 1:length(unique_id)){
  cov$site = gsub(unique_id[i], paste('site', i, sep='_'), cov$site)
}

cov$batch = as.character(cov$batch)
unique_id = unique(cov$batch)

unique_id = unique_id[!is.na(unique_id)]
unique_id

for (i in 1:nrow(cov)){
  idx = match(cov[i,'batch'], unique_id)
  if (length(idx) != 0){
    cov[i,'batch'] =  paste('batch', idx, sep='_')
  }
  cat(i, '\n')
}

#age^2
cov$age_2 = cov$age ^ 2
cov = cov[, c(1:3,6,4,5,ncol(cov), 7:(ncol(cov)-1))]

for (i in c(5:ncol(cov))){
  cov[is.na(cov[,i]), i] = mean(cov[,i], na.rm=T)
}
write.table(cov, 'phenotype/UKBB/ukb.cov_matrix.txt', sep="\t", row.names = F, col.names = T, quote = F)



####################### extract brain imaging matrix (white matter mean FA) #########################
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)

file = fread('/share/inspurStorage/home2/UKB_Preprocessed_Data/Behavior/ukb40844.csv', 
             data.table = F, sep = ',', header = T, stringsAsFactors = F)

col_name =  c('eid', paste(c(25056:25103), '2.0', sep = '-'))  
image = file[, which(names(file) %in% col_name)]
image = image[!is.na(image[,2]),]
nrow(image)   #33512 
#exclude samples not pass QC
ex_id = read.table('phenotype/UKBB/ukb.sample_exclude.txt', sep="\t", header = F, stringsAsFactors = F)
image = subset(image, !(image$eid %in% ex_id$V1))

for (i in c(2:ncol(image))){
  image[is.na(image[,i]), i] = mean(image[,i], na.rm=T)
}
write.table(image, 'phenotype/UKBB/ukb.image_matrix.wm_mean_FA.txt', sep="\t", row.names = F, col.names = T, quote = F)

sample_extract = data.frame('FID'=image$eid, 'IID'=image$eid)
write.table(sample_extract, 'phenotype/UKBB/ukb.sample_keep.wm_meanFA.txt', sep="\t", row.names = F, col.names = F, quote = F)
