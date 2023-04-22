############## make phenotype matrix #############
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)

pheno_info = fread('true_data/phenotype/UKBB/ukb.394_pheno_info_nsample.select.txt', data.table = F, header = T, stringsAsFactors = F)
pheno_info = pheno_info[order(pheno_info$csv_file), ]

pheno_matrix = data.frame(stringsAsFactors = F)
for (csv_id in 1:length(unique(pheno_info$csv_file))){
  csv = unique(pheno_info$csv_file)[csv_id]
  csv_pheno = subset(pheno_info, pheno_info$csv_file == csv)
  
  file = fread(paste('/share/inspurStorage/home2/UKB_Preprocessed_Data/Behavior/', csv, sep=''), data.table = F, sep = ',', header = T, stringsAsFactors = F)
  file_subset = file[, c('eid', csv_pheno$FieldID)]
  
  if(csv_id == 1) {
    pheno_matrix = file_subset
  } else{
    pheno_matrix = merge(pheno_matrix, file_subset, by = 'eid', all.y = T)
  }
  cat(nrow(csv_pheno), '\n')
  cat('csv_file:', csv, 'is OK.', '\n')
}
write.table(pheno_matrix, 'true_data/phenotype/UKBB/ukb.394_pheno_matrix.allsamples.txt', sep="\t", row.names = F, col.names = T, quote = F)



#####Replace missing value coding with missing values and mapping categorical multiple to continuous variables###
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR")
library(data.table)

pheno_matrix0 = fread('true_data/phenotype/UKBB/ukb.394_pheno_matrix.allsamples.txt', data.table = F, header = T, stringsAsFactors = F)
pheno_info = fread('true_data/phenotype/UKBB/ukb.394_pheno_info_nsample.select.txt', data.table = F, header = T, stringsAsFactors = F)
#missing value coding
coding_info = read.table('true_data/phenotype/UKBB/ukb.coding-missing.txt', sep = '\t', header = T, stringsAsFactors = F)

for (i in 1:nrow(pheno_info)){
  fieldID = pheno_info[i, 'FieldID']
  coding = pheno_info[i, 'Coding']
  if (is.na(coding)) next
  
  missing_val = coding_info[which(coding_info$coding == coding), 'missing']
  if (length(missing_val) == 0) next
  missing_val = as.numeric(unlist(strsplit(missing_val, split = ",")))

  for (val in missing_val){
    pheno_matrix0[, fieldID] = gsub(val, NA, pheno_matrix0[, fieldID])  
    if (pheno_info[i, 'ValueType'] == "Continuous" | pheno_info[i, 'ValueType'] == "Integer"){
      pheno_matrix0[, fieldID] = as.numeric(pheno_matrix0[, fieldID])
    }
  }
  cat(i, 'is OK.\n')
}


pheno_matrix = pheno_matrix0

pheno_info$ValueType.after_mapping = pheno_info$ValueType
mapping_info = read.table('true_data/phenotype/UKBB/ukb.394_pheno.mapping.txt', sep = '\t', header = T, stringsAsFactors = F)
mapping_info = subset(mapping_info, mapping_info$mapping != "" )
mapping_info$mapping = gsub('[{]', '', mapping_info$mapping)
mapping_info$mapping = gsub('[}]', '', mapping_info$mapping)

for (i in 1:nrow(mapping_info)){
  field = mapping_info[i, 'column']
  fieldID = names(pheno_matrix)[grep(field, names(pheno_matrix))]
  cat(unique(pheno_matrix[, fieldID]), '\n')
  
  mapping0 = mapping_info[i, 'mapping']
  mapping0 = unlist(strsplit(mapping0, ', '))
  mapping0 = strsplit(mapping0, ': ')
  mapping = unlist(lapply(mapping0, function(x) x[2]))
  names(mapping) = unlist(lapply(mapping0, function(x) x[1]))  
  mapping
  
  nrow_selected = c()
  for(n in 1:length(mapping)){
    nrow = which(pheno_matrix[, fieldID] == as.numeric(names(mapping)[n]))
    nrow = setdiff(nrow, nrow_selected) 
    pheno_matrix[nrow, fieldID] = mapping[n]  
    nrow_selected = c(nrow_selected, nrow)
  }
  pheno_matrix[, fieldID] = as.numeric(pheno_matrix[, fieldID])
  pheno_info[which(pheno_info$FieldID == fieldID), 'ValueType.after_mapping'] = 'Continuous'
  cat(unique(pheno_matrix[, fieldID]), '\n')
  
  cat(i, '\n')
}

Categ_multi = subset(pheno_info, pheno_info$ValueType == 'Categorical multiple')
Categ_multiID = Categ_multi$FieldID
pheno_matrix = pheno_matrix[, -which(names(pheno_matrix) %in% Categ_multiID)]


write.table(pheno_info, 'true_data/phenotype/UKBB/ukb.394_pheno_info_nsample.select.txt', sep="\t", row.names = F, col.names = T, quote = F)
write.table(pheno_matrix, 'true_data/phenotype/UKBB/ukb.394_pheno_matrix.processed.allsamples.txt', 
            sep="\t", row.names = F, col.names = T, quote = F)
