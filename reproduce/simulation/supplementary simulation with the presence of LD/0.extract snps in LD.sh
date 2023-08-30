#perform quality control on UKB's chr22 file
/home1/yanganyi/Desktop/Software/plink1/plink \
--bfile /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2 \
--maf 0.01 --hwe 1e-6 --geno 0.2 --mind 0.2 \
--make-bed \
--out /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc1

#7249 snps
wc -l /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc1.bim

#pruing
/home1/yanganyi/Desktop/Software/plink1/plink --bfile /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc1 \
--indep-pairwise 50 5 0.6 \
--out /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/indepSNP_pruned
#5995 snps
wc -l /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/indepSNP_pruned.prune.in

#文件中extract snps from .prune.in file
/home1/yanganyi/Desktop/Software/plink1/plink --bfile /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc1 \
--extract /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/indepSNP_pruned.prune.in \
--make-bed \
--out /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc2


#Remove samples having brain imaging data from previously quality controlled data, and randomly select 10,000 samples for simulation
rm(list=ls())
setwd("/home1/yanganyi/Desktop/MR/")

#read eid of all samples
all_sample = read.table('true_data/genotype/UKBB/ukb_chr22_v2_qc2.fam', sep = ' ', header = F, stringsAsFactors = F)
all_sample = all_sample$V1

#samples not passing QC 
sample_exclude = read.table('true_data/phenotype/UKBB/ukb.sample_exclude.txt', header = F, stringsAsFactors = F)
sample_exclude = sample_exclude$V1

#samples having brain imaging data
image_sample = read.table('true_data/phenotype/UKBB/ukb.sample_keep.wm_meanFA.txt', sep = ' ', header = F, stringsAsFactors = F)
image_sample = image_sample$V1

#remove samples not passing QC and having brain imaging data
sample_extract = all_sample[which(!(all_sample%in%sample_exclude) & !(all_sample%in%image_sample))]

#randomly select 10,000 samples
set.seed(2022)
n_sampled = sample(sample_extract, 10000)
sample_sampled = data.frame('FID'=n_sampled, 'IID'=n_sampled)
write.table(sample_sampled, 'true_data/genotype/UKBB/ukb.sampled10000.txt', sep="\t", row.names = F, col.names = F, quote = F)


 
#Extract the genotype data of the 10000 samples on chr22 using GCTA
/home1/yanganyi/Desktop/Software/plink1/plink \
--bfile /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc2 \
--keep /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb.sampled10000.txt \
--make-bed \
--out /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc3_sampled10000


#convert binary files to readable file format using GCTA
/share/inspurStorage/home1/yanganyi/Desktop/Software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
--bfile /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/extract_ld_snp/ukb_chr22_v2_qc3_sampled10000 \
--recode \
--out /share/inspurStorage/home1/yanganyi/Desktop/MR/true_data/genotype/UKBB/UKBB_extract_ld_snp_chr22_r20.6_sampled10000


