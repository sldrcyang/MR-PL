###transform BGEN to plink format and extract snps and samples
#chr1
/share/inspurStorage/home1/yanganyi/Desktop/Software/plink2 \
--bgen /share/inspurStorage/home1/UKB_Gene_v3/disk10t_2/batch/imp/ukb_imp_chr1_v3.bgen ref-first \
--sample /share/inspurStorage/home1/UKB_Gene_v3/disk10t_2/batch/sample/ukb19542_imp_chr1_v3_s487378.sample \
--extract /home1/yanganyi/Desktop/MR/select_snps/ukb_wm_meanFA.txt \
--keep /home1/yanganyi/Desktop/MR/phenotype/UKBB/ukb.sample_keep.wm_meanFA.txt \
--make-bed \
--out /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1


############################## merge plink file of 22 chromosomes #############################
for i in $(seq 2 21);do j="$(expr $i + 1)";/home1/yanganyi/Desktop/Software/plink1/plink \
--bfile /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-$i \
--bmerge /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr$j \
--make-bed \
--out /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-$j;\
rm /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-$i.bed;\
rm /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-$i.bim;\
rm /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-$i.fam;\
rm /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-$i.log;\
done



############################# extract readable matrix using GCTA ########################
gcta64 --bfile /home1/yanganyi/Desktop/MR/genotype/UKBB/extract_wm_meanFA/ukb.extract.wm_meanFA.chr1-22 \
--recode --out /home1/yanganyi/Desktop/MR/genotype/UKBB/UKBB_extract_wm_meanFA


