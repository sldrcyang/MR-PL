# MR-PL
MR-PL is a one-sample multivariable Mendelian Randomization method, which could take into consideration the correlation among numerous exposures (e.g., imaging features) followed by correcting for the winner’s curse bias.
## Prepare dependencies
Make sure you have the R package dependencies below installed and accessible in your $PATH.   

`BiocManager::install('pls')`  
`BiocManager::install('glmnet')`  
`BiocManager::install('hdi')`  
## Input files
* __g_matrix.txt__:  the genotype matrix of n samples and m snps.   
* __exposure_matrix.txt__:  the exposure matrix of n samples and k exposures.    
* __outcome_matrix.txt__:  the outcome vector of n samples.    
* __gwas_assoc.txt__:  the gwas summary statistics of m snps and k exposures, which is usually derived from a previously published study.   

Example data are listed in the `example/` folder.
## Usage
Examples describe in detail how to perform MR-PL can be found in example.R in the `example/` folder.  
   
Briefly, if the gwas summary statistics used to selected instrumental snps is based on the same dataset used for MR analyses, the following can be run
```
mr_result = mr_pl(g_matrix0, 
                  exposure_matrix0, 
                  outcome, 
                  gwas_assoc, 
                  cutoff = 5e-8, 
                  wcc = TRUE, 
                  c = 20, 
                  pleiotropy_test = TRUE)
```
```
if the gwas summary statistics used to selected IVs is based on another independent dataset, the following can be run  
mr_result = mr_pl(g_matrix0, 
                  exposure_matrix0, 
                  outcome, 
                  gwas_assoc, 
                  cutoff = 5e-8, 
                  wcc = FALSE, 
                  c = 0, 
                  pleiotropy_test = TRUE)
```    
## Output
mr_result$main_results: exposures with non-zero causal estimate and its corresponding p-value from lasso projection method
```    
    exposure_name causal_estimate lasso_proj_p
1   exposure_1      -0.1902849 5.204576e-19
2  exposure_10      -0.1676381 3.987234e-14
3  exposure_11      -0.2026717 1.558848e-22
4  exposure_12      -0.2092596 2.271563e-31
5  exposure_14       0.1878917 4.603154e-19
6  exposure_18      -0.1930250 3.658456e-18
7   exposure_2       0.1902596 2.348510e-20
8   exposure_4      -0.1685691 1.807560e-15
9   exposure_6      -0.1814735 4.848565e-19
10  exposure_9      -0.1850430 3.273956e-16
```    
mr_result$pleiotropy_test.p: p value of pleiotropy test; if p < 0.05, there exists horizontal pleiotropy, then this MR result should be discarded.

mr_result$iv_include: the instrumental snps used for MR analysis (after winner's curse correction)

mr_result$exposure_include: the exposures used for MR analysis (after winner's curse correction)




