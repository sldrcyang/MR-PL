# MR-PL
MR-PL is a one-sample multivariable Mendelian Randomization method, which could take into consideration the correlation among numerous exposures (e.g., imaging features) followed by correcting for the winner’s curse bias.
## Prepare dependencies
Make sure you have the R package dependencies below installed and accessible in your $PATH.   

`BiocManager::install('pls')`  
`BiocManager::install('glmnet')`  
`BiocManager::install('hdi')`  
## Input files
* __g_matrix.txt__: * the genotype matrix of n samples and m snps.   
* __exposure_matrix.txt__: * the exposure matrix of n samples and k exposures.    
* __outcome_matrix.txt__: * the outcome vector of n samples.    
* __gwas_assoc.txt__: * the gwas summary statistics of m snps and k exposures, which is usually derived from a previously published study.   
Example data are listed in the `example/` folder.
## Usage
Examples describe how to perform MR-PL can be found in example.R in the `example/` folder.  
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
