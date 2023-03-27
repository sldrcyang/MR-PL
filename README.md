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
Example data are listed in the example folder.
## Usage
There is an example describe how to perform MR-PL in example.R.
