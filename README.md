# MR-PL
MR-PL is a one-sample multivariable Mendelian Randomization method, which could take into consideration the correlation among numerous exposures (e.g., brain imaging features) followed by correcting for the winner’s curse bias.
## Prepare dependencies
Make sure you have the R package dependencies below installed and accessible in your $PATH.   

`BiocManager::install('pls')`  
`BiocManager::install('glmnet')`  
`BiocManager::install('hdi')`  
## Input files
* __g_matrix.txt__:  the genotype matrix of *n* samples and *m* snps.   
* __exposure_matrix.txt__:  the exposure matrix of *n* samples and *k* exposures.    
* __outcome_matrix.txt__:  the outcome vector of *n* samples.    
* __gwas_assoc.txt__:  the gwas summary statistics of *m* snps and *k* exposures, which is usually derived from a previously published study.  
 
Where the genotype matrix, exposure matrix and outcome are from the same dataset.  

Example data are listed in the `example/` folder.
## Usage
Detail information describes how to perform MR-PL can be found in example.R in the `example/` folder.  
   
Briefly, if the gwas summary statistics used to selected instrumental snps is based on the same dataset used for MR analyses, the following can be run
```
mr_result = mr_pl(g_matrix0, 
                  exposure_matrix0, 
                  outcome, 
                  gwas_assoc, 
                  cutoff = 5e-8, 
                  wcc = TRUE, 
                  c = 20,  #c = 15/25
                  pleiotropy_test = TRUE)
```
If the gwas summary statistics used to selected IVs is based on another independent dataset, the following can be run  
```
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
mr_result$main_results: exposures with non-zero causal estimate and its corresponding p-value from lasso projection method.
```  
   exposure_name causal_estimate lasso_proj_p
   exposure_1      -0.1908863 1.319257e-15
   exposure_10     -0.1690168 1.085111e-12
   exposure_11     -0.2010710 2.174622e-19
   exposure_12     -0.2086834 5.944588e-30
   exposure_14      0.1851574 3.035465e-18
   exposure_18     -0.1917043 1.215523e-14
   exposure_2       0.1868070 9.948648e-18
   exposure_4      -0.1706325 1.148394e-13
   exposure_6      -0.1796765 1.577770e-16
   exposure_9      -0.1832399 5.285689e-13
```    
mr_result$pleiotropy_test.p: p value of pleiotropy test; if p < 0.05, there exists horizontal pleiotropy, then this MR result should be discarded.
```    
1
```    
mr_result$iv_include: the instrumental snps used for MR analysis (after winner's curse correction).
```    
"snp_256" "snp_208" "snp_795" "snp_76"  "snp_962" ...
```    
mr_result$exposure_include: the exposures used for MR analysis (after winner's curse correction).
```    
"exposure_1"  "exposure_10" "exposure_11" ...
```
mr_result$R2: the prediction R2 of the outcome
```
0.304
```
## Reproduce
The the `reproduce/` folder contains all the codes to reproduce our simulation results, including the baseline simulation and supplementary simulation.
