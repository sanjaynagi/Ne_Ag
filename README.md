## Ag1000g - effective population size based on linkage disequilibrium

Analysis of Ag1000g populations looking at effective pop sizes.
Initially, the idea was just to get used to using LDNe to estimate changes in effective population size, thinking about future genomic data from the LLINEUP trial.

However, once the snakemake pipeline was written, it was not much more effort to apply this to all phase 2 populations, and then phase 3, this time split by collection location.  

the repo also contains some old code for running IBDNe. However, it is difficult to quantify short tracts of IBD as would be found in most mosquito populations, and thus I stopped this work.
#### input

- Ag1000g Zarrs 

## Code 

- snakefile
- analysis/scripts/Zarr_to_LDNe.py

## Output

- Effective pop sizes .png
- analysis/Ne_manifest.tsv - contains all results for chrom 3L, other chromosomes still need parsing with the jupyter notebook 'Ne_Ag_LDNeOutput.ipynb' 



## Negative NE estimates from LDNe 
Note - Why do you get Negative Ne estimates sometimes?
(Waples,Do , 2010) evol appl. 

As shown in eqn (2a), before estimating Ne in the LD method, the expected contribution of sampling error is subtracted from the empirical XX 
 If Ne is large, or if only limited data are available, by chance mean can be smaller than the sample size correction, in which case the estimate of Ne will be negative. A related phenomenon can occur with the standard temporal method (Nei and Tajima 1981; Waples 1989) and with unbiased estimators of genetic differentiation (Nei 1978; Weir and Cockerham 1984). Negative estimates occur when the genetic results can be explained entirely by sampling error without invoking any genetic drift, so the biological interpretation is An external file that holds a picture, illustration, etc. Object name is eva0003-0244-mu75.jpg = ∞ (Laurie-Ahlberg and Weir 1979; Nei and Tajima 1981). In this situation, the user can conclude that the data provide no evidence that the population is not ‘very large’. However, even if the point estimate is negative, if adequate data are available the lower bound of the CI generally will be finite and can provide useful information about plausible limits Ne. 
