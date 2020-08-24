## Ag1000g - effective population size based on linkage disequilibrium

This repo contains analyses and code for LDNe estimates of Ag1000g phase 3 populations.

Initially, the idea was just to get used to using LDNe to estimate changes in effective population size, thinking about future genomic data from the LLINEUP trial.

<<<<<<< HEAD
However, once the snakemake pipeline was written, it was not much more effort to apply this to all phase 2 populations, and then subsequently phase 3, this time split by specific collection location, rather than the phase 2 'populations'. I used a minimum of 9 samples per location/species/year combination. Really I think a sample size of at least ~15 is necessary, but in phase 3 there was a couple of Kenyan populations of n=9 which I wanted to try and include.

I thought it would be nice to try and plot the output in an interactive map with Bokeh and Geopandas.

At some point, I may run the pipeline varying the number of input SNPs, and input sample size, as well as trying out the Temporal estimator on the Burkina 2012/14 collections. Any thoughts and feedback is welcome :) 
=======
However, once the snakemake pipeline was written, it was not much more effort to apply this to all phase 2 populations, and then phase 3, this time split by collection location, with a minimum of 9 samples per location/species/year combination. Really I think a sample size of at least ~15 is necessary, but in phase 3  there was a couple of Kenyan populations of n=9 which I wanted to include. 
>>>>>>> 35b4e858c64c921c4273a174f6a13d9effa2fcc6

### input

- Ag1000g Zarrs
- Agam.P4.12 gff3

### Code 

- snakefile
- analysis/scripts/Zarr_to_LDNe.py

### Output

- Effective pop sizes .png
- 
- analysis/Ne_manifest.tsv - contains all Ne_estimate data along with other sample metadata

### Methods

LDNe was run on groups of individuals for collection location, species, and year that has greater than 8 samples. We extracted bi-allelic SNPs from each population, in non-coding and non-pericentromeric regions, at a MAF threshold of 0.05. Given that LDNe expects unlinked loci, we performed one iteration of LD-pruning (500 SNP window, 250 SNP step, r**=0.2). We then randomly sampled 20,000 SNPs from the remaining SNPs, and supplied this in FSTAT format to LDNe. 


#### Notes - Negative NE estimates from LDNe 
Why do you get Negative/infinity Ne estimates sometimes? 

(Waples, Do , 2010) evol appl. 

"As shown in eqn (2a), before estimating Ne in the LD method, the expected contribution of sampling error is subtracted from the empirical XX 
If Ne is large, or if only limited data are available, by chance mean can be smaller than the sample size correction, in which case the estimate of Ne will be negative. A related phenomenon can occur with the standard temporal method (Nei and Tajima 1981; Waples 1989) and with unbiased estimators of genetic differentiation (Nei 1978; Weir and Cockerham 1984). Negative estimates occur when the genetic results can be explained entirely by sampling error without invoking any genetic drift, so the biological interpretation is  =∞ (Laurie-Ahlberg and Weir 1979; Nei and Tajima 1981). In this situation, the user can conclude that the data provide no evidence that the population is not ‘very large’. However, even if the point estimate is negative, if adequate data are available the lower bound of the CI generally will be finite and can provide useful information about plausible limits Ne."
