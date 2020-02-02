---
output:
  pdf_document: default
  html_document: default
---
## Ag1000g - IBD and effective population size work

Analysis of Ag1000g populations looking at effective pop sizes and ibdne. The idea is to use compare Ne estimates between Kenya, outgroup and island populations.

Kenya expected to be smaller due to ROH etc

Populations 

- Kenya
- Mayotte
- Bioko
- Uganda


#### input

Data files:

* 'ag1000g/vcf.gz' for each chromosome
* 'ag1000g/list_samples' list of samples for each population within the VCFs

#### output

We aim to estimate recent effective pop sizes using three methods - LDNe , NB , IBDNe





## Code 



## Negative NE estimates from LDNe 
(Waples,Do , 2010) evol appl. 

As shown in eqn (2a), before estimating Ne in the LD method, the expected contribution of sampling error is subtracted from the empirical An external file that holds a picture, illustration, etc. Object name is eva0003-0244-mu73.jpg. If Ne is large, or if only limited data are available, by chance mean An external file that holds a picture, illustration, etc. Object name is eva0003-0244-mu74.jpg can be smaller than the sample size correction, in which case the estimate of Ne will be negative. A related phenomenon can occur with the standard temporal method (Nei and Tajima 1981; Waples 1989) and with unbiased estimators of genetic differentiation (Nei 1978; Weir and Cockerham 1984). Negative estimates occur when the genetic results can be explained entirely by sampling error without invoking any genetic drift, so the biological interpretation is An external file that holds a picture, illustration, etc. Object name is eva0003-0244-mu75.jpg = ∞ (Laurie-Ahlberg and Weir 1979; Nei and Tajima 1981). In this situation, the user can conclude that the data provide no evidence that the population is not ‘very large’. However, even if the point estimate is negative, if adequate data are available the lower bound of the CI generally will be finite and can provide useful information about plausible limits Ne. 
