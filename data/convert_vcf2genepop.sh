#!/bin/bash
#convert_vcf2genepop
set -e ; set -u ; set -o pipefail ;


for pop in KE UGgam FRgam GNgam; do
	for chrom in 3L 3R;
do echo initiating ${pop} ${chrom} ;perl vcf2genepop.pl vcf=noncoding/downsample/${pop}_${chrom}_random.vcf pop=${pop} > noncoding/downsample/genepops/${pop}_${chrom}.genepop ; done; done 
