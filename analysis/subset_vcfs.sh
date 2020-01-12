#!/bin/bash
# subset chromosomes, calculate segments of IBD and calculate Ne using IBDNe

for i in GNgam; 
	do
		~/apps/bcftools/bin/bcftools view -c1 /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.3L.vcf.gz -S ../data/list_samples/${i}_sample_list -r 3L:2000000-41963435 > ../data/${i}_3L.vcf 2> ${i}_3L_vcf.stderr;
done

