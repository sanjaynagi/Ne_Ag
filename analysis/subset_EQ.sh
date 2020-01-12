#!/bin/bash
# subset chromosomes, calculate segments of IBD and calculate Ne using IBDNe

~/apps/bcftools/bin/bcftools view -c1 -O z /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.3R.vcf.gz -S ../data/list_samples/GQgam_sample_list -r 3R:1-50000000 > ../data/GQgam_3R.vcf 2> GQgam_3R_vcf.stderr ; 

~/apps/bcftools/bin/bcftools view -c1 -O z /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.3L.vcf.gz -S ../data/list_samples/GQgam_sample_list -r 3L:2000000-41963435 > ../data/GQgam_3L.vcf 2> GQgam_3L_vcf.stderr;


