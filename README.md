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


Subset the vcfs, and remove variants that are no longer segregating (-c1), compress (-z) and exclude pericentromeric regions.

3R: 1-5000000
3L: 2000000-4196nnnn

```bash

for i in KE UGgam FRgam GQgam; 
	do
		~/apps/bcftools/bin/bcftools view -c1 -O z /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.3R.vcf.gz -S ../data/list_samples/${i}_sample_list -r 3R:1-50000000 > ../data/${i}_3R.vcf 2> ${i}_3R_vcf.stderr ; 

		~/apps/bcftools/bin/bcftools view -c1 -O z /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.3L.vcf.gz -S ../data/list_samples/${i}_sample_list -r 3L:2000000-41963435 > ../data/${i}_3L.vcf 2> ${i}_3L_vcf.stderr;
done


```

Then loop again over the populations with ibdseq to calculate shared pairwise segments of Identity by descent. 

```bash
for i in KE UGgam FRgam EQgam; 
	do
		java -jar ~/apps/ibdseq.r1206.jar gt=../data/${i}_3R_vcf out=${i}_ibd_3R 2> ${i}_3R_ibd.stderr ; 

		java -jar ~/apps/ibdseq.r1206.jar gt=../data/${i}_3L_vcf out=${i}_ibd_3L 2> ${i}_3L_ibd.stderr ; 

		cat ${i}_ibd_3R ${i}_ibd_3L > ${i}_ibd
done


```
We then take the .ibd file and run it through IBDNe, supplying the genetic map produced in Miles (2017). The output are then used in RStudio to plot and explore. 

```bash
for i in KE UGgam FRgam EQgam; 
	do
		cat ${i}_ibd | java -jar ~/apps/ibdne.19Sep19.268.jar map=recombination_maps/Ag_genetic.map out=${i}_ibdne 2> ${i}_ibdne.stderr; 
done

```




