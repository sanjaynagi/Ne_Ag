#!/bin/bash
# subset chromosomes, calculate segments of IBD and calculate Ne using IBDNe

for i in GNgam; 
	do
		java -jar ~/apps/ibdseq.r1206.jar gt=../data/${i}_3R.vcf out=${i}_ibd_3R 2> ${i}_3R_ibd.stderr && java -jar ~/apps/ibdseq.r1206.jar gt=../data/${i}_3L.vcf out=${i}_ibd_3L 2> ${i}_3L_ibd.stderr
done

