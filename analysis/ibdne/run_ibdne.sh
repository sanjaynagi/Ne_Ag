#!/bin/bash
# subset chromosomes, calculate segments of IBD and calculate Ne using IBDNe

for i in GNgam; 
	do
		cat ibd/${i}.ibd | java -jar ~/apps/ibdne.19Sep19.268.jar map=~/recombination_maps/Ag_genetic.map out=${i}_ibdne 2> ${i}_ibdne.stderr; 
done

