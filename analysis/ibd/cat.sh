#!/bin/bash
# subset chromosomes, calculate segments of IBD and calculate Ne using IBDNe

for i in GNgam; 
	do
		cat ${i}_ibd_3R.ibd ${i}_ibd_3L.ibd > ${i}.ibd
done

