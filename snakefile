import pandas as pd

chroms=['3L', '3R', '2L', '2R', 'X']
#read in names
LDNe_list = pd.read_csv("data/Phase3.LDNe.list")
LDNe_pops = LDNe_list['pop']

rule all:
	input:
		expand("analysis/LDNe/Ag_LDNe_{name}.{chrom}.out", name=LDNe_pops, chrom=chroms)

############## LDNe ############
############### noncoding, downsample to genepop and dat using the Zarr ############

rule Zarr2Dat:
	output:
		expand("analysis/LDNe/batch/{name}.{{chrom}}.batch.txt", name=LDNe_pops),
		expand("data/dat/{name}.{{chrom}}.dat", name=LDNe_pops)
	log:
		"logs/Zarr2Dat/Zarr_to_LDNe_{chrom}.log"
	params:
		nSNPs = 20000,
		localpath = "/home/sanj/ag1000g/data/phase3/",
		manifest = "data/Ag1000g.phase3.manifest.full.tsv",
		gff = "data/An.gambiae-PEST-BASEFEATURES_agamP4.12.gff3.gz"
	shell:
		"python analysis/scripts/Zarr_to_LDNe.py --n {params.nSNPs} --chroms {wildcards.chrom} --localpath {params.localpath} --gff {params.gff} --manifest {params.manifest} 2> {log}"

rule RunLDNe:
	input:
		dat="data/dat/{name}.{chrom}.dat",
		batch="analysis/LDNe/batch/{name}.{chrom}.batch.txt"
	output:
		"analysis/LDNe/Ag_LDNe_{name}.{chrom}.out"
	log:
		"logs/LDNe/{name}.{chrom}.log"
	shell:
		"analysis/scripts/Ne2-1L c:{input.batch} 2> {log}"
















############## IBDNE ######################

############# subset vcfs #############

rule subset_vcfs:
	output:
		"data/vcfs/{pop}_{chrom}.vcf.gz"
	log:
		"logs/subset/{pop}_{chrom}.log"
	params:
		region=lambda wildcards: '3L:2000000-41963435' if wildcards.chrom == '3L' else '3R:1-50000000',
		ac='1',
		compression='z',
		samples='data/list_samples/{pop}_sample_list'
	threads:8
	shell:
		"bcftools view --threads {threads} -r {params.region} -c {params.ac} -O {params.compression} -S {params.samples} /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.{wildcards.chrom}.vcf.gz > {output} 2> {log}"		

rule unzip_vcfs:
	input:
		"data/vcfs/{pop}_{chrom}.vcf.gz"
	output:
		temp("data/vcfs/{pop}_{chrom}.vcf")
	group:
		"ibd"
	shell:
		"gzip -d {wildcards.pop}_{wildcards.chrom}.vcf.gz > {output}"

rule ibd_segments:
	input:
		"data/vcfs/{pop}_{chrom}.vcf"
	output:
		"analysis/ibd/{pop}_{chrom}.ibd"
	log:
		"logs/ibd/{pop}_{chrom}_ibd.log"
	group:
		"ibd"
	threads:4
	shell:
		"java -jar ~/apps/ibdseq.r1206.jar gt={input} nthreads={threads} out=analysis/ibd/{wildcards.pop}_{wildcards.chrom} 2> {log}"

rule run_ibdne:
	input:
		gen_map = "data/Ag_genetic.map",
		ibd = expand("analysis/ibd/{{pop}}_{chrom}.ibd", chrom=chroms)
	output:
		"analysis/ibdne/{pop}_ibdne.ne"
	log:
		"logs/ibdne/{pop}.ibdne.log"
	threads:8
	shell:
		"cat {input.ibd} | java -jar ~/apps/ibdne.19Sep19.268.jar map={input.gen_map} out=analysis/ibdne/{wildcards.pop}_ibdne nthreads={threads} 2> {log}"
