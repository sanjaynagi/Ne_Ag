import pandas as pd
samples = pd.read_csv("~/ag1000g/data/samples.meta.txt", sep="\t")
pops = samples.population.unique()
chroms=['2L', '2R','3L', '3R', 'X']

configfile:"config.yaml"

rule all:
	input:
		expand("analysis/LDNe/Ag_LDNe_{pop}_{chrom}.out", pop=pops, chrom=chroms)

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

############## LDNe ############
############### noncoding, downsample to genepop and dat using the Zarr ############

rule Zarr_to_LDNe:
	input:
		zarr = config['zarr'],
		gff = config['gff'],
		samples = config['samples']
	output:
		expand("data/dat/{pop}_{chrom}.dat", pop=pops, chrom=chroms)
	log:
		"logs/Zarr_to_LDNe/Zarr_to_LDNe.log"
	params:
		nSNPs = 20000,
	shell:
		"python analysis/scripts/Zarr_to_LDNe.py --n {params.nSNPs} --chroms {chroms} --zarr {input.zarr} --gff {input.gff} --samples {input.samples} 2> {log}"

rule run_ldne:
	input:
		dat="data/dat/{pop}_{chrom}.dat",
		batch="analysis/LDNe/batch/ag_batch_{pop}_{chrom}.txt"
	output:
		"analysis/LDNe/Ag_LDNe_{pop}_{chrom}.out"
	log:
		"logs/ldne/{pop}_{chrom}.log"
	shell:
		"~/apps/NeEstimator/Ne2-1L c:{input.batch} 2> {log}"
#move NeEstimator to analysis/scripts 
