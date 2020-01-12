pops=['KE', 'UGgam', 'FRgam', 'GNgam', 'GHcol']
chroms=['3L', '3R']


rule all:
	input: 
		expand("analysis/ibdne/{pop}_ibdne.ne", pop=pops)

rule all_ldne:
	input:
		expand("Ag_LDNe_{pop}_{chrom}.out", pop=pops, chrom=chroms)


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
		"~/apps/bcftools/bin/bcftools view --threads {threads}-r {params.region} -c {params.ac} -O {params.compression} -S {params.samples} /kwiat/vector/ag1000g/release/phase2.AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.{wildcards.chrom}.vcf.gz > {output} 2> {log}"		

################ IBD ####################
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

################ Run IBDne #####################
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
		"cat {input.ibd} | java -jar ~/apps/ibdne.19Sep19.268.jar map={input.gen_map} out={wildcards.pop}_ibdne nthreads={threads} 2> {log}"


############### noncoding, downsample to genepop and dat ############
rule restrict_noncoding:
	input:
		"data/vcfs/{pop}_{chrom}.vcf.gz"
	output:
		pipe("data/vcfs/noncoding/{pop}_{chrom}_noncoding.vcf.gz")
	log:
		"logs/gatk_noncoding/{pop}_{chrom}.log"
	group:
		'downsample'
	shell:
		"gatk SelectVariants -V {input} -O {output} -XL ~/reference/regions/Ag_coding_and_regulatory.intervals 2> {log}"

rule downsample:
	input:
		"data/vcfs/{pop}_{chrom}_noncoding.vcf.gz"
	output:
		"data/vcfs/noncoding/downsample/{pop}_{chrom}_random.vcf"
	log:
		"logs/shuf_downsample/{pop}_{chrom}.log"
	group:
		'downsample'
	shell:
		"gzip -d {input} | shuf -n 10000 > {output} 2> {log}"

rule vcf2genepop:
	input:
		"data/vcfs/noncoding/downsample/{pop}_{chrom}_random.vcf"
	output:
		"data/genepops/{pop}_{chrom}.gen"
	log:
		"logs/vcf2genepop/{pop}_{chrom}.log"
	group:
		"convert"
	shell:
		"perl analysis/vcf2genepop.pl vcf={input} pop={wildcards.pop} > {output} 2> {log}"

rule genepop2dat:
	input:
		"data/genepops/{pop}_{chrom}.gen"
	output:
		"data/dat/{pop}_{chrom}.dat",
		temp("loci_{pop}_{chrom}"),
		temp("fline_{pop}_{chrom}"),
		temp("gen_{pop}_{chrom}")
	group:
		"convert"
	log:
		"logs/genepop2dat/{pop}_{chrom}.log"
	shell:
		"(Rscript analysis/genepop2dat.R {input} {wildcards.pop}_{wildcards.chrom} {output}) 2> {log}"


####################################
############## Run LDNe ############
rule create_batch_file:
	output:
		"analysis/LDNe/ag_batch_{pop}_{chrom}.txt"
	group:
		"LDNe"
	run:
		with open(f'analysis/LDNe/ag_batch_{wildcards.pop}_{wildcards.chrom}.txt', 'w') as batch_file:
			batch_file.write(f'1\t0\n3\n0.05\t0.02\t0.01\n15\t0\t1\n1\n0\n0\n0\n0\nanalysis/LDNe/Ag_LDNe_{wildcards.pop}_{wildcards.chrom}.out\n')
			batch_file.write(f"data/dat/{wildcards.pop}_{wildcards.chrom}.dat\n")
			batch_file.write("*")
			batch_file.close()

rule run_ldne:
	input:
		dat="data/dat/{pop}_{chrom}.dat",
		batch="analysis/LDNe/ag_batch_{pop}_{chrom}.txt"
	output:
		"Ag_LDNe_{pop}_{chrom}.out"
	group:"LDNe"
	log:
		"logs/ldne/{pop}_{chrom}.log"
	shell:
		"~/apps/NeEstimator/Ne2-1L c:{input.batch} 2> {log}"
