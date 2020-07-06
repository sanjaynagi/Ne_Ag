#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import zarr
import allel
import argparse

parser = argparse.ArgumentParser(description='Take Zarr, downsample and convert to .dat for LDNe')
parser.add_argument('--n', type=int, action='store', default=10000, help='Number of snps to random sample')
parser.add_argument('--zarr', type=str, action='store',help='Path to zarr file')
parser.add_argument('--samples', type=str, action='store', help='Tab-delimited Samples metadata file')
parser.add_argument('--gff', type=str, action='store', help='Ag gff3 file')
parser.add_argument('--chroms', type=str, action='store', nargs='*', default=['3L','3R'], help='Which chromosomes to use')
args=parser.parse_args()

print("------------ Zarr to LDNe -------------")
# As the current method of running LDNe in the snakemake pipeline requires subsetting the VCF (very slow), 
# I will instead write a script to convert the zarr to .dat format for LDNe.
samples = pd.read_csv(args.samples, sep="\t")
pops = samples.population.unique()
chroms = args.chroms
n = args.n  # number of SNPs to choose randomly

for chrom in chroms:

    #open arrays
    print(f"Opening arrays {chrom}")
    Ag_store = zarr.open_array(f"{args.zarr}/{chrom}/calldata/GT/", mode = 'r')
    positions = allel.SortedIndex(zarr.open_array(f"{args.zarr}/{chrom}/variants/POS", mode='r'))
    ag_geno = allel.GenotypeChunkedArray(Ag_store)

    #filter the gff3 to be coding and regulatory regions
    df = allel.gff3_to_dataframe(f"{args.gff}")
    coding_reg_df = df[~df.type.isin(['chromosome', 'three_prime_UTR','five_prime_UTR',
                      'mRNA', 'CDS', 'exon'])].drop(columns=['source', 'strand', 'phase', 'score'])
    coding_reg_df = coding_reg_df[coding_reg_df.seqid == chrom]

    #get non-centromeric regions. currently chosen by eye based on ag1000g phase1 paper fig1.
    if chrom == '2L':
        centromere = (positions > 3000000)
    elif chrom == '2R':
        centromere = (positions < 57000000)
    elif chrom == '3L':
        centromere = (positions > 2000000)
    elif chrom == '3R':
        centromere = (positions < 50000000)
    elif chrom == 'X':
        centromere = (positions < 21000000)

    positions = positions[centromere]
    #get boolean array for positions that are coding - allel.locate_ranges so fast!
    coding = positions.locate_ranges(coding_reg_df.start, coding_reg_df.end, strict=False)

    #compress to get noncoding SNPs and remove centromeric regions of low recombination
    #TODO currently using all non-coding regions, alternate option may be to take SNPs x distance from coding regions 
    ag_geno = ag_geno.compress(centromere, axis=0)
    ag_geno = ag_geno.compress(~coding, axis=0) #we want noncoding regions so '~' to get inverse of boolean
    print(f"Filtering out coding regions {chrom}")

    # Converting to .dat format for LDNe
    for pop in pops:
        pop_bool = samples.population == pop
        pop_geno = ag_geno.compress(pop_bool, axis=1)

        # MAF 0.05 filter
        ac = pop_geno.count_alleles()
        freqs = ac.to_frequencies()[:]
        ALT1 = freqs[:,1] > 0.05
        ALT2 = freqs[:,2] > 0.05
        ALT3 = freqs[:,3] > 0.05
        maf_flt = np.logical_or(ALT1, ALT2, ALT3)

        pop_geno = pop_geno.compress(maf_flt, axis=0)

        #take random sample of n SNPs 
        snp_sample = np.random.choice(pop_geno.shape[0], n, replace=False)
        snp_sample.sort()
        gnr = pop_geno.take(snp_sample, axis=0)
        gnr = np.array(gnr[:])
        gnr = gnr.astype(str)

        pos = positions[~coding]
        pos = pos[seg]
        pos = pos[vidx]
        prefix = f'{chrom}_'
        pos_string = [prefix + p for p in pos.astype(str)]

        gnr[gnr == '-1'] = '00' #convert missing alleles 
        dat = np.empty([gnr.shape[0], gnr.shape[1]])

        #join genotypes in same individual 
        print(f"Converting to .dat {pop} {chrom}")
        for x in range(gnr.shape[0]):
            for y in range(gnr.shape[1]):
                dat[x,y] = ''.join(gnr[x,y])

        #convert to .dat format genotypes 
        dat = dat.astype(str)
        dat[dat == '0.0'] = '0101'
        dat[dat == '1.0'] = '0102'
        dat[dat == '2.0'] = '0103'
        dat[dat == '3.0'] = '0104'
        dat[dat == '10.0'] = '0201'
        dat[dat == '11.0'] = '0202'
        dat[dat == '12.0'] = '0203'
        dat[dat == '13.0'] = '0204'
        dat[dat == '20.0'] = '0301'
        dat[dat == '21.0'] = '0302'
        dat[dat == '22.0'] = '0303'
        dat[dat == '23.0'] = '0304'
        dat[dat == '30.0'] = '0401'
        dat[dat == '31.0'] = '0402'
        dat[dat == '32.0'] = '0403'
        dat[dat == '33.0'] = '0404'

        #stack population name and transposed genotypes 
        popnames = np.repeat(f"{pop}_{chrom}", gnr.shape[1])
        dat = np.column_stack((popnames, dat.T)) #
        
        #write out .dat file for LDNe 
        with open(f'data/dat/{pop}_{chrom}.dat', 'w') as datfile:
            datfile.write(f'{gnr.shape[1]}\t{gnr.shape[0]}\t4\t2\n')
            datfile.write("\n".join("".join(map(str, x)) for x in pos_string)) 
            datfile.write("\n")
            datfile.write("\n".join("\t".join(map(str, x)) for x in dat))    
