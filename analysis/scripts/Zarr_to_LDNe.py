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
parser.add_argument('--samples', type=str, action='store', help='Tab-delimited Samples metadate file')
parser.add_argument('--gff', type=str, action='store', help='Ag gff3 file')
parser.add_argument('--chroms', type=str, action='store', default=['3L','3R'], help='Tab-delimited Samples metadate file')
parser.add_argument('--pops', type=str, action='store', help='Tab-delimited Samples metadate file')
args=parser.parse_args()

print("------------ Zarr to LDNe -------------")
# As the current method of running LDNe in the snakemake pipeline requires subsetting the VCF (very slow), 
# I will instead write a script to convert the zarr to .dat format for LDNe.
samples = pd.read_csv(args.samples, sep="\t")
pops = args.pops
chroms = args.chroms
n = args.n  # number of SNPs to choose randomly

for chrom in chroms:

    #open arrays
    print(f"Opening arrays {pop} {chrom}")
    Ag_store = zarr.open_array(f"{args.zarr}/{chrom}/calldata/GT/", mode = 'r')
    positions = allel.SortedIndex(zarr.open_array(f"{args.zarr}/{chrom}/variants/POS", mode='r'))
    ag_geno = allel.GenotypeChunkedArray(Ag_store)

    #filter the gff3 to be coding and regulatory regions
    df = allel.gff3_to_dataframe(f"{args.gff}")
    coding_reg_df = df[~df.type.isin(['chromosome', 'three_prime_UTR','five_prime_UTR',
                      'mRNA', 'CDS', 'exon'])].drop(columns=['source', 'strand', 'phase', 'score'])
    coding_reg_df = coding_reg_df[coding_reg_df.seqid == chrom]

    centromere = (positions > 2000000) if chrom == '3L' else (positions > 50000000) #filter to remove centromeres
    positions = positions[centromere]
    #get boolean array for positions that are coding - allel.locate_ranges so fast!
    coding = positions.locate_ranges(coding_reg_df.start, coding_reg_df.end, strict=False)

    #compress to get noncoding SNPs
    ag_geno = ag_geno.compress(centromere, axis=0)
    ag_geno = ag_geno.compress(~coding, axis=0) #we want noncoding regions so '~' to get inverse of boolean

    # Remove centromeric regions of low recombination
    # Converting to .dat format for LDNe

    for pop in pops:
        pop_bool = samples.population == pop
        pop_geno = ag_geno.compress(pop_bool, axis=1)

        ac = pop_geno.count_alleles()
        seg = ac.is_segregating()
        pop_geno = pop_geno.compress(seg, axis=0)

        vidx = np.random.choice(pop_geno.shape[0], n, replace=False)
        vidx.sort()
        gnr = pop_geno.take(vidx, axis=0)
        gnr = np.array(gnr[:])
        gnr = gnr.astype(str)

        positions = positions[~coding]
        positions = positions[seg]
        positions = positions[vidx]
        prefix = f'{chrom}_'
        pos_string = [prefix + pos for pos in positions.astype(str)]

        gnr[gnr == '-1'] = '00' #convert missing alleles 
        dat = np.empty([gnr.shape[0], gnr.shape[1]])

        print(f"Converting to .dat {pop} {chrom}")
        for x in range(gnr.shape[0]):
            for y in range(gnr.shape[1]):
                dat[x,y] = ''.join(gnr[x,y])

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

        popnames = np.repeat(f"{pop}_{chrom}", n)
        dat = np.column_stack((popnames, dat)) #
        
        with open(f'{pop}_{chrom}.dat', 'w') as datfile:
            datfile.write(f'{gnr.shape[1]}\t{gnr.shape[0]}\t4\t2\n')
            datfile.write("\n".join("".join(map(str, x)) for x in pos_string)) 
            datfile.write("\n")
            datfile.write("\n".join(" ".join(map(str, x)) for x in dat))    
