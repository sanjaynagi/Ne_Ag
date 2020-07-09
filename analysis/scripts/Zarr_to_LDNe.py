#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import zarr
import allel
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Take Zarr, downsample and convert to .dat for LDNe')
parser.add_argument('--n', type=int, action='store', default=10000, help='Number of snps to random sample')
parser.add_argument('--localpath', type=str, action='store',help='Path to zarrs + metadata')
parser.add_argument('--manifest', type=str, action='store', help='Tab-delimited Samples metadata file')
parser.add_argument('--gff', type=str, action='store', help='Ag gff3 file')
parser.add_argument('--chroms', type=str, action='store', nargs='*', default=['3L','3R'], help='Which chromosomes to use')
args=parser.parse_args()

print("-------------------- Zarr to LDNe ----------------------")
manifest = pd.read_csv(args.manifest, sep="\t")
manifest.location = [loc.replace(" ", "") for loc in manifest.location] #remove whitespace
all_sets = manifest.sample_set.unique()

chroms = args.chroms
n = args.n  # number of SNPs to choose randomly

local_path = Path(f"{args.localpath}").expanduser() 

##### functions #####
def load_arrays_noncoding_and_centromeres(local_path, _set, chrom, coding_reg_df, sitefilter='gamb_colu'):
    
    """
    This function reads and filters a genotyping array to the noncoding, noncentromeric regions, and applys a filter depending on 
    whether the samples are arabiensis (arab) or gambiae/coluzzii (gamb_colu)
    """
    
    Ag_array = zarr.open_array(f"{local_path}/snp_genotypes/all/{_set}/{chrom}/calldata/GT/", mode = 'r')
    filters = zarr.open(f"{local_path}/site_filters/dt_20200416/{sitefilter}/{chrom}/variants/filter_pass", mode="r")
    positions = zarr.open_array(f"{local_path}/snp_genotypes/all/sites/{chrom}/variants/POS/", mode='r')
    positions = positions[:][filters[:]]    
    geno = allel.GenotypeDaskArray(Ag_array)
    geno = geno[filters[:]]
    
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
        
    positions = allel.SortedIndex(positions[centromere])
    #get boolean array for positions that are coding - allel.locate_ranges so fast!
    coding = positions.locate_ranges(coding_reg_df.start, coding_reg_df.end, strict=False)
    #compress to get noncoding SNPs and remove centromeric regions of low recombination
    #get non-centromeric regions. currently chosen by eye based on ag1000g phase1 paper fig1.
  
    #TODO currently using all non-coding regions, alternate option may be to take SNPs x distance from coding regions 
    geno = geno.compress(centromere, axis=0)
    geno = geno.compress(~coding, axis=0) #we want noncoding regions so '~' to get inverse of boolean
    positions = positions[~coding]
    
    return(geno, positions)

def convert2dat(gn, positions, n, _set, loc, yr, sp, chrom):
    
    """
    This function takes a genotyping array and positions, along with info on the sample_set, location , year, species, chrom,
    and converts the genotyping array into a .dat file, suitable for use with LDNe. It randomly downsamples the genotyping array
    by n SNPs
    """
    
    # MAF 0.05 filter
    ac = gn.count_alleles()
    freqs = ac.to_frequencies().compute()
    ALT1 = freqs[:,1] > 0.05
    ALT2 = freqs[:,2] > 0.05
    ALT3 = freqs[:,3] > 0.05
    maf_flt = np.logical_or(ALT1, ALT2, ALT3)
    gn = gn.compress(maf_flt, axis=0)

    #take random sample of n SNPs 
    snp_sample = np.random.choice(gn.shape[0], n, replace=False)
    snp_sample.sort()
    gnr = gn.take(snp_sample, axis=0)
    gnr = np.array(gnr[:])
    gnr = gnr.astype(str)

    pos = positions[maf_flt]
    pos = pos[snp_sample]
    prefix = f'{chrom}_'
    pos_string = [prefix + p for p in pos.astype(str)]

    gnr[gnr == '-1'] = '00' #convert missing alleles 
    dat = np.empty([gnr.shape[0], gnr.shape[1]])

    #join genotypes in same individual 
    print(f"Converting to .dat {_set} {loc} {yr} {sp} {chrom}")
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
    popnames = np.repeat(f"{_set}{loc}{yr}{sp}{chrom}", gnr.shape[1])
    dat = np.column_stack((popnames, dat.T)) #

    #write out .dat file for LDNe 
    with open(f'data/dat/{_set}.{loc}.{yr}.{sp}.{chrom}.dat', 'w') as datfile:
        datfile.write(f'{gnr.shape[1]}\t{gnr.shape[0]}\t4\t2\n')
        datfile.write("\n".join("".join(map(str, x)) for x in pos_string)) 
        datfile.write("\n")
        datfile.write("\n".join("\t".join(map(str, x)) for x in dat))
        
    

################

#### main #####


for chrom in chroms:

    print(f"------------------------ LDNe_Ag -------------------------------\n")
    print(f"Storing metadata in data/Phase3.LDNe.tsv")
    with open('data/Phase3.LDNe.tsv', 'w') as metafile:
        metafile.write("sample_set\tlocation\tyear\tspecies\tchromosome\n")
 
        
    #filter the gff3 to be coding and regulatory regions
    df = allel.gff3_to_dataframe(f"{args.gff}")
    coding_reg_df = df[~df.type.isin(['chromosome', 'three_prime_UTR','five_prime_UTR',
                    'mRNA', 'CDS', 'exon'])].drop(columns=['source', 'strand', 'phase', 'score'])
    coding_reg_df = coding_reg_df[coding_reg_df.seqid == chrom]

    for _set in all_sets:
    
        #subset metadata
        metadata = manifest[manifest.sample_set == _set].reset_index(drop=True)

        ### loop through combos 
        for loc in metadata.location.unique():

            nmeta = metadata[metadata.location == loc]

            for yr in nmeta.year.unique():

                nmeta2 = nmeta[nmeta.year == yr]
                
                #have edited .species_gambiae_coluzzii column to contain arabiensis instead of NA 
                for sp in nmeta2.species_gambiae_coluzzii.unique():
                    
                    # file exists so ignore and skip
                    myfile = Path(f"data/dat/{_set}.{loc}.{yr}.{sp}.{chrom}.dat")
                    if myfile.is_file():
                        continue

                    #if there is less than 9 samples than skip
                    if (nmeta2.species_gambiae_coluzzii == sp).sum() <= 8:
                        continue
                    
                    if sp == 'arabiensis':  
                        geno, positions = load_arrays_noncoding_and_centromeres(local_path,_set, chrom, coding_reg_df, sitefilter='arab')
                        #filter to species 
                        nmeta3 = nmeta2[nmeta2.species_gambiae_coluzzii == sp]
                        flt = np.array(nmeta3.index)
                        #filter to correct loc, year, species individuals
                        gn = geno.take(flt, axis=1)
                        print(f"\nProducing LDNe input for {loc}, {yr}, {sp}, {chrom}. {nmeta3.shape[0]} individuals, (arabiensis filter)")

                        convert2dat(gn, positions, n, _set, loc, yr, sp, chrom)
                
                    else:
                        geno, positions = load_arrays_noncoding_and_centromeres(local_path, _set, chrom, coding_reg_df, sitefilter='gamb_colu')
                        #filter to species 
                        nmeta3 = nmeta2[nmeta2.species_gambiae_coluzzii == sp]
                        flt = np.array(nmeta3.index)
                        #filter to correct loc, year, species individuals
                        gn = geno.take(flt, axis=1)
                        print(f"\nProducing LDNe input for {loc}, {yr}, {sp}, {chrom}. {nmeta3.shape[0]} individuals, (gambcolu filter)")

                        convert2dat(gn, positions, n, _set, loc, yr, sp, chrom)

                    #write metadata file for samples that are included
                    with open('data/Phase3.LDNe.tsv', 'a') as metafile:
                        metafile.write(f'{_set}\t{loc}\t{yr}\t{sp}\t{chrom}\n')
                    
                    print("Writing .batch.txt file for LDNe...")
                    with open(f'analysis/LDNe/batch/{_set}.{loc}.{yr}.{sp}.{chrom}.batch.txt', 'w') as batch_file:
                        batch_file.write(f'1\t0\n1\n0.05\t-1\n15\t0\t1\n1\n0\n0\n0\n0\nanalysis/LDNe/Ag_LDNe_{_set}.{loc}.{yr}.{sp}.{chrom}.out\n')
                        batch_file.write(f"data/dat/{_set}.{loc}.{yr}.{sp}.{chrom}.dat\n")
                        batch_file.write("*")
                        batch_file.close()