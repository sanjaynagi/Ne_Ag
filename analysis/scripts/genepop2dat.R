#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glue)
library(data.table)
library(tidyverse)

#This is to supply arguments to for snakemake and reproducibility 

gen = read.table(args[1], skip=1, colClasses = 'character')

gen %>% colnames(.) %>% 
  str_remove(., "X") %>% 
  str_remove(., "[.]") %>% 
  list() %>% 
  fwrite(., glue("loci_{args[2]}"), sep = "\t")

gen = gen %>% rownames_to_column('sample')
gen$sample = glue("{args[2]}")

fwrite(gen, glue("gen_{args[2]}"), sep="\t", col.names = FALSE)
gen = gen %>% select(-`sample`) 

paste(nrow(gen), ncol(gen), 04, 2) %>% list() %>% fwrite(., glue("fline_{args[2]}"), sep="\t")

system(glue("cat fline_{args[2]} loci_{args[2]} gen_{args[2]} > {args[3]}"))



