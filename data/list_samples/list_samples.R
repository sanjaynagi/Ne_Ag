#phase2metadata
library(data.table)
library(tidyverse)
library(glue)

samples = fread("samples.meta.txt")

#write out each pop individually
for (pop in unique(samples$population)){
  samples %>% filter(`population` == pop) %>% 
    select(1) %>% 
    fwrite(., glue("list_samples/{pop}_sample_list"), col.names = FALSE)
}

table(samples$population)

    colnames(samples)

#for UMAP
samples %>% select(1,3) %>% 
  fwrite(., "samples_pops", col.names = FALSE, sep = '\t')

samples %>% 
  select(4,3) %>% 
  mutate('country' = paste0(country,'_', population)) %>% 
  unique() %>% 
  fwrite(., "samples_pops_desc", col.names = FALSE, sep = '\t')


  
