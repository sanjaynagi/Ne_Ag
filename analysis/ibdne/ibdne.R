#### IBDNe ####
library(tidyverse)
library(glue)
library(data.table)
library(scales)

####### geometric mean function ###########
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
###############
samples = fread("~/ag1000g/data/samples.meta.txt")
pops = unique(samples$population)

#################
pops = pops[!pops %in% c('GQgam','GNcol','GHgam','CMgam')]

ibdne = list()

for (pop in pops){
  ibdne[[pop]] = fread(glue("ibdne/{pop}_ibdne.ne"))
}

all_ibdne = bind_rows(ibdne, .id = "pop")

hm_mean = function(a){
  b = 1/mean(1/a)
  return(b)
}
#### calculate geometric means for IBDNe estimates
IBDNe_harmonic_means = c()
for (pop in pops){
  Ne =hm_mean(ibdne[[pop]]$NE)
  IBDNe_harmonic_means = c(IBDNe_harmonic_means, Ne)
}
names(IBDNe_harmonic_means) = names(ibdne)

IBDNe_harmonic_means = data.frame(IBDNe_harmonic_means) %>% rownames_to_column('Population')
IBDNe_harmonic_means = IBDNe_harmonic_means[order(IBDNe_harmonic_means$IBDNe_harmonic_means),]

fwrite(IBDNe_harmonic_means, "ibdne/IBDNe_Harmonic_means.txt", sep="\t", row.names = FALSE)




### plots ####
ns = c(50,200)
pdf("ibdne_plots.pdf")
for (population in pops){
  for (n in ns){
    
    a = all_ibdne %>% filter(`pop` == population) %>% filter(`GEN` >= 4 & `GEN` <= n)
  
    print(ggplot(a, aes(x=`GEN`, y=`NE`))+ 
      geom_point() +
      scale_y_continuous(breaks = round(seq(0, max(a$`UPR-95%CI`), by = max(a$`UPR-95%CI`)/10),0) ,labels = comma) +
      geom_ribbon(aes(ymin=a$`LWR-95%CI`, ymax=a$`UPR-95%CI`),fill='blue', alpha=0.3) +
      ggtitle(glue("{population} 4 to {n} generations ago")) +
      theme_light())
  }
}
dev.off()
  
### plots ####
ns = c(50,200)
pdf("ibdne_plots_no_CIs.pdf")
for (population in pops){
  for (n in ns){
    
    a = all_ibdne %>% filter(`pop` == population) %>% filter(`GEN` >= 4 & `GEN` <= n)
    
    print(ggplot(a, aes(x=`GEN`, y=`NE`))+ 
            geom_point() +
            scale_y_continuous(breaks = round(seq(0, max(a$NE), by = max(a$NE)/10),0) ,labels = comma) +
            ggtitle(glue("{population} 4 to {n} generations ago")) +
            theme_light())
  }
}
dev.off()

  
all_ibdne %>% filter(`pop` == 'FRgam') %>% 
  ggplot(., aes(x=GEN, y=NE)) + 
  geom_line() + 
  scale_y_continuous(labels = comma) +
  geom_ribbon(aes(ymin=`LWR-95%CI`, ymax=`UPR-95%CI`),fill='blue', alpha=0.3) +
  theme_light()



all_ibdne %>% filter(`pop` == 'AOcol') %>% 
  ggplot(., aes(x=GEN, y=NE)) + 
  geom_line() + 
  scale_y_continuous(labels = comma) +
  geom_ribbon(aes(ymin=`LWR-95%CI`, ymax=`UPR-95%CI`),fill='blue', alpha=0.3) +
  theme_light()

all_ibdne %>% filter(`pop` == 'GHcol') %>% 
  ggplot(., aes(x=GEN, y=NE)) + 
  geom_line() + 
  scale_y_continuous(labels = comma) +
  geom_ribbon(aes(ymin=`LWR-95%CI`, ymax=`UPR-95%CI`),fill='blue', alpha=0.3) +
  theme_light() +
  ggtitle("GHcol")

all_ibdne %>% filter(`pop` == 'GW') %>% 
  ggplot(., aes(x=GEN, y=NE)) + 
  geom_line() + 
  scale_y_continuous(labels = comma) +
  geom_ribbon(aes(ymin=`LWR-95%CI`, ymax=`UPR-95%CI`),fill='blue', alpha=0.3) +
  theme_light() + ggtitle("GW")
ggsave("ibdne/plots/GW.png")
  #samples[samples$population == 'GW']

all_ibdne %>% filter(`pop` == 'UGgam') %>% 
  ggplot(., aes(x=GEN, y=NE)) + 
  geom_line() + 
  scale_y_continuous(labels = comma)+
  geom_ribbon(aes(ymin=`LWR-95%CI`, ymax=`UPR-95%CI`),fill='blue', alpha=0.3) +
  theme_light() +
  ggtitle("UGgam")
ggsave("UGgam.png")


all_ibdne %>% filter(`pop` == 'BFcol') %>% 
  ggplot(., aes(x=GEN, y=NE)) + 
  geom_line() + 
  scale_y_continuous(labels = comma) +
  geom_ribbon(aes(ymin=`LWR-95%CI`, ymax=`UPR-95%CI`),fill='blue', alpha=0.3) +
  theme_light() + ggtitle("BFcol")
ggsave("ibdne/plots/CIcol.png")