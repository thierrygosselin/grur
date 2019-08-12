library(strataG)
library(tidyverse)

files <- dir("test", pattern = "_genotypes_", full.names = TRUE)

ne <- do.call(rbind, lapply(files, function(f) {
  load(f)
  do.call(rbind, lapply(sim.list, function(rep) {
    g <- df2gtypes(rep$rms, ploidy = 2)
    as.data.frame(ldNe(g)) %>% 
      mutate(scenario = sc$scenario, replicate = rep$rep) %>% 
      select(scenario, replicate, stratum, Ne, param.lci, param.uci)
  }))
}))
