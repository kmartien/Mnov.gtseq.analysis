library(vcfR)
library(vcftoolsR)
library(tidyverse)
library(dplyr)
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC")
source("R/functions/mplot2tgt.R")
source("R/functions/Compare.replicates.R")
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.analysis")

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 20
num.locs <- 322
min.genos.per.ind <- num.locs * 0.6

#data_sets <- list('RunMS58', 'RunMS51', 'RunMS45.90s', 'RunMS54', 'GTseq.val', 'GTseq.prod')
data_sets <- list('RunMS58.replicates.merged', 
                  'RunMS58.MS51.merged', 
                  'RunMS58.MS51.MS45.90s.merged',
                  'RunMS58.MS51.MS45.90s.MS54.merged',
                  'RunMS58.MS51.MS45.90s.MS54.GTseq.val.merged',
                  'RunMS58.MS51.MS45.90s.MS54.GTseq.val.GTseq.prod.merged')
geno.cols <- do.call(c, data_sets)
res <- lapply(data_sets, function(d){
  temp <- read.csv(paste0('results-raw/', d, '.20readsMin.num.genos.per.ind.csv')) |> 
    select(labID, genos)
  names(temp)[2] <- d
  return(temp)
}) |> 
  reduce(full_join, by = 'labID') |> 
  mutate(max.genos = pmax(!!!rlang::syms(geno.cols), na.rm = TRUE),
         MS58.max = ifelse(max.genos == RunMS58.replicates.merged,1,0),
         MS51.max = ifelse(max.genos == RunMS58.MS51.merged,1,0),
         MS45.max = ifelse(max.genos == RunMS58.MS51.MS45.90s.merged,1,0),
         MS54.max = ifelse(max.genos == RunMS58.MS51.MS45.90s.MS54.merged,1,0),
         GTseq.val.max = ifelse(max.genos == RunMS58.MS51.MS45.90s.MS54.GTseq.val.merged,1,0),
         GTseq.prod.max = ifelse(max.genos == RunMS58.MS51.MS45.90s.MS54.GTseq.val.GTseq.prod.merged,1,0))
res$run.to.use <- do.call(c, lapply(1:nrow(res), function(i){
  data_sets[[which(res[i,9:14] == 1)[1]]]
}))

write.csv(res, file = 'results-raw/genos.per.ind.across.runs.csv')
save(res, file = 'results-R/genos.per.ind.across.runs.rda')

all.tgts <- do.call(rbind, lapply(data_sets, function(d){
  readRDS(paste0('results-R/tgt.all.inds', d, '.rds')) |> 
    mutate(Indiv = paste0(Indiv, '_', d))
}))
LABIDs <- unique(all.tgts$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- all.tgts[grep(substr(r, start = 1, stop = 8), all.tgts$Indiv),]
  #  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches = compare.replicates(rep.tgt)
})) 
mismatches.to.check <- mismatches.to.check |> 
  mutate(labID = substr(mismatches.to.check$Indiv, start = 1, stop = 8)) |> 
  left_join(select(res, c(labID, run.to.use)), by = 'labID') |> 
  distinct()
if(nrow(mismatches.to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/all.tgts.genotype.mismatches.rda"))
  save(mismatches.to.check, file = paste0("results-R/all.tgts.genotype.mismatches.rda"))
}
write.csv(mismatches.to.check, file = paste0('results-raw/mismatches.all.runs.csv'))

runs.to.use <- read.csv('data-raw/sam.files/Runs.to.use.by.ind.csv') |> 
  select(c(labID, run.to.use))
for(i in 1:nrow(runs.to.use)){
  file.copy(from = paste0('data-raw/sam.files/', runs.to.use$run.to.use[i], '/', runs.to.use$labID[i], '.merged.sam'), 
              to = paste0('data-raw/sam.files/final.sams/', runs.to.use$labID[i], '.merged.sam'))
}
