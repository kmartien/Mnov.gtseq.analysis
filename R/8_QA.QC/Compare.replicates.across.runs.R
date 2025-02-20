library(dplyr)
library(tidyverse)
library(swfscMisc)
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/Compare.replicates.R")

projects <- list("GTseq.prod", "GTseq.val", "RunMS43", "RunMS45.90s", "RunMS51", "RunMS54", "RunMS58")
tgt.list <- lapply(projects, function(p){
  load(paste0("results-R/", p, ".10readsMin.geno.eval.rda"))
  tgt$Indiv <- paste0(tgt$Indiv, ".", p)
  return(tgt)
})
names(tgt.list) <- projects

tgt <- do.call(rbind, tgt.list)

unique.files <- unique(tgt$Indiv)

LABIDs <- unique.files %>% substr(start = 1, stop = 8)
replicates <- unique(LABIDs[duplicated(LABIDs)])

samp.check <-do.call('rbind',lapply(replicates, function(r){
  print(r)
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  unique.ids <- unique(rep.tgt$Indiv)
  locs <- unique(rep.tgt$locus)
  pairs.df <- data.frame(t(combn(unique.ids, 2)))

  do.call(rbind, lapply(1:nrow(pairs.df), function(p){
    overlapping.genos <- mismatches <- 0
    genos <- filter(rep.tgt, Indiv %in% pairs.df[p,1:2])
    for (l in 1:length(locs)){
      gts <- filter(genos, locus == locs[l]) %>% select(gt) %>% na.omit()
      if(nrow(gts) == 2) overlapping.genos <- overlapping.genos+1
      if(nrow(unique(gts)) > 1) mismatches <- mismatches+1
    }
    return(c(r, pairs.df[p,], overlapping.genos, mismatches))
  }))
})) %>% data.frame()
names(samp.check) <- c("LABID", "Rep1", "Rep2", "overlapping.genos", "mismatches")
samp.check$pct.mismatched <- samp.check$mismatches/samp.check$overlapping.genos

samps.with.mismatches <- filter(samp.check, mismatches > 0) %>% select(LABID) %>% unique()
samp.check <- filter(samp.check, LABID %in% samps.with.mismatches$LABID)

write.csv(samp.check, file = paste0("results-raw/all.runs.comparisons.csv"))

save(samp.check, file = "results-R/num.mismatches.btwn.runs.rda")
mismatches <- do.call('rbind',lapply(samps.with.mismatches$LABID, function(r){
  print(r)
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  compare.replicates(rep.tgt)
}))

write.csv(mismatches, file = paste0("results-raw/allruns.mismatches.csv"))
save(samp.check, mismatches, file = "results-R/mismatches.btwn.runs.rda")

samps.with.mismatches <- filter(samp.check, pct.mismatched > 0.1) %>% select(LABID) %>% unique()
samp.check.gt.10pct <- filter(samp.check, LABID %in% samps.with.mismatches$LABID)

write.csv(samp.check.gt.10pct, file = paste0("results-raw/allruns.gt.10pct.comparisons.csv"))

mismatches <- do.call('rbind',lapply(samps.with.mismatches$LABID, function(r){
  print(r)
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
  compare.replicates(rep.tgt)
}))

write.csv(mismatches, file = paste0("results-raw/all.runs.gt.10pct.mismatches.csv"))

