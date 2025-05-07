library(vcfR)
library(vcftoolsR)
library(tidyverse)
library(dplyr)
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC")
source("R/functions/mplot2tgt.R")
source("R/functions/Compare.replicates.R")
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.analysis")

project <- "final.sams"

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 20
num.locs <- 322
min.genos.per.ind <- num.locs * 0.6

tgt <- mplot2tgt(project = project, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo,
                 min.read.depth = min.read.depth)
project <- 'final.sams.no.dupes'
saveRDS(tgt, file = paste0('results-R/tgt.no.dupes.', project, '.rds'))

# compare replicates
LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
#  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt)
}))
if(nrow(mismatches.to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/", project, ".genotype.mismatches.rda"))
  save(mismatches.to.check, file = paste0("results-R/", project, ".genotype.mismatches.rda"))
}

# Identify samples with Ns or Xs in their haplotypes or more than 2 haplotypes
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)

############################################################################
###### NEED TO DECIDE HOW TO DEAL WITH QUESTIONABLE HAPLOTYPES AT THIS POINT

# Remove one locus with an indel that's screwing up genotypes and two samples
# with evidence of contamination
tgt <- filter(tgt, locus != 'Mnov_gtseq_526') |> 
  filter(Indiv != 'z0037611' & Indiv != 'z0053845' & Indiv != 'z0195647')

# Change genotypes for questionable haps
genos_to_change <- read.csv('data-raw/genos_to_change.csv')
for (i in 1:nrow(genos_to_change)){
  idx <- which(tgt$locus == genos_to_change$locus[i] & tgt$Indiv == genos_to_change$Indiv[i])
  tgt$gt[idx] <- genos_to_change$gt[i]
}


############################################################################

# summarize individual data
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
inds.2.keep <- filter(missing.data.ind, genos >= min.genos.per.ind) |> 
  pull(labID)
tgt <- filter(tgt, Indiv %in% inds.2.keep)
num.inds <- length(unique(tgt$Indiv))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))

# summarize locus data
tgt_long <- tgt |> 
  select(locus, Indiv, gt, depth.1, depth.2) |> 
  separate_wider_delim(
    cols = gt, 
    delim = '/', 
    names = c('haplo.1', 'haplo.2'), 
    cols_remove = FALSE
  ) |> 
  mutate(depth.2 = ifelse(haplo.2 == haplo.1, NA, depth.2))
loc.sum <- tgt_long %>%
  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
  mutate(tmp = strsplit(as.character(value), "")) %>%
  unnest(tmp) %>%
  group_by(locus, Indiv, hap) %>%
  mutate(name = 1:n()) %>%
  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
  ungroup() |> 
  filter(!is.na(gt)) |> 
  group_by(locus) |> 
  summarise(
    inds.genoed = n() / 2,
    num.unique.genos = length(unique(gt)),
    num.alleles.pos1 = length(unique(snp1)),
    num.alleles.pos2 = length(unique(snp2)),
    num.alleles.pos3 = length(unique(snp3)),
    num.alleles.pos4 = length(unique(snp4)),
    num.alleles.pos5 = length(unique(snp5))
  )
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.csv"))
tgt <- filter(tgt,
              locus %in% 
                (filter(loc.sum,
                       inds.genoed >= num.inds * 0.5) |> 
                pull(locus)))
geno.table <- tgt.2.geno.table(tgt) 

save(geno.table, tgt, loc.sum, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.geno.eval.rda"))
