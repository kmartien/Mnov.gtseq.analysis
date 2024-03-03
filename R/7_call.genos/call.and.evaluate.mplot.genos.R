library(tidyverse)
library(dplyr)
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC")
source("R/functions/mplot2tgt.R")
source("R/functions/Compare.replicates.R")
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.analysis")

project <- "RunMS58"

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 10
num.locs <- 368
min.genos.per.ind <- num.locs * 0.8

tgt <- mplot2tgt(project = project, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo,
                 min.read.depth = min.read.depth)
project <- "RunMS58.allSNPs"
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
if(nrow(genos.to.check > 0)) {
  print("Some samples with Ns or Xs in their haplotypes or more than 2 haplotypes")
  print(paste0("Questionable genotypes saved to results-R/", project, ".genos.to.check.rda"))
  save(genos.to.check, file = paste0("results-R/", project, ".genos.to.check.rda"))
}

# summarize individual data
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))
num.inds <- nrow(missing.data.ind)

# summarize locus data
geno.table <- tgt.2.geno.table(tgt)
geno.table$num.genos <- do.call(rbind, lapply(1:nrow(geno.table), function(i){
  length(which(!is.na(geno.table[i,2:ncol(geno.table)])))
}))

loc.sum <- data.frame(colnames(geno.table)[-c(1,ncol(geno.table))], do.call(rbind, lapply(2:(ncol(geno.table)-1), function(l){
  inds.genoed <- length(which(!is.na(geno.table[1:nrow(geno.table),l])))
  num.haps <- sum(!is.na(unique(geno.table[,l])))
  c(inds.genoed, num.haps)
})))
names(loc.sum) <- c("locus", "genos", "num.haps")
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.csv"))

save(geno.table, tgt, loc.sum, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.geno.eval.rda"))

