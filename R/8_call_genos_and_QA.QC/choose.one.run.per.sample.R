library(dplyr)
library(tidyverse)
library(swfscMisc)
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/tgt.2.geno.table.R")

projects <- list("RunMS58", "GTseq.prod", "GTseq.val", "RunMS51", "RunMS54", "RunMS43", "RunMS45.90s")
tgt.list <- lapply(projects, function(p){
  load(paste0("results-R/", p, ".10readsMin.geno.eval.rda"))
  tgt$Indiv <- paste0(tgt$Indiv, "-", p)
  return(tgt)
})
names(tgt.list) <- projects

tgt <- do.call(rbind, tgt.list)

geno.table <- tgt.2.geno.table(tgt)

geno.table <- cbind(do.call(rbind, strsplit(geno.table$Indiv, split = "-")), geno.table)
names(geno.table)[1:2] <- c("LABID", "Run")
geno.table$LABID <- substr(geno.table$LABID, start = 1, stop = 8)
geno.table <- geno.table %>% mutate(Run = recode(Run, "RunMS58" = 1, "GTseq.prod" = 2, "GTseq.val" = 3, "RunMS51" = 4, "RunMS54" = 5, "RunMS45.90s" = 6, "RunMS43" = 7))

geno.table$num.genos <- sapply(1:nrow(geno.table), function(i){
  sum(!is.na(geno.table[i,-1]))
})

samps <- unique(geno.table$LABID)

best.runs <- do.call(c, lapply(samps, function(i){
  i.samps <- filter(geno.table, LABID == i)
  most.genos <- i.samps[which(i.samps$num.genos == max(i.samps$num.genos)),]
  most.genos <- filter(i.samps, num.genos == max(num.genos)) %>%
    filter(Run == min(Run)) %>% select(Indiv)
  return(most.genos$Indiv[])
}))

final.geno.table <- filter(geno.table, Indiv %in% best.runs)
