library(dplyr)
library(tidyverse)
library(swfscMisc)
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/tgt.2.geno.table.R")

#mismatched.samples <- c("z0196894.RunMS54", "z0107078.RunMS51", "z0197663.RunMS54", "z0196841.RunMS54", "z0053845.RunMS51")
samps.without.replicates <- read.csv("data-raw/samps.wo.replicates.csv")

projects <- list("GTseq.prod", "GTseq.val", "RunMS43", "RunMS45.90s", "RunMS51", "RunMS54", "RunMS58")
tgt.list <- lapply(projects, function(p){
  load(paste0("results-R/", p, ".10readsMin.geno.eval.rda"))
  tgt$Indiv <- paste0(tgt$Indiv, ".", p)
  return(tgt)
})
names(tgt.list) <- projects

tgt <- do.call(rbind, tgt.list)

geno.table <- tgt.2.geno.table(tgt)

#genos.to.find <- geno.table[which(geno.table$Indiv %in% mismatched.samples),]
genos.to.find <- geno.table[which(geno.table$Indiv %in% samps.without.replicates$LABID),]

#geno.table <- geno.table[-which(geno.table$Indiv %in% mismatched.samples),]
geno.table <- geno.table[-which(geno.table$Indiv %in% samps.without.replicates$LABID),]

num.geno.matches <- function(genos.to.find, geno.table){
  num.matches <- lapply(1:nrow(genos.to.find), function(x){
    target.geno <- genos.to.find[x,-1]
    diffs <- do.call(bind_rows, lapply(1:nrow(geno.table), function(y){
      id <- geno.table[y,1]
      gt <- geno.table[y,-1]
      overlap <- sum(!is.na(gt[1,which(!is.na(target.geno))]))
      matches <- sum(gt[1,] == target.geno, na.rm = TRUE)
      return(c(id = id, overlap = overlap, matches = matches))
    })) 
  })
  names(num.matches) <- genos.to.find$Indiv
  return(num.matches)
}

num.matches <- num.geno.matches(genos.to.find, geno.table)
max.matches <- do.call(rbind, lapply(num.matches, function(i){
  return(i[which(i$matches == max(i$matches)),])
}))

write.csv(max.matches, file = "results-raw/matches.samps.without.reps.csv")
save(num.matches, file = "results-R/Identifying.mixed.up.samples.rda")
