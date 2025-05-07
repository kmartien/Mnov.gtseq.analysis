library(tidyverse)
library(strataG)
library(Mnov.GTseq.data)
library(swfscMisc)

data("id.key")
id.key$LABID <- paste0("z0", zero.pad(id.key$LABID))

project <- "final.sams.no.dupes"
min.reads <- 20

load(file = paste0("data/gtypes_", project, "_minReads.", min.reads, ".rda"))
load(file = paste0("results-R/", project, ".", min.reads, "readsMin.geno.eval.rda"))

dupes <- dupGenotypes(g)
dupe.mat <- pivot_wider(select(dupes, c(ids.1, ids.2, prop.loci.shared)), names_from = ids.2, values_from = prop.loci.shared)

dupe.genos <- lapply(1:nrow(dupes), function(i){
  filter(geno.table, Indiv == dupes$ids.1[i] | Indiv == dupes$ids.2[i])
})

ind.summary <- summarizeInds(g)
possible.dupe.inds <- mutate(ind.summary,
                      dupe.identified = id %in% dupes$ids.1 | id %in% dupes$ids.2) |> 
  filter(dupe.identified)

dupes.to.merge <- read.csv('data-raw/QA.QC/dupes.to.merge.csv')
ids.to.merge <- do.call(rbind, lapply(1:nrow(dupes.to.merge), function(i){
  filter(id.key, LABID %in% c(dupes.to.merge$ids.1[i], dupes.to.merge$ids.2[i])) |> cbind(i)
}))
write.csv(dupes.to.merge, 
          file = '/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.GTseq.data/data-raw/ids.to.merge.csv',
          row.names = FALSE)
for (i in 1:nrow(dupes.to.merge)){
  err <- system(paste0("samtools merge -o data-raw/sam.files/final.sams/", dupes.to.merge$ids.1[i], "B.merged.sam data-raw/sam.files/final.sams/", dupes.to.merge$ids.1[i], ".merged.sam data-raw/sam.files/final.sams/", dupes.to.merge$ids.2[i], ".merged.sam"))
  if(err == 0){
    ids.to.merge[[i]] <- filter(id.key, LABID %in% c(dupes.to.merge$ids.1[i], dupes.to.merge$ids.2[i]))
    file.remove(c(paste0("data-raw/sam.files/final.sams/", dupes.to.merge$ids.1[i], ".merged.sam"), paste0("data-raw/sam.files/final.sams/", dupes.to.merge$ids.2[i], ".merged.sam")))
    file.rename(from = paste0("data-raw/sam.files/final.sams/", dupes.to.merge$ids.1[i], "B.merged.sam"), to = paste0("data-raw/sam.files/final.sams/", dupes.to.merge$ids.1[i], ".merged.sam"))
  }
}

# update id.key as appropriate!!!!! - DONE 01.29.25
