library(tidyverse)
library(strataG)
library(dplyr)
library(Mnov.GTseq.data)
library(swfscMisc)
data("GTseq.samps.final")
data("Mnov.strata")
data("id.key")

project <- "all"
min.reads <- 20

Mnov.strata$LABID <- paste0("z0", zero.pad(as.numeric(Mnov.strata$LABID)))
all.samps <- select(GTseq.samps.final, c(AnimalID, LABID, SEX, HAP, mito.haps)) %>%
  left_join(select(Mnov.strata, c(LABID, wint.area, feed.area, herd, HI.v.HI.SEAK, MnMx.v.MnMx.CA.OR.WA)))

load(paste0("results-R/", project,".", min.reads, "readsMin.geno.eval.rda"))

genos <- filter(geno.table, num.genos >= 170) %>% select(-num.genos) %>%
  column_to_rownames(var = "Indiv")

num.g.per.loc <- do.call(c, lapply(1:ncol(genos), function(i){
  sum(!is.na(genos[, i]))
})) %>% data.frame()
num.g.per.loc$locus <- colnames(genos)
names(num.g.per.loc) <- c("num.genos", "locus")
num.g.per.loc <- filter(num.g.per.loc, num.genos > nrow(genos) * 0.5)

genos <- select(genos, all_of(num.g.per.loc$locus))

split.genos <- alleleSplit(genos, sep= "/") %>% 
  data.frame() %>%
  rownames_to_column(var = "LABID")

df <- right_join(all.samps, split.genos, by = "LABID")
missing.animalid <- filter(df, is.na(AnimalID)) %>% select(LABID)
df <- df[-which(is.na(df$AnimalID)),]

df.strata <- select(df, c(AnimalID, wint.area, feed.area, herd, HI.v.HI.SEAK, MnMx.v.MnMx.CA.OR.WA)) %>%
  column_to_rownames(var = "AnimalID")
loc.col <- grep("Mnov_gtseq", names(df))[1]

g <- df2gtypes(df, ploidy = 2, id.col = 1, 
               strata.col = NULL, loc.col = 9, schemes = df.strata)

save(g, df, file = "data/data.for.PSRG.2024.rda")
