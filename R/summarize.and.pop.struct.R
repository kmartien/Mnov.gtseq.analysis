library(strataG)
library(genepop)
library(tidyverse)
library(ggplot2)
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")

n.reps.pvals <- 1000

strat.scheme <- "herd"

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)
ind.smry <- summarizeInds(g.stratified)

plot <- ggplot(ind.smry, aes(x = num.loci.missing.genotypes, y = pct.loci.homozygous)) +
  geom_point()
plot

#only do pairwise tests on strata with 10 or more samples
g.9 <- g.stratified[,,which(getNumInd(g.stratified, by.strata = TRUE)$num.ind >= 9)]
pws.struct <- pairwiseTest(g.9, nrep = n.reps.pvals)
pws.sum <- pairwiseSummary(pws.struct)
chi2.mat <- pairwiseMatrix(pws.struct, stat = 'CHIsq')
fst.mat <- pairwiseMatrix(pws.struct, stat = 'Fst') 

write.csv(chi2.mat, file = paste0("results-raw/pairwise_chi2_", strat.scheme, ".csv"))
write.csv(fst.mat, file = paste0("results-raw/pairwise_fst_", strat.scheme, ".csv"))
save(g.stratified, ind.smry, pws.struct, file = paste0("results-R/pop.struct.", strat.scheme, ".rda"))
