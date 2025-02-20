library(strataG)
library(genepop)
library(tidyverse)
library(ggplot2)
load("data/gtypes_all_minReads20.rda")

#n.reps.pvals <- 10

strat.scheme <- "wint.area"

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)
ind.smry <- summarizeInds(g.stratified)

plot <- ggplot(ind.smry, aes(x = num.loci.missing.genotypes, y = pct.loci.homozygous)) +
  geom_point()
plot

#only do pairwise tests on strata with 10 or more samples
g.10 <- g.stratified[,,which(getNumInd(g.stratified, by.strata = TRUE)$num.ind >= 10)]
#pws.struct <- pairwiseSummary(pairwiseTest(g.10, nrep = n.reps.pvals))

g.pop.infile <- genepopWrite(g.10)

g.pop.outfile <- Fst(g.pop.infile$fname, pairs = TRUE)


genepop.res <- read.csv(paste0("results-raw/genepop.results.", strat.scheme, ".csv"))
genepop.mat <- pivot_wider(select(genepop.res, c(Pop1, Pop2, p)), names_from = Pop2, values_from = p, names_sort = TRUE)
write.csv(genepop.mat, file = paste0("results-raw/genepop.", strat.scheme, ".res.mat.csv"))
save(g.stratified, ind.smry, genepop.res, file = paste0("results-R/pop.struct.", strat.scheme, ".rda"))
