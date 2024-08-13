library(assignPOP)
library(strataG)
library(genepop)
library(tidyverse)
library(ggplot2)
source("R/functions/genepopWrite.KKM.R")
load("data/data.for.PSRG.2024.rda")

strat.scheme <- "herd"

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

strata2keep <- c("CentAm-CA.OR.WA", "MnMx-CA.OR.WA")

#only do pairwise tests on strata with 10 or more samples
g.westcoast <- g.stratified[,,which(getStrataNames(g.stratified) %in% strata2keep)]
#pws.struct <- pairwiseSummary(pairwiseTest(g.10, nrep = n.reps.pvals))

#g.pop.infile <- genepopWrite.KKM(g.westcoast, path = "data-raw", filename = "genepop.westcoast.txt")
g.pop.infile <- genepopWrite(g.westcoast)

genepop.herd <- read.Genepop(g.pop.infile$fname, pop.names = strata2keep, haploid = FALSE)
genepop.herd$SampleID <- gsub(" ", "_", genepop.herd$SampleID)

genepop.herd.Rd <- reduce.allele(genepop.herd, p = 0.95)

mdl <- "randomForest"
res.dir <- paste0("results-raw/assignPOP/MC.",mdl,"/")
# MCMC Assignment
assign.MC( genepop.herd.Rd, dir=res.dir, train.inds=c(0.7, 0.8, 0.9),
           train.loci=c(0.5,1), loci.sample="fst", iterations=50,
           model=mdl )

accuMC <- accuracy.MC(dir = res.dir)

accuracy.plot( accuMC, pop=c("all", strata2keep)) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=1/length(strata2keep),yend=1/length(strata2keep),colour="red",size=1) +
  #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
ggtitle("Monte-Carlo cross-validation using genetic loci")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))

# K-fold assignment
assign.kfold( genepop.herd.Rd, k.fold=c(3,4,5), train.loci=c(0.8, 0.9, 1),
              loci.sample="fst", dir=paste0("results-raw/assignPOP/westcoast.kfold.lda/"), model="lda" )

accuMC.kfold.switch <- accuracy.MC(dir = "results-raw/assignPOP/westcoast.kfold.lda/")
names(accuMC.kfold.switch)[1:2] <- c("train.loci", "KF")
accuMC.kfold.switch <- accuMC.kfold.switch %>% relocate(KF, .before = train.loci)
accuracy.plot( accuMC.kfold.switch, pop=c("all", strata2keep)) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=1/length(strata2keep),yend=1/length(strata2keep),colour="red",size=1) +
#Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))

