library(tidyverse)
library(Mnov.GTseq.data)
library(randomForest)
library(strataG)

load('data/gtypes_final.sams.no.dupes_minReads.20.rda')

g_mat_complete <- imdo(as.matrix(g), plot = FALSE)$opt.mat |> data_frame()
schemes_complete <- 
  
locs.w.missing.data <- which(numMissing(g)$num.missing>0)
temp <- numMissing(g)

# This is how I did it when I was designing the panel, but eliminating all loci
# with missing data doesn't leave me with much
g_complete <- g[,-locs.w.missing.data,]

g_rfDF <- gtypes2rfDF(g_complete)

freq <- table(g_rfDF$stratum)
sampsize <- rep(ceiling(min(freq / 2)), length(freq))

#convert predictors to a matrix so as not to exceed memory limits
x <- sapply(g_rfDF[-1], function(l) as.numeric(l))
rf <- randomForest(
  x,
  g_rfDF$stratum,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = ntree,
  keep.forest = FALSE
)
plot(rf)
