library(tidyverse)
library(strataG)
library(dplyr)
#library(Mnov.GTseq.data)
#library(swfscMisc)

load(paste0("data/data.for.PSRG.2024.rda"))

strat.scheme <- "herd"
strata.col <- which(names(df) == strat.scheme)
loc.col <- grep("Mnov_gtseq", names(df))[1]

g <- df2gtypes(df, ploidy = 2, id.col = 1, strata.col = 8, loc.col = 9)
g.list <- strataSplit(g)

strats.to.analyze <- g.list[c(2,6,9)]

hwe.list <- lapply(1:length(strats.to.analyze), function(s){
  x <- hweTest(strats.to.analyze[[s]]) %>% data.frame() %>% rownames_to_column()
  names(x) <- c("locus", names(strats.to.analyze)[s])
  return(x)
})
hwe.res <- full_join(hwe.list[[1]], hwe.list[[2]], by = "locus")
hwe.res <- full_join(hwe.res, hwe.list[[3]], by = "locus")

hwe.res$num.significant <- do.call(rbind, lapply(1:nrow(hwe.res), function(i){
  sum(hwe.res[i,] < 0.05)
}))

write.csv(hwe.res, file = "data-raw/QA.QC/hwe.results.csv")

locs2exclude <- filter(hwe.res, num.significant > 1) %>% select(locus)

save(locs2exclude, hwe.res, file = "results-R/HWE.test.results.rda")

df <- select(df, -c(paste0(locs2exclude$locus,".1"),paste0(locs2exclude$locus,".2")))
save(df, file = "data/data.for.PSRG.2024.rda")

ld.list <- lapply(1:length(strats.to.analyze), function(s){
  x <- LDgenepop(strats.to.analyze[[s]], iterations = 10)# %>% data.frame() %>% rownames_to_column()
  #names(x) <- c("locus", names(strats.to.analyze)[s])
  return(x)
})
