library(tidyverse)
library(strataG)
library(dplyr)
#library(Mnov.GTseq.data)
#library(swfscMisc)

project <- 'final.sams.no.dupes'
min.reads <- 20

load(paste0('data/gtypes_', project, '_minReads.', min.reads, '.rda'))

### REMOVE INDIVIDUALS THAT ARE OUTLILERS IN TERMS OF HIGH HOMOZYGOSITY ####
ind.summary <- summarizeInds(g)
high.homo.samps <- filter(ind.summary, pct.loci.homozygous > 0.7) |> 
  pull(id)
g <- g[-which(getIndNames(g) %in% high.homo.samps),,]

### HARDY-WEINBERG EQUILIBRIUM ###############################
herd.g <- stratify(g, 'herd')

hwe.list <- lapply(c("CentAm-CA.OR.WA", "MnMx-CA.OR.WA", "HI-SEAK.NBC"), function(s){
    x <- hweTest(herd.g[,,s]) %>% data.frame() %>% rownames_to_column()
    names(x) <- c("locus", s)
  return(x)
})
hwe.res <- hwe.list |> reduce(full_join, by = 'locus') |> 
  rowwise() %>%
  mutate(
    num.sig = sum(
      c_across('CentAm-CA.OR.WA':'MnMx-CA.OR.WA') < 0.05
    ),
    num.sig.after.correction = sum(
      c_across('CentAm-CA.OR.WA':'MnMx-CA.OR.WA') < (0.05/3)
    )
  ) %>%
  ungroup()

write.csv(hwe.res, file = paste0("data-raw/QA.QC/",project, ".hwe.results.csv"))

# Identify and exclude loci that were significantly out of HWE in one or more herd, after Bonferroni correction
locs2exclude <- filter(hwe.res, num.sig.after.correction > 0) %>% pull(locus)

g <- g[,-which(getLociNames(g) %in% locs2exclude),]
herd.g <- stratify(g, 'herd')
save(g, file = paste0('data/gtypes_', project, '_minReads.', min.reads, '.rda'))

### SUMMARIZE REMAINING LOCI ################################################

loc.sum <- summarizeLoci(g)
loc.sums.by.strat <- summarizeLoci(herd.g, by.strata = TRUE) |> 
  filter(stratum %in% c("CentAm-CA.OR.WA", "MnMx-CA.OR.WA", "HI-SEAK.NBC")) |> 
  mutate(obs_minus_exp = exptd.het - obsvd.het)

loc.sum <- left_join(loc.sum, 
g@data |> select(c(locus, allele)) |> 
  distinct() |> 
  filter(!is.na(allele)) |> 
  mutate(num.snps = nchar(allele)) |> 
  select(c(locus, num.snps)) |> 
  distinct()
)
write.csv(loc.sum, file = 'results-raw/final.loc.sum.csv')
save(loc.sums.by.strat, loc.sum, file = 'data/final.loc.sum.rda')

### SUMMARIZE REMAINING INDIVIDUALS #######################################

ind.sum <- summarizeInds(g)
samp.size.by.stratum <- table(paste(df.strata$wint.area, df.strata$feed.area, sep='_|_'))
write.csv(samp.size.by.stratum, file = 'results-raw/samp.size.by.herd.csv')

### LINKAGE DISEQUILIBRIUM ################################################

# Not doing LD. I know these loci aren't physically linked, so result would not be useful
# ld.overall <- LDgenepop(g)
# 
# ld.list <- lapply(1:length(strats.to.analyze), function(s){
#   print('next')
#   x <- LDgenepop(strats.to.analyze[[s]])# %>% data.frame() %>% rownames_to_column()
#   #names(x) <- c("locus", names(strats.to.analyze)[s])
#   return(x)
# })
# names(ld.list) <- names(strats.to.analyze)
# 
# ld.sig.res <- lapply(ld.list, function(s){
#   filter(s, p.value < 0.05) |> 
#     select(c(Locus.1, Locus.2, p.value))
# }) 
# for (i in 1:length(ld.sig.res)){
#   names(ld.sig.res[[i]])[3] <- paste0('p.val.',names(ld.list[i]))
# }
#  
# ld.sig.res <- ld.sig.res |> reduce(full_join) |> 
#   rowwise() %>%
#   mutate(
#     num.sig = sum(
#       c_across('p.val.CentAm-CA.OR.WA':'p.val.MnMx-CA.OR.WA') < 0.05, na.rm = TRUE
#     ),
#     num.sig.after.correction = sum(
#       c_across('p.val.CentAm-CA.OR.WA':'p.val.MnMx-CA.OR.WA') < (0.05/3), na.rm = TRUE
#     )
#   ) %>%
#   ungroup()
# 
# ld.list$all <- LDgenepop(g)
