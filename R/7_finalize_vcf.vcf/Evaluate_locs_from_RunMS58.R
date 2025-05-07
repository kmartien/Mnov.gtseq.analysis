library(vcfR)
library(vcftoolsR)
library(tidyverse)
library(dplyr)
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC")
source("R/functions/mplot2tgt.R")
source("R/functions/Compare.replicates.R")
setwd("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.analysis")

project <- "RunMS58"

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 20
num.locs <- 368
min.genos.per.ind <- num.locs * 0.8

tgt <- mplot2tgt(project = project, AB.min.het = AB.min.het, AB.max.homo = AB.max.homo,
                 min.read.depth = min.read.depth)

# Mnov_gtseq_155 and Mnov_gtseq_139 are bad loci;they appear to have
# off-target amplicons 
tgt <- filter(tgt, locus != 'Mnov_gtseq_155' & locus != 'Mnov_gtseq_139')

# compare replicates
LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),]
#  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt)
}))
if(nrow(mismatches.to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/", project, ".genotype.mismatches.rda"))
  save(mismatches.to.check, file = paste0("results-R/", project, ".genotype.mismatches.rda"))
}

# Identify samples with Ns or Xs in their haplotypes or more than 2 haplotypes
questionable.hap <- sapply(1:nrow(tgt), function(i){
  ifelse(length(grep("N", tgt$gt[i],)) > 0 || length(grep("X", tgt$gt[i],)) > 0 
         || tgt$num.haps[i] > 2, TRUE, FALSE) 
})
genos.to.check <- filter(tgt, questionable.hap == TRUE)
if(nrow(genos.to.check > 0)) {
  print("Some samples with Ns or Xs in their haplotypes or more than 2 haplotypes")
  print(paste0("Questionable genotypes saved to results-R/", project, ".genos.to.check.rda"))
  save(genos.to.check, file = paste0("results-R/", project, ".genos.to.check.rda"))
}

# summarize individual data
missing.data.ind <- data.frame(table(tgt$Indiv[!is.na(tgt$gt)])) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))
num.inds <- nrow(missing.data.ind)

# summarize locus data
tgt_long <- tgt |> 
  select(locus, Indiv, gt, depth.1, depth.2) |> 
  separate_wider_delim(
    cols = gt, 
    delim = '/', 
    names = c('haplo.1', 'haplo.2'), 
    cols_remove = FALSE
  ) |> 
  mutate(depth.2 = ifelse(haplo.2 == haplo.1, NA, depth.2))
loc.sum <- tgt_long %>%
  pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
  mutate(tmp = strsplit(as.character(value), "")) %>%
  unnest(tmp) %>%
  group_by(locus, Indiv, hap) %>%
  mutate(name = 1:n()) %>%
  pivot_wider(id_cols = c(locus, Indiv, gt, hap), values_from = tmp, names_prefix = 'snp') |> 
  ungroup() |> 
  filter(!is.na(gt)) |> 
  group_by(locus) |> 
  summarise(
    inds.genoed = n() / 2,
    num.unique.genos = length(unique(gt)),
    num.alleles.pos1 = length(unique(snp1)),
    num.alleles.pos2 = length(unique(snp2)),
    num.alleles.pos3 = length(unique(snp3)),
    num.alleles.pos4 = length(unique(snp4)),
    num.alleles.pos5 = length(unique(snp5))
  )
write.csv(loc.sum, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.locus.summary.csv"))

locs <- unique(tgt_long$locus)[order(unique(tgt_long$locus))]
allele.freqs <- lapply(
  locs,
  function(loc){
    filter(tgt_long, locus == loc) |> 
      pivot_longer(cols = c(haplo.1, haplo.2), names_to = 'hap') |> 
      group_by(value) |> 
      summarise(freq = n())
  }
)
names(allele.freqs) <- locs

# this code is to look at specific loci
focal.loc <- 'Mnov_gtseq_56'
allele.freqs[[focal.loc]]
View(filter(tgt_long, locus == focal.loc))

# I went through all of the loci and eliminated ones that were questionable - multiple
# SNPs, but only 2 alleles, loci where one SNP had a MAC < 10, etc. I made notes
# in the .locus.summary.csv document generated ~20 lines above and changed the name
# to 'results_raw/RunMS58.locus.notes.csv. I then printed out all locus names
# from the vcf file used to generate the tgt object above, and deleted SNPs that
# I had deemed untrustworthy (see results_raw/RunMS58.locus.notes.csv for reason)
# The next lines remove all SNPs not deemed trustworthy re-do the microhaplot calling

vcf <- read.vcfR(paste0("vcf/", project, ".filtered.recode.vcf"), convertNA = T)
tidy.vcf <- vcfR2tidy(vcf, single_frame = TRUE, toss_INFO_column = FALSE,
                      info_fields = c("DP","RO", "AO"), format_fields = c("GT", "RO", "AO"))$dat %>%
  mutate(locus = paste(CHROM, POS, sep= "_")) %>%
  relocate(locus, .after = POS)
write.csv(
  tidy.vcf |> select(locus) |> distinct(), 
  file = 'data-raw/locs.to.keep.csv',
  row.names = FALSE
)

# loc.to.remove file should have one column = locus, one column = position, and 
# should be a text file. The loci in it are the ones to remove
locs.2.keep <- read.csv('data-raw/locs.to.keep.csv')
write.table(
  tidy.vcf |> 
    filter(!locus %in% locs.2.keep$locus) |> 
    select(c(CHROM, POS)) |> 
    distinct(), 
  file = 'data-raw/locs.to.remove.txt',
  col.names = FALSE,
  row.names = FALSE, 
  quote = FALSE
)
vcftools.rmPOS(vcf.fn = paste0("vcf/", project, ".filtered.recode"),
               res.fn = paste0("vcf/", project, ".locs.to.keep"),
               CHROM.POS.file.name = 'data-raw/locs.to.remove.txt')
  
# Use the vcf generated in the previous step for future microhaplot runs

# geno.table <- tgt.2.geno.table(tgt)
# geno.table$num.genos <- do.call(rbind, lapply(1:nrow(geno.table), function(i){
#   length(which(!is.na(geno.table[i,2:ncol(geno.table)])))
# }))
# 
# 
# save(geno.table, tgt, loc.sum, allele.freqs, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.geno.eval.rda"))
