library(seqinr)
library(vcfR)
library(tidyr)
library(dplyr)
library(tidyverse)
library(pegas)
source('R/functions/filter_tgt_depth_allelic_balance.R')
source("R/functions/Compare.replicates.R")
#source('R/functions/mplot2tgt.R')

run.label <- "RunMS58"
vcf <- read.vcfR(paste0("vcf/",run.label, ".filtered.recode.vcf"), convertNA = T)

# convert vcf to a tibble and add a column 'locus' that combines CHROM and POS
# of the SNPs
tidy.vcf <- vcfR2tidy(vcf, single_frame = TRUE, toss_INFO_column = FALSE,
                      info_fields = c("DP","RO", "AO"), format_fields = c("GT", "RO", "AO"))$dat %>%
  mutate(locus = paste(CHROM, POS, sep= "_")) %>%
  relocate(locus, .after = POS)

# filter out genotypes that don't meet minimum read depth and allelic balance thresholds
tgt <- filter_tgt(tidy.vcf, min.read.depth = 10)

num.locs <- length(unique(tgt$locus))
num.inds <- length(unique(tgt$Indiv))
imiss <- tgt |> 
  group_by(Indiv) |> 
  summarise(num.genos = n()) |> 
  mutate(prop.missing = 1 - (num.genos / num.locs))
lmiss <- tgt |> 
  group_by(locus) |> 
  summarise(num.inds.genod = n()) |> 
  mutate(prop.missing = 1 - (num.inds.genod / num.inds))

# identify loci that three genotypes (both homozygotes and the heterozygote)
loci.to.keep <- tgt |> 
  group_by(locus) |> 
  distinct(gt_GT_alleles) |> 
  summarise(num.genotypes = n()) |> 
  filter(num.genotypes > 2)
tgt <- filter(tgt, locus %in% loci.to.keep$locus)

# compare replicates
LABIDs <- unique(tgt$Indiv) %>% substr(start = 1, stop = 8)
replicates <- LABIDs[duplicated(LABIDs)]
mismatches.to.check <-do.call('rbind',lapply(replicates, function(r){
  rep.tgt <- tgt[grep(substr(r, start = 1, stop = 8), tgt$Indiv),] |> 
    rename(gt = gt_GT_alleles)
  #  rep.tgt <- filter(tgt, Indiv %in% c(r,paste0(r,"b")))
  mismatches <- compare.replicates(rep.tgt) |> mutate(LabID = r)
}))
if(nrow(mismatches.to.check > 0)) {
  print("Some replicates have mismatched genotypes")
  print(paste0("Mismatches saved to results-R/", run.label, ".genotype.mismatches.rda"))
  save(mismatches.to.check, file = paste0("results-R/", run.label, ".genotype.mismatches.rda"))
}
write.csv(mismatches.to.check, file = paste0('results-raw/', run.label, 'replicate.mismatches.csv'))
Need to compare replicates, then combine them, then filter by missingness