# Does minimal filtering of GTseq data before moving over to vcf2genotypes.R

source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/Summarize.vcf.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/Visualized.filtered.SNPs.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/basic.vcf.filtering.R")
library(vcftoolsR)

PROJECT <- "RunMS58"
fname <- PROJECT
vcf.dir <- "vcf"
results.dir <- "results-raw/vcf_filtering"

# Filter out low-confidence SNP calls
# basic.vcf.filtering defaults: minDP < 5, minQ < 20, meanDP < 15, mac < 3, remove monomorphic sites
res <- basic.vcf.filtering(vcf.dir, fname, paste0(fname,".basicFilters"))
filter.res <- res$filter.res
old.fname <- res$fname
fname <- paste0(fname,".gtseqFilters")
file.rename(from = paste0(old.fname, ".recode.vcf"), to = paste0("vcf/", fname, ".recode.vcf"))

#  Summarize again and plot individual and locus summaries
summarize.vcf(vcf.dir, results.dir, fname = paste0(fname, ".recode"), res.name = fname)

# Decompose variants and retain only SNPs.
#########################################################################
# Execute the next line and paste the result into a terminal window
paste0("vcfallelicprimitives vcf/", fname, ".recode.vcf --keep-info --keep-geno > vcf/", PROJECT, ".SNPs.vcf")
#########################################################################

fname <- paste0(PROJECT, ".SNPs")
filter.res$SNPs <- vcftools.removeIndels(paste0("vcf/", fname), paste0("vcf/", fname))

# Summarize again
summarize.vcf(vcf.dir, results.dir, fname = paste0(fname,".recode"), res.name = fname)

# Missing data iterative filtering
imiss <- read.table(paste(results.dir,"/",fname,".imiss",sep=""), header = TRUE, stringsAsFactors = FALSE)
lmiss <- read.table(paste(results.dir,"/",fname,".lmiss",sep=""), header = TRUE, stringsAsFactors = FALSE)
print(paste0("Before iterative filtering, max_missing per individual = ", max(imiss$F_MISS), 
             "; max_miss per locus = ", max(lmiss$F_MISS)))
LQ_indv <- imiss %>% filter(F_MISS > 0.75) %>% select(INDV)
write.table(LQ_indv, paste0(results.dir,"/LQ_Ind_", 75),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
filter.res$final <- vcftools.removeInd(paste0("vcf/", fname, ".recode"), paste0("vcf/", PROJECT, ".imiss"),
                                       ind.file.name = paste0(results.dir,"/LQ_Ind_",75))
fname <- paste0(PROJECT, ".imiss")

filter.res$lmiss75 <- vcftools.maxMiss(paste0("vcf/", fname, ".recode"),
                            paste0("vcf/", PROJECT, ".filtered"), max.miss = 0.75)

fname <- paste0(PROJECT, ".filtered")

summarize.vcf(vcf.dir, results.dir, fname = paste0(fname,".recode"), res.name = fname)
loc_stats_raw <- read.loc.stats(dir = results.dir, fname) 
snps.per.contig <- loc_stats_raw %>%
  dplyr::group_by(CHR) %>% dplyr::summarise(Count = dplyr::n())

save(filter.res, file = paste0("data/", fname, "new.rda"))
