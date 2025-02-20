# Implements filtering scheme 5 from O'Leary et al 2018

source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/Summarize.vcf.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/Visualized.filtered.SNPs.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/functions/basic.vcf.filtering.R")
library(vcftoolsR)

PROJECT <- "RunMS58"
fname <- PROJECT
vcf.dir <- "vcf"
results.dir <- "results-raw/vcf_filtering"

### Summarize raw stats
summarize.vcf(vcf.dir, results.dir, fname, res.name = fname)

# Filter out low-confidence SNP calls
# basic.vcf.filtering defaults: minDP < 5, minQ < 20, meanDP < 15, mac < 3, remove monomorphic sites
res <- basic.vcf.filtering(vcf.dir, fname, paste0(fname,".basicFilters"), mac = 1)
#fname <- paste0(fname,".basicFilters")

filter.res <- res$filter.res
old.fname <- res$fname
fname <- paste0(fname,".basicFilters")
file.rename(from = paste0(old.fname, ".recode.vcf"), to = paste0("vcf/", fname, ".recode.vcf"))

#  Summarize again and plot individual and locus summaries
summarize.vcf(vcf.dir, results.dir, fname = paste0(fname, ".recode"), res.name = fname)

# Remove sites with excess reads
# Includes only loci with mean read depth across all individuals less than or 
# equal to the 95th percentile of the distribution (the value is noted on the locus summary plot). 

loc_stats_raw <- read.loc.stats(dir = results.dir, fname) 
names(loc_stats_raw) <- sapply(names(loc_stats_raw), function(x) strsplit(x,fname))
reads.95pctile <- sort(loc_stats_raw$MEAN_DEPTH_)[.95*length(loc_stats_raw$MEAN_DEPTH_)]

filter.res$excessReads <- vcftools.excessReads(paste0("vcf/", fname, ".recode"), 
                                  paste0("vcf/", fname,  ".excessReads"), reads.95pctile)
fname <- paste0(fname,".excessReads")

# Decompose variants and retain only SNPs.
#########################################################################
# Execute the next line and paste the result into a terminal window
paste0("vcfallelicprimitives vcf/", fname, ".recode.vcf --keep-info --keep-geno > vcf/", PROJECT, ".SNPs.vcf")
#########################################################################

fname <- paste0(PROJECT, ".SNPs")

filter.res$SNPs <- vcftools.removeIndels(paste0("vcf/", fname), paste0("vcf/", fname))

# Summarize again
summarize.vcf(vcf.dir, results.dir, fname = paste0(fname,".recode"), res.name = fname)

#########################################################################

GET RID OF EVERYTHING FROM HERE TO FINAL FILTERING BY MISSING DATA (NEXT SET OF HASTAGS)


# INFO filters
# Filter based on allelic balance (AB), strandednesss, mapping quality ratio, and quality/depth ratio
# For flashed data, also filter based on proper pairing

###  Execute these paste commands, paste each result into a terminal window,
###  delete the forward slashes (\) before the double quotes, and press enter to execute
vcffilter.AB <- paste0('vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.1" vcf/', fname, '.recode.vcf > vcf/', fname, '.AB.vcf')
vcffilter.AB
# Excluding strandedness filter since these are single-end data
#vcffilter.strand <- paste0('vcffilter -s -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" vcf/', fname, '.AB.vcf > vcf/', fname, '.AB.strand.vcf')
#vcffilter.strand
vcffilter.mapQual <- paste0('vcffilter -f "MQM / MQMR > 0.25 & MQM / MQMR < 1.75" vcf/', fname, '.AB.vcf > vcf/', fname, '.AB.mapQual.vcf')
vcffilter.mapQual

# High depth/quality ratio filtering (remove sites with depth > mean + 1 stddev and qual < 2*depth)

# calculate depth and quality per locus 
vcftools.siteDepth(paste0("vcf/", fname, ".AB.mapQual"),
                   paste0(results.dir, "/", fname, ".AB.mapQual"))
vcftools.siteQual(paste0("vcf/", fname, ".AB.mapQual"),
                   paste0(results.dir, "/", fname, ".AB.mapQual"))
fname <- paste0(fname, ".AB.mapQual")
site_qual <- read.table(paste(results.dir, "/",fname,".lqual",sep=""), header = TRUE, stringsAsFactors = FALSE)
meanDP <- read.table(paste(results.dir, "/",fname,".ldepth.mean",sep=""), header = TRUE, stringsAsFactors = FALSE)
meanDP$QUAL <- site_qual$QUAL
high_depth <- subset(meanDP, MEAN_DEPTH > (mean(MEAN_DEPTH) + sqrt(var(MEAN_DEPTH))))
loc.to.remove <- subset(high_depth, QUAL < 2*MEAN_DEPTH)
write.table(loc.to.remove, paste0(results.dir, '/loci.to.remove.txt'), col.names = FALSE,
            row.names = FALSE, quote = FALSE)

# If any loci have high mean depth and quality < 2 * mean depth, remove them
if(length(loc.to.remove$CHROM) > 0) {
  write.csv(loc.to.remove[,c(1,2)], file="loc.to.remove.csv", row.names=FALSE)
  vcftools.rmPOS(paste0("vcf/", fname), 
                   paste0("vcf/", fname, ".highDepth"),
                   paste0(results.dir, "/loci.to.remove.txt"))
  fname <- paste0(fname, ".highDepth")
}
##############################################################################

# Remove sites with more than 75% missing data and individuals with more than 75% missing
vcftools.imiss(paste0("vcf/", fname, '.recode'), paste0(results.dir, "/", fname))
imiss <- read.table(paste(results.dir, "/",fname,".imiss",sep=""), header = TRUE, stringsAsFactors = FALSE)
LQ_indv <- imiss %>% filter(F_MISS > 0.75) %>% select(INDV)
write.table(LQ_indv, paste0(results.dir, "/LQ_Ind_", 75),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
filter.res$imiss75 <- vcftools.removeInd(paste0("vcf/", fname, '.recode'), paste0("vcf/", fname, ".imiss75"),
                            ind.file.name = paste0(results.dir, "/LQ_Ind_",75))
fname <- paste0(fname, ".imiss75")
filter.res$lmiss75 <- vcftools.maxMiss(paste0("vcf/", fname, ".recode"), paste0("vcf/", PROJECT, ".final"), max.miss = 0.25)
fname <- paste0(PROJECT, ".final")

#3.1 Summarize final, filtered dataset
summarize.vcf(vcf.dir, results.dir, paste0(fname, ".recode"), res.name = fname)
imiss <- read.table(paste(results.dir, "/",fname,".imiss",sep=""), header = TRUE, stringsAsFactors = FALSE)
lmiss <- read.table(paste(results.dir, "/",fname,".lmiss",sep=""), header = TRUE, stringsAsFactors = FALSE)
print(paste0("After final filtering, max_missing per individual = ", max(imiss$F_MISS), 
             "; max_miss per locus = ", max(lmiss$F_MISS)))

save(filter.res, file = paste0("results-R/", PROJECT, ".filter.res.rda"))
