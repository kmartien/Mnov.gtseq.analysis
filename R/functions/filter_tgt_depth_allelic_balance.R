source("R/functions/tgt.2.geno.table.R")
source("R/functions/haps.w.ns.R")

filter_tgt <- function(tidy.vcf, AB.min.het = 3/7, AB.max.homo = 2/8, 
                      min.read.depth = 10) {
  
  tgt <- tidy.vcf |> 
    mutate(AO = as.numeric(AO),
           gt_AO = as.numeric(gt_AO),
           read_depth = gt_RO + gt_AO,
           allelic_balance = gt_AO / gt_RO)

  # remove haplotypes that are NA
  tgt <- filter(tgt, gt_GT_alleles != '.')

  # remove genotypes that don't meet the minimum read depth criterion
  tgt <- filter(tgt, read_depth >= min.read.depth)

  # remove genotypes for individuals whose read depth of its major haplotype at a locus is less than min.read.depth/2
  tgt <- filter(tgt, gt_RO >= min.read.depth/2)
  
  # re-call genos based on allelic_balance thresholds
  tgt <- tgt |> 
    mutate(new.geno = '.',
      new.geno = ifelse(allelic_balance <= AB.max.homo, paste(REF, REF, sep= '/'), 
                        ifelse (allelic_balance >= (1/AB.max.homo), paste(ALT, ALT, sep='/'), 
                                ifelse(allelic_balance >= AB.min.het & allelic_balance <= 1/AB.min.het, paste(REF, ALT, sep='/'), '.')
                        )
      )
    ) |> 
    relocate(new.geno, .after = gt_GT_alleles) |> 
    mutate(gt_GT_alleles = new.geno) |> 
    select(-new.geno) |> 
    filter(gt_GT_alleles != '.')
  
  return(tgt)
}
