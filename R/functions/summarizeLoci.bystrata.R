summarizeLoci.bystrata <- function(g, by.strata = TRUE){
  by.cols <- if (by.strata) 
    c("stratum", "locus")
  else "locus"
  smry <- numGenotyped(g, by.strata) %>% dplyr::left_join(numMissing(g, 
          by.strata), by = by.cols) %>% dplyr::mutate(prop.genotyped = .data$num.genotyped/(.data$num.genotyped + 
          .data$num.missing)) %>% dplyr::left_join(numAlleles(g, 
          by.strata), by = by.cols) %>% dplyr::left_join(allelicRichness(g, 
          by.strata), by = by.cols) %>% dplyr::left_join(propUniqueAlleles(g, 
          by.strata), by = by.cols) %>% dplyr::left_join(heterozygosity(g, 
          by.strata, "expected"), by = by.cols) %>% dplyr::left_join(heterozygosity(g, 
          by.strata, "observed"), by = by.cols)
  if (getPloidy(g) == 1) {
    smry <- smry %>% dplyr::rename(num.haplotypes = .data$num.alleles, 
                                 num.unique.haplotypes = .data$num.unique, 
                                 haplotypic.diversity = .data$exptd.het.x) %>% dplyr::select(-.data$exptd.het.y, -num.genotyped.y)
  }
  smry
}
