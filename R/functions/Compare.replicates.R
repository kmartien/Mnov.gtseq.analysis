# Compares the genotypes of two or more replicate samples. Argument rep.tgt is
# the tgt object created by create.tgt.from.microhaplot.R, but with only the rows
#for the replicates being compared

compare.replicates <- function(rep.tgt){
#  rep.tgt <- unite(rep.tgt, gt, haplo.1, haplo.2, sep = '/', remove = FALSE)
  locs <- unique(rep.tgt$locus)
  mismatches <- sapply(locs, function(l){
    filter(rep.tgt, locus == l) %>% select(gt) %>% unique() %>% na.omit() %>% nrow()
  })
  filter(rep.tgt, locus %in% names(which(mismatches > 1))) %>% arrange(locus)

}