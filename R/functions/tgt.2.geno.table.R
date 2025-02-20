# Creates a data frame where rows are IDs and columns are genotypes.

tgt.2.geno.table <- function(tgt){
  tmp <- select(tgt, c(Indiv, locus, gt)) %>%
    pivot_wider(names_from = locus, values_from = gt, 
                names_sort = TRUE, names_vary = "slowest")
  return(tmp)
}