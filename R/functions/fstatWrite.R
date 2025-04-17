fstatWrite <- function (g, label = NULL) 
{
  if (getPloidy(g) != 2) 
    stop("'g' must be a diploid object")
  .convertAlleles <- function(x) {
    x <- as.numeric(factor(x))
    x[is.na(x)] <- 0
    max.width <- formatC(x, width = max(2, nchar(x)), flag = "0")
  }
  genepop.fmt <- g@data %>% dplyr::group_by(.data$locus) %>% 
    dplyr::mutate(allele = .convertAlleles(.data$allele)) %>% 
    dplyr::ungroup() %>% dplyr::group_by(.data$stratum, .data$id, 
                                         .data$locus) %>% dplyr::summarize(genotype = paste(.data$allele, 
                                                                                            collapse = "")) %>% dplyr::ungroup() %>% tidyr::spread(.data$locus, 
                                                                                                                                                   .data$genotype) %>% dplyr::mutate(stratum = gsub(",", 
                                                                                                                                                                                                    "_", .data$stratum), id = gsub(",", "_", .data$id)) %>% 
    dplyr::ungroup()
  locus.names <- stats::setNames(colnames(genepop.fmt)[-(1:2)], 
                                 paste("LOC", 1:(ncol(genepop.fmt) - 2), sep = ""))
  fname <- paste0(label, "_Fstat_data.txt")
  fname <- gsub(" ", "_", fname)
  
  write(
    paste(
      length(getStrataNames(g)), 
      length(locus.names), 
      max(summarizeLoci(g)$num.alleles), 
      2, 
      sep= " "
    ),
    file = fname)
  for(i in 1:length(locus.names)){
    write(locus.names[i], append = TRUE, file = fname)
  }
  genepop.fmt <- split(genepop.fmt, genepop.fmt$stratum)
  for (p in 1:length(genepop.fmt)) {
    pop <- genepop.fmt[[p]]
    for (i in 1:nrow(pop)) {
      stratum <- p
      genotypes <- paste(pop[i, -(1:2)], collapse = " ")
      write(paste(stratum, genotypes, sep = " , "), file = fname, 
            append = TRUE)
    }
  }
  invisible(list(fname = fname, locus.names = locus.names))
}