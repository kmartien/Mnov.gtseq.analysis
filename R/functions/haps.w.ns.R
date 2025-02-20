haps.w.ns <- function(tgt){
  
  tgt$contains.ns <- sapply(1:nrow(tgt), function(i){
    ifelse(length(grep("N", tgt$haplo.x[i],)) > 0 || length(grep("N", tgt$haplo.y[i],)) > 0, TRUE, FALSE) 
  })
  
  length(which(tgt$contains.ns == TRUE))
  locs.w.ns <- filter(tgt, contains.ns == TRUE)
  inds.w.ns <- data.frame(table(locs.w.ns$Indiv))
  names(inds.w.ns)[1] <- "Indiv"
  locs.w.ns <- data.frame(table(locs.w.ns$locus))
  names(locs.w.ns)[1] <- "locus"
  
  return(list(inds.w.ns = inds.w.ns, locs.w.ns = locs.w.ns))
}
