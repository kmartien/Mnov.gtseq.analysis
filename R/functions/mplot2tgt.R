source("R/functions/tgt.2.geno.table.R")
source("R/functions/haps.w.ns.R")

mplot2tgt <- function(project, AB.min.het = 3/7, AB.max.homo = 2/8, 
                      min.read.depth = 10) {
  
  why.removed <- list()
  genos <- readRDS(paste0("/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot/", project, ".rds"))[,-1]
  genos <- mutate(genos, id.loc = paste0(id, "-", locus))
  genos$to.remove <- FALSE
  
  # remove haplotypes that are NA
  genos$to.remove[which(is.na(genos$haplo))] <- TRUE
  why.removed$na <- length(which(genos$to.remove))
  
  # remove genotypes that don't meet the minimum read depth criterion
  inds <- unique(genos$id)
  loci <- unique(genos$locus)
  num.locs <- length(loci)
  tot.depth <- data.frame(do.call(rbind, lapply(inds, function(i){
    df <- filter(genos, id == i) %>% filter(to.remove == FALSE)
    do.call(rbind, lapply(loci, function(l){
      haps <- filter(df, locus == l)
      data.frame(i, l, sum(haps$depth))
    }))
  })))
  names(tot.depth) <- c("id", "locus", "tot.depth")
  
  genos.to.drop <- filter(tot.depth, tot.depth < min.read.depth) %>% mutate(id.loc = paste0(id, "-", locus))
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$minDepth <- length(which(genos$to.remove))
  print("Done filtering on read depth")
  
  # remove genotypes for individuals whose read depth of its major haplotype at a locus is less than min.read.depth/2
  genos.to.drop <- filter(genos, depth < min.read.depth/2) %>% filter(rank == 1) %>% filter(to.remove == FALSE)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$majorAlleleDepth <- length(which(genos$to.remove))
  
  # remove haplotypes with rank > 1 and AB < AB.max.hom
  # this removes minor allele(s) from genos where the major allele frequency
  # is high enough to call a homozygote
  genos$to.remove[which(genos$rank > 1 & genos$allele.balance < AB.max.homo)] <- TRUE
  why.removed$AB.max.hom <- length(which(genos$to.remove))
  
  # remove genotypes with MAF < AB.min.het
  # this removes genotypes with minor allele frequency not high enough to call a heterozygote
  genos.to.drop <- filter(genos, rank > 1) %>% filter(to.remove == FALSE) %>% 
    filter(allele.balance < AB.min.het) %>% select(id.loc)
  genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
  why.removed$AB.min.het <- length(which(genos$to.remove))
  
  # remove genotypes with rank > 4
  genos$to.remove[which(genos$rank > 4)] <- TRUE
  why.removed$rank.gt.4 <- length(which(genos$to.remove))
  print("Done filtering on rank and allelic balance")
  
  # create tgt-like data structure
  tgt <- data.frame(do.call(rbind, lapply(unique(genos$id.loc), function(x){
    df <- filter(genos, id.loc == x)
    if (nrow(df) > 4) df <- df[1:4,]
    loc <- df$locus[1]
    ind <- df$id[1]
    haps <- df$haplo
    depth <- df$depth
    if(length(haps) < 4) {
      haps <- c(haps, rep(NA, (4-length(haps))))
      depth <- c(depth, rep(0, (4-length(depth))))
    }
    #create gt field
    if(df$to.remove[1]) gt <- NA else{
      if(nrow(df) == 1) gt <- paste(haps[1],haps[1], sep= "/") else {
        if(df$to.remove[2]) gt <- paste(haps[1],haps[1], sep= "/") else {
          gt <- ifelse(haps[2] < haps[1], paste(haps[2],haps[1], sep= "/"), paste(haps[1],haps[2], sep= "/"))
        }
      }
    }
    num.haps <- sum(!df$to.remove) 
    res <- c(loc, ind, gt, haps, depth, num.haps)
    names(res) <- c("locus", "Indiv", "gt", "haplo.1", "haplo.2", "haplo.3", "haplo.4", 
                    "depth.1", "depth.2", "depth.3", "depth.4", "num.haps")
    return(res)
  })))
  
  tgt$depth.1 <- as.integer(tgt$depth.1)
  tgt$depth.2 <- as.integer(tgt$depth.2)
  tgt$depth.3 <- as.integer(tgt$depth.3)
  tgt$depth.4 <- as.integer(tgt$depth.4)
#  return(list(tgt = tgt, why.removed = why.removed))
  return(tgt)
}
