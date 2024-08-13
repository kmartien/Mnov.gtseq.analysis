function (dir = NULL) 
{
  fileName_vec <- list.files(path = dir, pattern = "Out_*")
  fileName_vec <- sort(fileName_vec)
  noFiles <- length(fileName_vec)
  result01 <- read.table(paste0(dir, fileName_vec[1]), header = T, 
                         check.names = F, stringsAsFactors = T)
  pops <- names(result01)[4:length(names(result01))]
  noPops <- length(pops)
  Var1 <- NULL
  Var2 <- NULL
  train.inds <- NULL
  train.loci <- NULL
  iters <- NULL
  assign.rate.all <- NULL
  assign.rate.each <- as.data.frame(matrix(nrow = 0, ncol = noPops), 
                                    stringsAsFactors = F)
  for (i in 1:noFiles) {
    oneFileName <- unlist(strsplit(fileName_vec[i], split = "_"))
    train.inds[i] <- oneFileName[2]
    train.loci[i] <- oneFileName[3]
    iters[i] <- unlist(strsplit(oneFileName[4], split = ".txt"))
    df <- read.table(paste0(dir, fileName_vec[i]), header = T, 
                     stringsAsFactors = T)
    df$origin.pop <- factor(df$origin.pop, levels = levels(factor(pops)))
    df$pred.pop <- factor(df$pred.pop, levels = levels(factor(pops)))
    ctable <- table(df$origin.pop, df$pred.pop)
    ftable <- as.data.frame(ctable)
    totalSample <- sum(ftable$Freq)
    AllcorrectNo <- sum(subset(ftable, Var1 == Var2)$Freq)
    assign.rate.all[i] <- AllcorrectNo/totalSample
    popCorrectRate_vec <- NULL
    for (p in pops) {
      pop_size <- sum(subset(ftable, Var1 == p)$Freq)
      if (pop_size == 0) {
        popCorrectRate = 0
      }
      else {
        popCorrectNo <- subset(subset(ftable, Var1 == 
                                        Var2), Var1 == p)$Freq
        popCorrectRate <- popCorrectNo/pop_size
      }
      popCorrectRate_vec <- c(popCorrectRate_vec, popCorrectRate)
    }
    assign.rate.each[i, ] <- popCorrectRate_vec
  }
  assign_rate_df <- cbind(train.inds, train.loci, iters, assign.rate.all, 
                          assign.rate.each, stringsAsFactors = T)
  names(assign_rate_df)[5:ncol(assign_rate_df)] <- paste0("assign.rate.", 
                                                          pops)
  write.table(assign_rate_df, file = paste0(dir, "Rate_of_", 
                                            nrow(assign_rate_df), "_tests_", noPops, "_pops.txt"), 
              quote = F, row.names = F)
  cat("\n  Correct assignment rates were estimated!!")
  cat(paste0("\n  A total of ", nrow(assign_rate_df), " assignment tests for ", 
             noPops, " pops."))
  cat(paste0("\n  Results were also saved in a 'Rate_of_", 
             nrow(assign_rate_df), "_tests_", noPops, "_pops.txt' file in the directory."))
  return(assign_rate_df)
}
<bytecode: 0x7ff22c62d588>
  <environment: namespace:assignPOP>