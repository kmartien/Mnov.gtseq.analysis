library(dplyr)
library(Mnov.GTseq.data)
library(swfscMisc)
data("GTseq.samps.final")

Mnov.qa.qc.dat <- select(GTseq.samps.final, c(AnimalID, LABID, SEX, HAP, mito.haps))
Mnov.qa.qc.dat$AnimalID <- paste0("z0", zero.pad(as.numeric(Mnov.qa.qc.dat$AnimalID)))
Mnov.qa.qc.dat$LABID <- paste0("z0", zero.pad(as.numeric(Mnov.qa.qc.dat$LABID)))

Runs <- list("RunMS58", "GTseq.prod.all", "GTseq.val.all", "RunMS51", "RunMS54", "RunMS43", "RunMS45.90s")

x <- data.frame(do.call(cbind, lapply(Runs, function(r){
  files <- list.files(path = paste0("data-raw/bam.files/", r), pattern = "*.bam$") %>% substr(start = 1, stop = 8)
  in.run <- as.numeric(Mnov.qa.qc.dat$LABID %in% files) %>% case_match(0 ~ "", 1 ~ "1")
  #in.run <- mutate(in.run, recode("TRUE" = "1", "FALSE" = ""))
  in.run[which(Mnov.qa.qc.dat$LABID %in% files[duplicated(files)])] <- "2"
  return(in.run)
})))
names(x) <- Runs
x$LABID <- Mnov.qa.qc.dat$LABID

Mnov.qa.qc.dat <- full_join(Mnov.qa.qc.dat, x, by = "LABID")

write.csv(Mnov.qa.qc.dat, file = "data-raw/QA.QC/Mnov GTseq QA QA.csv")
