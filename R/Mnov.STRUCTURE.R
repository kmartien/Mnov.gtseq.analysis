library("swfscMisc")
library("strataG")

# Load data
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")
project <- 'final_sams.no.dupes'
minReads <- 20

strat.scheme <- 'herd'

# stratify and select strata to include
g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

# choose which strata I want to include
g.HI.CentAm <- g.stratified[,,which(getStrataNames(g.stratified) %in% c('CentAm-CA.OR.WA', 'HI-SEAK.NBC'))]

# Run STRUCTURE
sr.k1 <- structureRun(g.HI.CentAm, k.range = 1, num.k.rep = 1, label = strat.scheme, delete.files = FALSE, num.cores = 3, 
                    burnin = 1000, numreps = 10000, noadmix = FALSE, freqscorr = TRUE, 
                    pop.prior = NULL)

save(sr, file = paste0("results-R/", project, '_minReads.', minReads, '.', strat.scheme, "_sr.rda"))

# Calculate Evanno metrics
evno <- evanno(sr)
print(evno)

#haven't updated past here
strata.num <- na.omit(as.numeric(strata.df[2]))
strata.num <- strata.num[1:166] + 5

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp(sr, k = 2)
clumpp2$orig.pop <- strata.num
mean.assignment <- aggregate(clumpp2[,4:5],list(clumpp2$orig.pop),mean)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp2, sort.probs=F,label.pops=F,col=c("#009E73","#D55E00"))

# Run CLUMPP to combine runs for K = 3
clumpp3 <- clumpp.run(sr, k = 3)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp3)

# Run CLUMPP to combine runs for K = 4
clumpp4 <- clumpp(sr, k = 4)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp4)

# Run CLUMPP to combine runs for K = 5
clumpp5 <- clumpp.run(sr, k = 5)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp5)

# Run CLUMPP to combine runs for K = 6
clumpp6 <- clumpp.run(sr, k = 6)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp6)

save.image("AS177_all.rdata")

