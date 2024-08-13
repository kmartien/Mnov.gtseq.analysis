library("swfscMisc")
library("strataG")

# Load data
load("data/data.for.PSRG.2024.rda")

strat.scheme <- "wint.area"

# stratify and select strata to include
g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

#only do pairwise tests on strata with 10 or more samples
g.10 <- g.stratified[,,which(getNumInd(g.stratified, by.strata = TRUE)$num.ind >= 10)]

# Run STRUCTURE
sr <- structureRun(g.10, k.range = 2:7, num.k.rep = 5, label = strat.scheme, delete.files = FALSE, num.cores = 3, 
                    burnin = 1000, numreps = 10000, noadmix = FALSE, freqscorr = TRUE, 
                    pop.prior = NULL)

save(sr, file = paste0("results-R/", strat.scheme, "_sr.rda"))

# Calculate Evanno metrics
evno <- evanno(sr)
print(evno)

#haven't updated past here
strata.num <- na.omit(as.numeric(strata.df[2]))
strata.num <- strata.num[1:166] + 5

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp.run(sr, k = 2)
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
clumpp4 <- clumpp.run(sr, k = 4)
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

