library(adegenet)
library(strataG)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(viridis)
library(RColorBrewer)
source('R/functions/DAPC.fit.and.predict.R')

# Load data
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")
project <- 'final_sams.no.dupes'
minReads <- 20

#############################################
# stratified by herds
strat.scheme <- 'herd'
g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

# choose which strata I want to include
#g.stratified <- g.stratified[,,which(getStrataNames(g.stratified) %in% c('CentAm.CA.OR.WA', 'HI.SEAK.NBC', 'MnMx.CA.OR.WA'))]
g.stratified <- g.stratified[,,
                             filter(getNumInd(g.stratified, by.strata = TRUE), num.ind >= 9) |> 
                               pull(stratum)]

genind.strat <- gtypes2genind(g.stratified)

#exclude SMex because it is a mixed stratum
genind.no.SMex <- genind.strat[-which(pop(genind.strat) == 'SMx.CA.OR.WA.SBC')]
col <- seasun(length(levels(pop(genind.no.SMex))))[c(4,2,3,5,6,1,7,8)]

# cross-validation to choose number of PCs and execute DAPC with the chosen number
mat <- tab(genind.no.SMex, NA.method = 'mean')
grp <- pop(genind.no.SMex)
xval.herd <- xvalDapc(mat, grp, n.pca.max = 200, 
                 training.set = 0.9, result = 'groupMean', center = TRUE,
                 scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = TRUE)
res <- xval.herd$DAPC

dapc.loov <- predictAllIDsDAPC(genind.no.SMex, n.da = 8, n.pca = 100)
save(dapc.loov, file = 'results-R/dapc_loov_res.rda')

# print results
jpeg(filename = paste0('results-raw/DAPC_scatter_', strat.scheme, '.jpg'))
scatter.dapc(res, col = col, scree.da = FALSE)
dev.off()
jpeg(filename = paste0('results-raw/DAPC_assign_', strat.scheme, '.jpg'))
table.value(table(res$assign, grp), col.lab = levels(grp))
dev.off()
mean(as.character(res$assign) == grp)

# assign SMEX samples to herd

genind.SMex <- genind.strat[which(pop(genind.strat) == 'SMx.CA.OR.WA.SBC')]

pred.SMex <- predict.dapc(res, newdata = genind.SMex)
#col <- seasun(length(levels(pop(genind.no.SMex))))[c(4,2,3,5,6,1,7,8)]
#col <- brewer.pal(8, 'Set1')
col.points <- transp(col[as.integer(pop(genind.no.SMex))],.7)
jpeg(filename = 'results-raw/SMex_assignment_herd_DAPC.jpg')
scatter(res, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, #xlim=c(-10,10), 
        legend=TRUE, posi.leg = 'bottomleft')
par(xpd=TRUE)
points(res$ind.coord[,1], res$ind.coord[,2], pch=20,
       col=col.points, cex=1)
col.sup <- col[as.integer(pred.SMex$assign)]
points(pred.SMex$ind.scores[,1], pred.SMex$ind.scores[,2], pch=15,
       col=transp(col.sup,.7), cex=2)
#add.scatter.eig(res$eig,15,1,2, posi="bottomright", inset=.02)
dev.off()

#############################################
# stratified by winter areas
strat.scheme <- 'wint.area'
g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

# choose which strata I want to include
g.stratified <- g.stratified[,,
                             filter(getNumInd(g.stratified, by.strata = TRUE), num.ind >= 10) |> 
                               pull(stratum)]

genind.strat <- gtypes2genind(g.stratified)

#exclude SMex because it is a mixed stratum
genind.no.SMex <- genind.strat[-which(pop(genind.strat) == 'SMx')]

mat <- tab(genind.no.SMex, NA.method = 'mean')
grp <- pop(genind.no.SMex)
xval.wint <- xvalDapc(mat, grp, n.pca.max = 200, 
                 training.set = 0.9, result = 'groupMean', center = TRUE,
                 scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = TRUE)
res <- xval.wint$DAPC

# print results
jpeg(filename = paste0('results-raw/DAPC_scatter_', strat.scheme, '.jpg'))
wint.cols <- col[c(4, 2, 3, 8)]
scatter.dapc(res, scree.da = FALSE, col = wint.cols, legend = TRUE)
dev.off()
jpeg(filename = paste0('results-raw/DAPC_assign_', strat.scheme, '.jpg'))
table.value(table(res$assign, grp), col.lab = levels(grp))
dev.off()
mean(as.character(res$assign) == grp)

# assign SMEX samples to herd

genind.SMex <- genind.strat[which(pop(genind.strat) == 'SMx')]

pred.SMex <- predict.dapc(res, newdata = genind.SMex)
col.points <- transp(wint.cols[as.integer(pop(genind.no.SMex))],.7)
jpeg(filename = 'results-raw/SMex_assignment_wint_DAPC.jpg')
scatter(res, col=wint.cols, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, #xlim=c(-10,10), 
        legend=TRUE, posi.leg = 'bottomleft')
par(xpd=TRUE)
points(res$ind.coord[,1], res$ind.coord[,2], pch=20,
       col=col.points, cex=1)
col.sup <- wint.cols[as.integer(pred.SMex$assign)]
points(pred.SMex$ind.scores[,1], pred.SMex$ind.scores[,2], pch=15,
       col=transp(col.sup,.7), cex=2)
#add.scatter.eig(res$eig,15,1,2, posi="bottomright", inset=.02)
dev.off()

pws.dist <- dist(res$grp.coord)


#############################################
# stratified by feeding areas
strat.scheme <- 'feed.area'
g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)

# choose which strata I want to include
g.stratified <- g.stratified[,,
                             filter(getNumInd(g.stratified, by.strata = TRUE), num.ind >= 10) |> 
                               pull(stratum)]

genind.strat <- gtypes2genind(g.stratified)

mat <- tab(genind.strat, NA.method = 'mean')
grp <- pop(genind.strat)
xval.feed <- xvalDapc(mat, grp, n.pca.max = 200, 
                      training.set = 0.9, result = 'groupMean', center = TRUE,
                      scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = TRUE)
res <- xval.feed$DAPC

# print results
jpeg(filename = paste0('results-raw/DAPC_scatter_', strat.scheme, '.jpg'))
scatter.dapc(res)
dev.off()
jpeg(filename = paste0('results-raw/DAPC_assign_', strat.scheme, '.jpg'))
table.value(table(res$assign, grp), col.lab = levels(grp))
dev.off()
mean(as.character(res$assign) == grp)
