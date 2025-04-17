library(strataG)
library(genepop)
library(tidyverse)
library(ggplot2)
load("data/gtypes_final.sams.no.dupes_minReads.20.rda")
source('R/functions/fstatWrite.R')

strat.scheme <- 'herd'

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)
g.stratified <- g.stratified[,,which(getNumInd(g.stratified, by.strata = TRUE)$num.ind >= 9)]
g.stratified <- g.stratified[,,-which(getStrataNames(g.stratified) == 'SMx.CA.OR.WA.SBC')]
fstatWrite(g.stratified, label = strat.scheme)

g.westcoast <- g.stratified[,,which(getStrataNames(g.stratified) %in% c('CentAm.CA.OR.WA.SBC', 'MnMx.CA.OR.WA.SBC', 'HI.CA.OR.WA.SBC'))]
fstatWrite(g.westcoast, label = 'westcoast.herds')

strat.scheme <- 'wint.area'

g.stratified <- stratify(g, scheme = strat.scheme, drop = TRUE)
g.stratified <- g.stratified[,,which(getNumInd(g.stratified, by.strata = TRUE)$num.ind >= 9)]
g.stratified <- g.stratified[,,-which(getStrataNames(g.stratified) == 'SMx')]
fstatWrite(g.stratified, label = strat.scheme)
