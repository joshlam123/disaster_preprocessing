rm(list=ls())
setwd("/Users/soojun/Dropbox/Capstone/rainfall/rainfall model/") #set where to get dataset from

library(survival)
#install.packages("SPREDA")
#install.packages("extRemes")
library(SPREDA)
library(extRemes)


#The R functions for performing the extreme value analysis (for maxima)
#can be downloaded from: 
source("http://www.datall-analyse.nl/R/eva_max.R")


##data
#yearly maximum rainfall (in mm/day)
data <- read.csv("STASIUN METEOROLOGI MARITIM TANJUNG PRI.csv")
rainfall <- data$Rainfall

#explore data visually
hist(rainfall)

#Gumbel
gumbelmod <- fevd(x=rainfall, type="Gumbel", method="MLE")
summary(gumbelmod)

#extract MLEs (these are needed for the remaining part of the analysis)
muG <- gumbelmod$results$par[1]
sigmaG <- gumbelmod$results$par[2]


##probability plot for largest values (=maxima)
probplot(values=rainfall, model=gumbelmod, varname="Rainfall (mm)",
         alpha=1-.95, dist="gumbel")


##cumulative probability plot for largest values (including return period)
#note: for the rainfall data, a return period of 20 means that once every
#20 years the wind speed is (on average) expected to be larger than
#muG+qlev(1-1/20)*sigmaG=49.4 mph (in case of a Gumbel distribution)
#(this expected rainfall is also called the return level)

twentyyears <- muG+qlev(1-1/20)*sigmaG

#Gumbel distribution
cprobplot(values=rainfall, model=gumbelmod, varname="Rainfall (mm/day)",
          alpha=1-.95, dist="gumbel")

##diagnostics
#q-q and p-p plot for Gumbel
QQplot(values=rainfall, mu=muG, sigma=sigmaG, dist="gumbel")
PPplot(values=rainfall, mu=muG, sigma=sigmaG, dist="gumbel")
