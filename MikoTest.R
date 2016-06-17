#MikoINt
# Data from Dryad.

library(dplyr)
source('~/Documents/Repos/Magwene/Magwene2001.R')

dats<-read.csv('~/Desktop/MikoInt.csv')

albifrons<-filter(dats, species == 'albifrons')
albFish<-filter(albifrons, predator=="vertebrate")
albInv<-filter(albifrons, predator=="invertebrate")

use.albFish<-na.omit(albFish[,c(7:26)])
use.albInv<-na.omit(albInv[,c(7:26)])
head(use.albFish)
dim(use.albFish)
dim(use.albInv)

par(mfrow=c(1,2))
magwene.inversion(use.albFish, no.sample = 104, Out = TRUE,suppress=TRUE)
magwene.inversion(use.albInv, no.sample = 76, Out = TRUE,suppress=TRUE)
