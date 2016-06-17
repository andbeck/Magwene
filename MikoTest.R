#MikoINt
# Data from Dryad.

library(dplyr)
source('~/Documents/Repos/Magwene/Magwene2001.R')
source('~/Documents/Repos/Magwene/Magwene2001_Mistake.R')

dats<-read.csv('~/Desktop/MikoInt.csv')

albifrons<-filter(dats, species == 'albifrons')
albFish<-filter(albifrons, predator=="vertebrate")
albInv<-filter(albifrons, predator=="invertebrate")

use.albFish<-na.omit(albFish[,c(7:26)])
use.albInv<-na.omit(albInv[,c(7:26)])
head(use.albFish)
dim(use.albFish)
dim(use.albInv)

par(mfrow=c(2,2))
alb_1<-magwene.inversion(use.albFish, no.sample = 104, Out = TRUE,suppress=TRUE)
alb_2<-magwene.inversion(use.albInv, no.sample = 76, Out = TRUE,suppress=TRUE)

alb_3<-magwene.inversion.wrong(use.albFish, no.sample = 104, Out = TRUE,suppress=TRUE)
alb_4<-magwene.inversion.wrong(use.albInv, no.sample = 76, Out = TRUE,suppress=TRUE)

alb_1$ConditionalCorr
alb_2$ConditionalCorr
pcor(use.albFish)$estimate