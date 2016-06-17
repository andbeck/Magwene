# magwene test
library(ppcor)
library(igraph)

source("Magwene2001.R")

FowlCors<-c(1,NA,NA,NA,NA,NA,
			0.583,1,NA,NA,NA,NA,
			0.621,0.584,1,NA,NA,NA,
			0.603,0.526,0.937,1,NA,NA,
			0.569,0.515,0.877,0.878,1,NA,
			0.602,0.548,0.874,0.894,0.926,1)
FowlPrep<-matrix(FowlCors, nrow = 6, ncol = 6, byrow=TRUE)
FowlPrep
FowlPrep[upper.tri(FowlPrep)]<-t(FowlPrep)[upper.tri(t(FowlPrep))]
Fowl<-FowlPrep

# Step 1 - inverse correlation
invcm<-solve(Fowl)

pcor<--cov2cor(invcm)

invcm
pcor

magwene.inversion(Fowl)