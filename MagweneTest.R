# magwene test on Magwene Fowl Matrix.
library(ppcor)
library(igraph)

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

# checks against Magwene 2001
invcm
pcor

# magwene fowl example with graphing
source("Magwene2001.R")

# 276 chickens
# alternative layouts see ?layout_
# play with edge multiplier to see different edges
# solid = positive, dashed = negative
magwene.inversion(Fowl, no.sample = 276,
	edge.mult = 30, layout = layout.kamada.kawai,
	Out =TRUE)