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
pcorx<--cov2cor(invcm)

# checks against Magwene 2001
invcm
pcorx

# magwene fowl example with graphing
source("Magwene2001.R")
source("Magwene2001_Mistake.R")

# ppcor method
ppcorMethod<-function(covmat){
	icvx <- solve(covmat)
    pcor <- -cov2cor(icvx)
    diag(pcor) <- 1
    return(pcor)
}

# 276 chickens
# alternative layouts see ?layout_
# play with edge multiplier to see different edges
# solid = positive, dashed = negative
par(mfrow=c(1,2))
mag_wrong<-magwene.inversion.wrong(Fowl, no.sample = 276,
	edge.mult = 30, layout = layout_with_kk,
	Out =TRUE)

mag_right<-magwene.inversion(Fowl, no.sample = 276,
	edge.mult = 30, layout = layout_with_kk,
	Out =TRUE)

#COMPARE OLD wrong, new and ppcor methods
mag_wrong$PartialCorr
mag_right$PartialCorr
ppcorMethod(Fowl)


# PCOR package example
y.data <- data.frame(
				hl=c(7,15,19,15,21,22,57,15,20,18),
				disp=c(0.000,0.964,0.000,0.000,0.921,0.000,0.000,1.006,0.000,1.011),
				deg=c(9,2,3,4,1,3,1,3,6,1),
				BC=c(1.78e-02,1.05e-06,1.37e-05,7.18e-03,0.00e+00,0.00e+00,0.00e+00
              ,4.48e-03,2.10e-06,0.00e+00)
			)
y.data.scale<-scale(y.data)

# partial correlation
pcor(y.data)$estimate
magwene.inversion(y.data, no.sample=10, Out = TRUE)$PartialCorr

# stat/edge test
pcor(y.data)$statistic
magwene.inversion(y.data, no.sample=10, Out = TRUE)$EdgeSig

# p-value
pcor(y.data)$p.value
magwene.inversion(y.data, no.sample=10, Out = TRUE)$EdgeGraph

# estimate of partial is independent of raw or scaled cov/cor 
pcor(y.data)$estimate
pcor(y.data.scale)$estimate