# Magwene Methods
# May 19 2009
# Upated 21 Dec 2013
# UPDATED 25.4.2016

# Major error?  Reversed correlation sign by not using -cov2cor(invmatrix) ????
# see ppcor package and pcor() function.

magwene.inversion<-function(x, no.sample = NULL, 
	node.width=10, vlb = 1, vertex.cols = 'black',
	edge.mult=10, layout=layout.kamada.kawai, suppress=FALSE, 
	curve=FALSE, Out = FALSE){
	
	# function takes x = data frame of trait values
	# number in sample (needed for sig tests)
	
	# and node.width of vertex size (default is 10)
	# vlb is offset of vertex label
	# vertex.cols for altering by trait (needs work)
	# and edge.mult for multiplier on edges to make visible etc (default is 10)
	# and layout - based in plot.igraph options
	# alternative might be cirlce; see ?layout_
	# suppress = FALSE -> if TRUE, DON'T print the correlations on each edge

	# Out = TRUE returns the list of matrices etc to the screen.  
	# Otherwise, collect with assignment
	
	# check igraph
	if(require(igraph)==FALSE)
	stop("igraph not loaded")
	
	#--------------------------------------------
	# STEP 1
	# generate correlation matrix among traits
	#--------------------------------------------
	# Raw data in - scale - correlation
	if(is.data.frame(x)==TRUE){
		dd<-scale(x,scale=TRUE,center=TRUE)
		cc<-cov(dd)
		}
	else
	#--------------------------------------------
	# covariance in - cov2cor
	if(is.matrix(x)&all(diag(x)!=1)){
		cc<-cov2cor(x)
		}	
	else
	#--------------------------------------------
	# cor in.
	if(is.matrix(x)&all(diag(x)==1)){
		cc<-x
		}
	#--------------------------------------------
	
	# some stuff and checks for calculating signficance
	dimension<-dim(cc)
	if(is.null(no.sample))
	stop("You have not provided the number of observations (no.sample =) that were used to build the correlation matrix")

	# set number of observations
	subjects<-no.sample

	cat('This is your starting correlation matrix','\n')
	print(cc)

	#--------------------------------------------
	# STEP 2
	# generate inverse correlation matrix
	# values here close to 0 mean conditional independece from other traits
	invcc<-solve(cc)

	#--------------------------------------------	
	# STEP 3
	# generate a scaled inverse correlation matrix = 
	# partial correlation matrix
	# again, values here close to 0 mean conditional independence from other traits

	#--------------------------------------------	
	# UPDATED 25.4.2016
	# PCOR R package uses -cov2cor.....  reverses meaning!
	# I had it wrong - the cov2cor produces the negatives of the partial correlations
	# thus -cov2cor are the partial correlations
	# see Magwene 2001 1738
	#--------------------------------------------
	scaleinvcc <- -cov2cor(invcc)
	diag(scaleinvcc)<-1
	#print("conditional correlations")
	#print(scaleinvcc)

	#--------------------------------------------
	# Step 4 - Edge Exclusion Deviance to test for correlations = 0
	pcm<-round(scaleinvcc,3)
	pcm[upper.tri(pcm)]<-0
	pcm
	ppcm<-as.numeric(pcm)
	#print(pcm)

	# generate edge exclusion deviance table
	# AND CORRELATION STRENGTHS
	eed.sig<-numeric()
	eed.str<-numeric()
	count<-length(ppcm)

	#--------------------------------------------
	for(i in 1:count){
		# edge exclusion deviance matrix (Magwene 2001)
		eed.sig[i]<-(-subjects*log(1-ppcm[i]^2)) 
		# edge strength matrix (Magwene 2001)
		eed.str[i]<-(-0.5*log(1-ppcm[i]^2))
		}

	#--------------------------------------------
	eed.sig.mat<-matrix(eed.sig,dimension)
	#print("Edge Significnance Matrix")
	#print(eed.sig.mat)

	#--------------------------------------------
	eed.str.mat<-matrix(eed.str,dimension)
	diag(eed.str.mat)<-0
	#print("Edge Strength Matrix")
	#print(eed.str.mat)

	#--------------------------------------------
	# Significance on Chisq - THE TEST
	eed.graph.qual<-ifelse(eed.sig.mat<3.84,NA,"*")
	eed.graph<-ifelse(eed.sig.mat<3.84,0,1)
	dimnames(eed.graph)<-dimnames(cc)
	# cat("\n")
	# print("Edge Significance Graph")
	# print(eed.graph.qual)
	# cat("\n")
	# print(paste("CONNECTIVITY: ", 
		# round(sum(eed.graph)/dim(eed.graph)[1]^2, 3)))

	#--------------------------------------------
	# STRENGTHS
	strs<-eed.str.mat[eed.graph==1]
	# deals with perfect correlations -1 in pppm
	strs[strs==Inf]<-1.05+max(strs[strs!=Inf]) 

	#--------------------------------------------	
	# PLOT THE RESULTING GRAPH FROM THE MATRIX
	#--------------------------------------------
	g<-graph.adjacency(eed.graph, mode="undirected")
	g$layout<-layout

		
	# GET LABELS ONTO THE EDGES, BUT ONLY WHERE THERE ARE EDGES from tests
	# labels
	labels<-round(scaleinvcc[eed.graph==1],2)
	
	# set type to solid (1) for positives, dotted (3) for negs
	types<-ifelse(labels>0,1,2)
	E(g)$lty<-types
	E(g)$width = edge.mult*round(strs,3)+0.1

	# Suppress the edge labels if wanted, otherwise set to labels
	if(suppress==TRUE){
		edge.label=NA	# NO LABELS
		}
	else
		{E(g)$label<-labels} #set labels


	# get rid of edges that are vertex to vertex
	g<-delete.edges(g, which(labels==1))

	# make the plot
	plot(g,
		vertex.color = rgb(0,0,1, alpha=0.3),
		vertex.shape="circle",
		vertex.frame.color=NA,
		vertex.size = node.width,
		vertex.size2 = 15,
		vertex.label.dist = vlb,
		vertex.label.family = "Arial",
		vertex.label.cex = 0.8,
		edge.color='grey80',
		edge.curved = curve,
		edge.label.cex = 0.8,
		edge.label.family ="Arial")
		
	out<-list(ConditionalCorr = scaleinvcc, EdgeSig = eed.sig.mat, 
		EdgeStrength = eed.str.mat, EdgeGraph = eed.graph.qual, 
		Connectivity = round(sum(eed.graph)/dim(eed.graph)[1]^2, 3))	
		
	if(Out == TRUE){
		return(out)
	}
	else{
		print("Done - no output requested")
	}
}	