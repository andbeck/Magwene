# Magwene Methods
# May 19 2009
# Upated 21 Dec 2013
# UPDATED 25.4.2016

# Major error?  Reversed correlation sign by not using -cov2cor(invmatrix) ????
# see ppcor package and pcor() function.

magwene.inversion<-function(x,node.width=10, edge.mult=6, vertex.cols = 'black', layout=layout.circle, 
	alpha=FALSE, suppress=FALSE, curve=FALSE, Out = FALSE){
	
	# function takes x = data frame of trait values
	# and node.width of graph node (default is 10)
	# and edge.mult for multiplier on edges to make visible etc (default is 6)
	# and layout - based in plot.igraph options - default is layout.circle
	# alternative might be layout.kamada.kawai
	# alpha selects between using thickness of edges or alpha color for edge emphasis
	# alpha = FALSE is default and changes edge thickness via edge.multiplier
	# suppress = FALSE - if true, don't print the correlations on each edge
	# Out = TRUE returns the list of matrices etc to the screen.  Otherwise, collect with assignment
	
	require(igraph)
	
	# STEP 1
	# generate correlation matrix among traits
	if(is.data.frame(x)==TRUE){
		dd<-scale(x,scale=TRUE,center=TRUE)
		cc<-cov(dd)
		subjects<-dim(dd)[1]
		}
	else
	if(is.matrix(x)&all(diag(x)!=1)){
		cc<-cov2cor(x)
		}	
	else
	if(is.matrix(x)&all(diag(x)==1)){
		cc<-x
		}
	
	dimension<-dim(cc)
	
	cat('This is your starting correlation matrix','\n')
	print(cc)

	
	# STEP 2
	# generate inverse correlation matrix
	# values here close to 0 mean conditional independece from other traits
	invcc<-solve(cc)
	
	# STEP 3
	# generate a scaled inverse correlation matrix = 
	# partial correlation matrix
	# again, values here close to 0 mean conditional independence from other traits
	
	# UPDATED 25.4.2016
	# PCOR R package uses -cov2cor.....  reverses meaning!
	# I had it wrong - the cov2cor produces the negatives of the partial correlations
	# thus -cov2cor are the partial correlations
	# see Magwene 2001 1738
	
	scaleinvcc <- -cov2cor(invcc)
	#print("conditional correlations")
	#print(scaleinvcc)

	# Step 4 - Edge Exclusion Deviance to test for correlations = 0
	pcm<-round(scaleinvcc,3)
	pcm[upper.tri(pcm)]<-0
	pcm
	ppcm<-as.numeric(pcm)
	#print(pcm)

	# generate edge exclusion deviance table
	# AND CORRELATION STRENGTHS
	subjects<-dim(dd)[1]
	eed.sig<-numeric()
	eed.str<-numeric()
	count<-length(ppcm)

	for(i in 1:count){
		eed.sig[i]<-(-subjects*log(1-ppcm[i]^2))
		eed.str[i]<-(-0.5*log(1-ppcm[i]^2))
		}

	eed.sig.mat<-matrix(eed.sig,dimension)
	diag(eed.sig.mat)<-0
	#print("Edge Significnance Matrix")
	#print(eed.sig.mat)

	eed.str.mat<-matrix(eed.str,dimension)
	diag(eed.str.mat)<-0
	#print("Edge Strength Matrix")
	#print(eed.str.mat)

	# Significance on Chisq - THE TEST
	eed.graph<-ifelse(eed.sig.mat<3.84,0,1)
	dimnames(eed.graph)<-dimnames(cc)
	#print("Edge Significance Graph")
	#print(eed.graph)
	cat("\n")
	print(paste("CONNECTIVITY: ", round(sum(eed.graph)/dim(eed.graph)[1]^2, 3)))

	# STRENGTHS
	strs<-eed.str.mat[eed.graph==1]
	strs[strs==Inf]<-1.05+max(strs[strs!=Inf]) # deaks with perfect correlations -1 in pppm
	
	# TRYING TO GET LABELS ONTO THE EDGES, BUT ONLY WHERE THERE ARE EDGES
	potential<-scaleinvcc[lower.tri(scaleinvcc)]
	codes<-eed.graph[lower.tri(eed.graph)==1]
	el<-round(potential[codes==1],3)

	# PLOT THE RESULTING GRAPH FROM THE MATRIX
	g<-graph.adjacency(eed.graph, mode="undirected")
	g$layout<-layout
	E(g)$lty=ifelse(el>0,1,3)
	# COULD change the edge size to be constant but be in alpha by the strs
	# http://rdatamining.wordpress.com/2012/05/17/an-example-of-social-network-analysis-with-r-using-package-igraph/
	
	# alpha colour verses edge thickness setting
	if(alpha==TRUE){
		edge.color=rgb(0,0,1,alpha=round(strs,2)/10)
		edge.width=5}
	else
		{edge.color="grey80"
		edge.width=edge.mult*round(strs,3)}

	# Suppress the edge labels if wanted
	if(suppress==TRUE){
		edge.label=NA	# NO LABELS
		}
	else
		{edge.label=el}
	

	plot(g,vertex.label=rownames(eed.graph),
		vertex.color=vertex.cols, # colors by trait
		#vertex.color = rgb(0,0,0.8, alpha=0.6)
		vertex.shape="circle",
		vertex.frame.color=NA,
		edge.color=edge.color,
		edge.width=edge.width,
		vertex.size=node.width,
		vertex.size2=15,
		vertex.label.dist=0.5,
		vertex.label.family="Arial",
		vertex.label.cex=0.8,
		edge.curved=curve,
		edge.label.cex=0.8,
		edge.label.family="Arial",
		edge.label=edge.label)
		
	out<-list(ConditionalCorr = scaleinvcc, EdgeSig = eed.sig.mat, 
		EdgeStrength = eed.str.mat, EdgeGraph = eed.graph, 
		Connectivity = round(sum(eed.graph)/dim(eed.graph)[1]^2, 3))	
	if(Out == TRUE){
		return(out)
	}
	else{
		print("Done - no output requested")
	}
}	