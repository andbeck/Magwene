# Magwene Function COnditional Independence Graphs

This is the Magwene (2001) conditional independence calcs and plotting

Two files are provided.  

One is the magwene.inversion() function which takes a data frame, covariance matrix or correlation matrix, passes a correlation matrix onwards and generates partial correlation, edge exclusion deviance, tests and ultimately a graph/network visualisation using the package igraph.

The second is an application of the function to the data in Magwene 2001 (The Fowl dataset), reproducing his matrices and graph.