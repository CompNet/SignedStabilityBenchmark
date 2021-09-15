

###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.benchmark.input.network.folder.path <- function(n, l0, d, prop.mispl.int=NA, prop.mispl.ext=NA, prop.neg=NA)
{	
	result <- file.path(BENCHMARK.NETWORKS.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
    if(!is.na(prop.mispl.int) && !is.na(prop.mispl.ext))
        result <- file.path(result, paste0("propMisplInt=",sprintf("%.4f",prop.mispl.int), paste0("_propMisplExt=",sprintf("%.4f",prop.mispl.ext))))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	
	return(result)
}




###############################################################################
# Builds the path of a subfolder corresponding to one specific parameter set.
#
# n: number of nodes in the graph.
# l0: number of clusters in the graph.
# d: density of the graph.
# prop.mispls: proportion of misplaced links.
# prop.neg: proportion of negative links in the network.
# network.no
# algo.name
# graph.desc.name
# rep.no: repetition no
#
# returns: the folder path defined to store the network and the associated result files.
###############################################################################
get.benchmark.part.folder.path <- function(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg=NA, algo.name=NA, rep.no=NA)
{	
	result <- file.path(BENCHMARK.PARTITIONS.FOLDER,paste0(paste0("n=",n),paste0("_l0=",l0),paste0("_dens=",sprintf("%.4f",d))))
    if(!is.na(prop.mispl.int) && !is.na(prop.mispl.ext))
        result <- file.path(result, paste0("propMisplInt=",sprintf("%.4f",prop.mispl.int), paste0("_propMisplExt=",sprintf("%.4f",prop.mispl.ext))))
	if(!is.na(prop.neg))
		result <- file.path(result, paste0("propNeg=",sprintf("%.4f",prop.neg)))
	if(!is.na(algo.name))
		result <- file.path(result, algo.name)
	if(!is.na(rep.no))
	    result <- file.path(result, rep.no)
	
	return(result)
}

