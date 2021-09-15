
library(igraph)


#################################################################
#
#################################################################
generate.param.membership = function(n, k){
	nrest = rep(0,k) # init
	nc = rep(n%/%k, k) # cluster size vector
	rest = n%%k
	if(rest>0){
		nrest[1:rest] = rep(1,rest)
		nc = nc + nrest
	}
	membership <- rep(1:k, nc)
	return(membership)
}



#################################################################
#
#################################################################
generate.complete.signed.graph <- function(membership)
{	n <- length(membership)
	m <- n*(n-1)/2
	
	# generate links and weights
	weights <- rep(NA, m)
	el <- matrix(NA, nrow=m, ncol=2)
	r <- 1
	for(i in 1:(n-1))
	{	for(j in (i+1):n)
		{	el[r,] <- c(i,j)
			if(membership[i]==membership[j])
				weights[r] <- 1
			else
				weights[r] <- -1
			r <- r + 1
		}
	}
	
	# build graph
	g <- make_empty_graph(n,directed=FALSE)
	g <- add_edges(graph=g, edges=c(t(el)), attr=list(weight=weights))
	
	return(g)
}



###############################################################################
#
###############################################################################
generate.incomplete.signed.graph <- function(membership, dens, prop.neg)
{	# init proportions
	n <- length(membership)
	qw <- 1
	# init sign probas
	p.neg <- prop.neg * dens
	p.pos <- dens - p.neg
	tlog(8,"pneg=",p.neg," p.pos=",p.pos," (total=",p.neg+p.pos,")")
	# init position probas
	p.int <- sum(sapply(1:max(membership), function(c)
		{	n <- length(which(membership==c))
			n*(n-1)/2
		})) / (n*(n-1)/2)
	p.ext <- sum(apply(t(combn(x=max(membership),m=2,simplify=TRUE)), 1, function(r)
		{	n1 <- length(which(membership==r[1]))
			n2 <- length(which(membership==r[2]))
			n1 * n2
		})) / (n*(n-1)/2)
	tlog(8,"p.int=",p.int," p.ext=",p.ext," (total=",p.int+p.ext,")")
	# check prop.neg interval
	
	lower.bound <- max(0, (dens - (1 - p.ext)) / dens)
	upper.bound <- min(1, p.ext / dens)
	tlog(8,"prop.neg (",sprintf("%.4f",prop.neg),") bounds: [",sprintf("%.4f",lower.bound)," ; ",sprintf("%.4f",upper.bound),"]")
	if(prop.neg<lower.bound | prop.neg>upper.bound)
		stop("Parameter prop.neg (",sprintf("%.4f",prop.neg),") must be in [",sprintf("%.4f",lower.bound)," ; ",sprintf("%.4f",upper.bound),"]")	
		
	# init internal probas
	p.pos.int <- p.pos / p.int
	p.none.int <- 1 - p.pos.int
	tlog(8,"Internal probas: pos=",sprintf("%.4f",p.pos.int)," none=",sprintf("%.4f",p.none.int))
	# init external probas
	p.neg.ext <- p.neg / p.ext
	p.none.ext <- 1 - p.neg.ext
	tlog(8,"External probas: neg=",sprintf("%.4f",p.neg.ext)," none=",sprintf("%.4f",p.none.ext))
	
	# draw links
	el <- t(combn(x=n, m=2, simplify=TRUE))
	comembership <- sapply(1:nrow(el),function(r) membership[el[r,1]]==membership[el[r,2]])
	norm.int <- length(which(comembership))
	norm.ext <- length(which(!comembership))
	tlog(8,"Maximal link numbers: total=",nrow(el)," internal=",norm.int," external=",norm.ext)
	tlog(8,"Expected link numbers: total=",round(nrow(el)*dens)," positive=",round(nrow(el)*dens*(1-prop.neg))," negative=",round(nrow(el)*dens*prop.neg))
	weights <- rep(NA,nrow(el))
	weights[comembership] <- sample(x=c(+1,0),size=length(which(comembership)),replace=TRUE,prob=c(p.pos.int,p.none.int))
	obs.p.pos.int <- length(which(weights[comembership]>0))/norm.int
	obs.p.none.int <- length(which(weights[comembership]==0))/norm.int
	tlog(8,"Drawn internal probas: pos=",sprintf("%.4f",obs.p.pos.int)," none=",sprintf("%.4f",obs.p.none.int))
	weights[!comembership] <- sample(x=c(-1,0),size=length(which(!comembership)),replace=TRUE,prob=c(p.neg.ext,p.none.ext))
	obs.p.neg.ext <- length(which(weights[!comembership]<0))/norm.ext
	obs.p.none.ext <- length(which(weights[!comembership]==0))/norm.ext
	tlog(8,"Drawn external probas: neg=",sprintf("%.4f",obs.p.neg.ext)," none=",sprintf("%.4f",obs.p.none.ext))
	tlog(8,"Drawn link numbers: total=",length(which(weights!=0))," positive=",length(which(weights>0))," negative=",length(which(weights<0)))
	idx <- which(weights==0)
	if(length(idx)>0)
	{	el <- el[-idx,,drop=FALSE]
		weights <- weights[-idx,drop=FALSE]
	}
	
	# build graph
	g <- make_empty_graph(n,directed=FALSE)
	g <- add_edges(graph=g, edges=c(t(el)), attr=list(weight=weights))
	
	return(g)
}


#################################################################
# It partitions a given network based on the considered algorithm name and graph type (weighted or not, etc.).
#
# n: graph size
# l0: number of cluster
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.heur.algos: the names of correlation clustering heuristic algorithms to run
# heur.reps: a sequence of values containing repetition numbers, i.e. 1, 2, 3, 4, etc.
# cor.clu.exact.algos: the names of exact correlation clustering algorithms to run
# exact.reps: a sequence of values containing repetition numbers, i.e. 1, 2, 3, 4, etc.
# keep.algo.log.files: 
# plot.layout: plot layout, such as "kamada.kawai", "fruchterman.reingold", "bn_zheng" or "circle"
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
generate.perfectly.balanced.networks = function(n, l0, d, prop.neg)
{
	membership <- generate.param.membership(n, l0)
		
    # generate the graph
	tlog(6,"Generating the graph")
	if(d<1)
		g <- tryCatch(
			generate.incomplete.signed.graph(membership, d, prop.neg),
			error=function(e) NA
		)
	else
		g <- tryCatch(
				generate.complete.signed.graph(membership),
				error=function(e) NA
		)
		
		
	if(!all(is.na(g)))
	{	# possibly create the output folder
	    prop.mispl.int = 0
	    prop.mispl.ext = 0
		folder <- get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg)
		dir.create(folder,showWarnings=FALSE,recursive=TRUE)

		# record the graph
		file.graph <- file.path(folder,paste0(GRAPH.FILENAME,".graphml"))
		tlog(6,"Recording the graph in file ",file.graph)
		write_graph(graph=g,file=file.graph,format="graphml")
		
	    # export using a format compatible with pILS
	    t <- get.edgelist(graph=g) - 1	# start numbering nodes at zero
	    t <- cbind(t,E(g)$weight)		# add weights as the third column
	    file.graph <- file.path(folder,paste0(GRAPH.FILENAME,".G"))
	    write.table(data.frame(vcount(g),nrow(t)), file=file.graph, append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) # write header
	    write.table(t, file=file.graph, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE) # write proper graph
				
	}

}



#################################################################
# It is the starting method in the aim of partitioning the considered networks. 
#   It handles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.heur.algos: the names of correlation clustering heuristic algorithms to run
# heur.reps: a sequence of values containing repetition numbers, i.e. 1, 2, 3, 4, etc.
# cor.clu.exact.algos: the names of exact correlation clustering algorithms to run
# exact.reps: a sequence of values containing repetition numbers, i.e. 1, 2, 3, 4, etc.
# keep.algo.log.files: 
# plot.format: plot format(s), such as PDF, PNG or NA (it means 'no plotting')
# plot.layout: plot layout, such as "kamada.kawai", "fruchterman.reingold", "bn_zheng" or "circle"
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
generate.all.perfectly.balanced.networks = function(graph.sizes, d, l0.values, prop.negs)
{
    tlog("starts generating all perfectly balanced networks")
    for (n in graph.sizes) {
        tlog(4, "generating all perfectly balanced networks => n: ", n)
        
        for (l0 in l0.values) {
        tlog(4, "generating all perfectly balanced networks => l0: ", l0)
        
        
            my.prop.negs = prop.negs # if we do not do that, for each n value, prop negs will not be the initial value(s)
            if(is.na(my.prop.negs) && d == 1){
                my.prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
            }
            
            for (prop.neg in my.prop.negs) {
                tlog(12, "generating all perfectly balanced networks => prop.neg: ", prop.neg)
                
                generate.perfectly.balanced.networks(n, l0, d, prop.neg)
            }
            
        }

    }
}
