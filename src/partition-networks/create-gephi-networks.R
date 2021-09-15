


#################################################################
# It first changes the structure of the "edges" nodes for compatibility reasons in Gephi visalisation:
#   instead of showing the edge weights with with negative and positive sign as a real number,
#   it separates this information into two parts: the absolute value fo the real number, and its sign (positive or negative).
#   An example:
#  <edge source="n2" target="n3">
#    <data key="e_weight">1</data>
#    <data key="e_sign">1</data>
#  </edge>
#
# Then, it adds partition information into nodes: If there are 5 optimal solutions associated with a signed network,
#   then, all these 5 partitions will be stored. EXample:
#     <node id="n1">
#       <data key="v_ExCC-sol0">1</data>
#       <data key="v_ExCC-sol1">3</data>
#     </node>
#
#
# n: node number, i.e graph size
# l0: number of initial clusters
# d: density
# prop.mispl: proportion of misplaced links
# prop.neg: proportion of negative links
# network.no: network id (the identifiers start from 1)
# cor.clu.exact.algos: correlation clustering algorithms to be considered 
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
create.gephi.network = function(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, force)
{
	source("src/define-graphml-operations.R")
	
	net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
	
	tlog(16, "creating gephi networks => algo.name: ", cor.clu.exact.algo)
	
	for(graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)){
		tlog(20, "creating gephi networks => graph.desc.name: ", graph.desc.name)
		
		part.folder = get.part.folder.path(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, graph.desc.name)
		if(!dir.exists(part.folder))
			dir.create(path=part.folder, showWarnings=FALSE, recursive=TRUE)
		
		tlog(20, "creating gephi networks in ", part.folder)		
		graph.name = paste0(graph.desc.name, ".graphml")
		gephi.graph.name = paste0(graph.desc.name, "-gephi.graphml")
		
		f.path = file.path(net.folder,graph.name)
		new.f.path = file.path(part.folder,gephi.graph.name)
		newGraphFile = addSignAttrIntoGraphmlFile(f.path)
		saveXML(newGraphFile, file=new.f.path)
		print("DONE: added signed attributes")
		
		# -------------
		
		# add partition info into gephi file
		partitions = load.membership.files(part.folder)
		m = length(partitions) # nb partition
		for(i in 1:m){
			attr.name = paste(cor.clu.exact.algo,"-sol",i)
			newGraphFile = addPartitionInfoIntoNodes(new.f.path, attr.name, partitions[[i]])
			saveXML(newGraphFile, file=new.f.path)
		}
		print("DONE: added partitioning information")
	
	}
	
}





#################################################################
#
# It is the starting method in the aim of creating gephi network files containing all partition information. 
#   It hanles all networks by graph.sizes,  prop.mispls, my.prop.negs and in.rand.net.folders
#
# graph.sizes: a vector of values regarding graph sizes to be considered
# d: density (it is a single value)
# l0: number of clusters to be considered (it is a single value)
# prop.mispls: a vector of values regarding proportion of misplaced links
# prop.negs: a vector of values regarding proportion of negative links (for now, it is not operational)
# in.rand.net.folders: a vector of values regarding input random graph folders. Sequantial integers (1, .., 10)
# cor.clu.exact.algo: the name of correlation clustering algorithm to run
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
#
##################################################################
create.gephi.networks = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders, cor.clu.exact.algo, force)
{
	tlog("starts creating gephi networks")
	for(n in graph.sizes){
		tlog(4, "creating gephi networks => n: ", n)
		
		for(prop.mispl in prop.mispls){
			tlog(8, "creating gephi networks => prop.mispl: ", prop.mispl)
			
            my.prop.negs = prop.negs # if we do not do that, for each n value, prop negs will not be the initial value(s)
            if(is.na(my.prop.negs) && d == 1){
                my.prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
            }
            
            for (prop.neg in my.prop.negs) {
				tlog(12, "creating gephi networks => prop.neg: ", prop.neg)
				
				net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no=NA)
				if(dir.exists(net.folder)){
				
					for(network.no in in.rand.net.folders){
						tlog(16, "creating gephi networks => network.no: ", network.no)
						
						create.gephi.network(n, l0, d, prop.mispl, prop.neg, network.no, cor.clu.exact.algo, force)
					}
				}
				
			}
			
		}
		
	}
	
}
