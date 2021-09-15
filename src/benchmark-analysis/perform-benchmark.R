


#################################################################
# Applies one of the following methods onto perturbed signed graphs:
#   - CoNSCC(r_{max}=3)
#   - CoNSCC(r_{max}=3) without MVMO pruning
#   - CoNSCC(r_{max}=4)
#   - CoNSCC(r_{max}=4) without MVMO pruning
# The two input parameters 'maxNbEdit' and 'pruning.without.MVMO' determine which one to use.
#
# 
# n: number of vertices
# l0: number of modules
# d: graph density
# prop.mispl.int: proportion of misplaced internal links
# prop.mispl.ext: proportion of misplaced extenal links
# prop.neg: proportion of negative edges
# maxNbEdit: a distance value, which is the input parameter of the CoNS method.
# pruning.without.MVMO: whether or not the MVMO pruning strategies are used or not.
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
perform.benchmark = function(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, 
		maxNbEdit, pruning.without.MVMO, force)
{
    graph.desc.name = "signed-unweighted"
    in.folder = get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg)
	print(in.folder)
	
	#algo.name = get.RNSCC.code(maxNbEdit=maxNbEdit, pruning.without.MVMO)
	algo.name = get.CoNSCC.code(maxNbEdit=maxNbEdit, pruning.without.MVMO)
    part.folder.EnumCC = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, rep.no=NA)
    if(!dir.exists(part.folder.EnumCC))
		dir.create(part.folder.EnumCC, showWarnings = FALSE, recursive = TRUE)
	
	membership.filepath = file.path(part.folder.EnumCC, paste0(MBRSHP.FILE.PREFIX,0,".txt"))
	exec.time.filepath = file.path(part.folder.EnumCC,EXEC.TIME.FILENAME)
	print(membership.filepath)
	print(exec.time.filepath)
	if(!file.exists(membership.filepath) || !file.exists(exec.time.filepath) || force)
	{ # we are interested in the execution time
		file.create(file.path(part.folder.EnumCC,"allResults.txt")) # an empty file
	
		part.folder.ExCC = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, COR.CLU.ExCC)
		membership.filepath.from = file.path(part.folder.ExCC, paste0(MBRSHP.FILE.PREFIX,0,".txt"))
		# print(membership.filepath.from)
		
		if(file.exists(membership.filepath.from)) {
			file.copy(from=membership.filepath.from, to=membership.filepath, overwrite = TRUE)
			apply.partitioning.algorithm(part.folder.EnumCC, in.folder, algo.name, graph.name, plot.formats)
		}
	
	}
	
}



#################################################################
# Applies the following methods onto perturbed signed graphs:
#   - CoNSCC(r_{max}=3)
#   - CoNSCC(r_{max}=3) without MVMO pruning
#   - CoNSCC(r_{max}=4)
#   - CoNSCC(r_{max}=4) without MVMO pruning
#
# graph.sizes: number of vertices
# d: graph density
# l0.values: number of modules
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
perform.all.benchmark = function(graph.sizes, d, l0.values, force)
{
    	prop.mispl.int = 0
        prop.mispl.ext = 0

    	tlog("starts fast jump model stats")
    	for(n in graph.sizes){
    	    for(l0 in l0.values){
    	        tlog(8, "performing benchmark => n: ", n, ", l0: ",l0)
			
                in.folder = get.benchmark.input.network.folder.path(n, l0, d, NA, NA, NA)
                dirs = list.dirs(path = in.folder, full.names = FALSE, recursive = FALSE)
                print(dirs)
    				
                for(dir in dirs){
                    # dir = 'propMisplInt=0.0000_propMisplExt=0.0000'
                    parts = unlist(strsplit(dir, "_"))
                    prop.mispl.int = as.numeric(unlist(strsplit(parts[1],"="))[2])
                    prop.mispl.ext = as.numeric(unlist(strsplit(parts[2],"="))[2])
                    in.folder = get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, NA)
                    dirs2 = list.dirs(path = in.folder, full.names = FALSE, recursive = FALSE)
                    
                    for(dir2 in dirs2){
						prop.neg = as.numeric(unlist(strsplit(dir2,"="))[2])
                        #print(force)
					
						maxNbEdit = 3
						tlog("maxNbEdit: ",maxNbEdit)
						pruning.without.MVMO = FALSE
			            perform.benchmark(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, maxNbEdit, pruning.without.MVMO, force)
						pruning.without.MVMO = TRUE
						perform.benchmark(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, maxNbEdit, pruning.without.MVMO, force)
					
						maxNbEdit = 4
						tlog("maxNbEdit: ",maxNbEdit)
						pruning.without.MVMO = FALSE
						perform.benchmark(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, maxNbEdit, pruning.without.MVMO, force)
						pruning.without.MVMO = TRUE
						perform.benchmark(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, maxNbEdit, pruning.without.MVMO, force)
						
					
					}
    	        }
    	    }
    	}
}
