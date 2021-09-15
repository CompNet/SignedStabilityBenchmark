


##################################################################
## Reads the execution time result of a given partitioning method.
##
## n: number of vertices
## l0: number of modules
## d: graph density
## prop.mispl.int: proportion of misplaced internal links
## prop.mispl.ext: proportion of misplaced extenal links
## prop.neg: proportion of negative edges
## algo.name: a partitioning algo name
## force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
##
###################################################################
#read.exec.time.result = function(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, force)
#{
#	exec.time = NA
#    part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, NA)
#    exec.filepath = file.path(part.folder, EXEC.TIME.FILENAME)
#    if(file.exists(exec.filepath)){
#       exec.time = read.table(exec.filepath)$V1
#    }
#    return(exec.time)
#}



#################################################################
# Reads the execution time result of a given partitioning method. This execution time is not the whole time spent by the CoNS(r_{max}).
#   Instead, we focus an the execution time of a specific value of r \in {1,2,.. r_{max}}, i.e. CoNS(r).
#
# n: number of vertices
# l0: number of modules
# d: graph density
# prop.mispl.int: proportion of misplaced internal links
# prop.mispl.ext: proportion of misplaced extenal links
# prop.neg: proportion of negative edges
# algo.name: a partitioning algo name
# nb.edit.of.interest: a specific value r of nb edit, where r <= r_{max}. 
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
read.exec.time.result = function(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, nb.edit.of.interest, force)
{
	exec.time = NA
    part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, NA)
    exec.filepath = file.path(part.folder, paste0("execTimes_nbEdit",nb.edit.of.interest,".txt"))
    if(file.exists(exec.filepath)){
        df = read.table(exec.filepath,sep=",")
        exec.time = as.numeric(unlist(strsplit(df[,2],":"))[2])
    }
    return(exec.time)
}



#################################################################
# Retrieves the execution times of the following methods for the given values of the input parameters:
#   - CoNSCC(r=3)
#   - CoNSCC(r=3) without MVMO pruning
#   - CoNSCC(r=4)
#   - CoNSCC(r=4) without MVMO pruning

# graph.sizes: number of vertices
# d: graph density
# l0.values: number of modules
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
##################################################################
collect.all.benchmark.results = function(graph.sizes, d, l0.values, force)
{
        col.names = c(paste(CoNSCC,"3edit"), paste(CoNSCC,"3edit without MVMO pruning"), paste(CoNSCC,"4edit"), paste(CoNSCC,"4edit without MVMO pruning"))
        data = c()
    	tlog("starts collecting the benchamrk results")
    	for(n in graph.sizes){
    	    for(l0 in l0.values){
    	        tlog(8, "fast jump model stats => n: ", n, ", l0: ",l0)
			
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
                        algo.name = get.CoNSCC.code(maxNbEdit=3, pruning.without.MVMO=FALSE)
			            res.CoNSCC3 = read.exec.time.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, 3, force)
			            #res.EnumCC3 = collect.all.benchmark.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, nb.repetitions, force, plot.formats)
                    
			            algo.name = get.CoNSCC.code(maxNbEdit=3, pruning.without.MVMO=TRUE)
			            res.CoNSCC3.pruning.without.MVMO = read.exec.time.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, 3, force)
			            #res.EnumCC3.brute.force = collect.all.benchmark.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, nb.repetitions, force, plot.formats)

						prop.neg = as.numeric(unlist(strsplit(dir2,"="))[2])
                        algo.name = get.CoNSCC.code(maxNbEdit=4, pruning.without.MVMO=FALSE)
			            res.CoNSCC4 = read.exec.time.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, 4, force)
			            #res.EnumCC4 = collect.all.benchmark.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, nb.repetitions, force, plot.formats)

			            algo.name = get.CoNSCC.code(maxNbEdit=4, pruning.without.MVMO=TRUE)
			            res.CoNSCC4.pruning.without.MVMO = read.exec.time.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, 4, force)
			            #res.EnumCC4.brute.force = collect.all.benchmark.result(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, algo.name, nb.repetitions, force, plot.formats)
		                line = matrix(c(res.CoNSCC3, res.CoNSCC3.pruning.without.MVMO, res.CoNSCC4, res.CoNSCC4.pruning.without.MVMO), nrow=1)

			            rownames(line) = paste0("n=",n,",l0=",l0,",propMisplInt=",prop.mispl.int,",propMisplExt=",prop.mispl.ext,",propNeg=",prop.neg)

			            data = rbind(data, line)
		            }
    	        }
    	    }
    	}
    	
    	colnames(data) = col.names
    	#print(data)

		if(length(data)>0){
			if(!dir.exists(OUTPUT.CSV.FOLDER))
				dir.create(OUTPUT.CSV.FOLDER, recursive=FALSE, showWarnings=FALSE)
			csv.path = file.path(OUTPUT.CSV.FOLDER,paste0("benchmark-results_d=",d,".csv"))
			write.csv(file=csv.path, x=data)
		}

}

