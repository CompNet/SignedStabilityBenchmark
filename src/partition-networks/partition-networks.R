

#############################################################################################
# It loads all partitions from the raw result files whose the name is structured as "sol<XX>.txt" where XX
#   indicates the partition number (e.g. 0, 1, etc.), it converts the raw result into vector, 
#   and it finally stores as a list of vectors. Note that each partition info, i.e. membership, is stored as a vector.
#   Note that the files "sol0.txt" might be like this (each line corresponds to a cluster):
#     [0,2,4,5,9] ----> cluster 1
#     [1,3,6,7,8] ----> cluster 2
#         .
#         .
#
# part.folder: the partitioning result folder
# algo.name: the name of the partitioning algorithm which produced the results
#############################################################################################
load.result.files = function(part.folder, algo.name) {
    mbrshps = list()
    result.files = list.files(path = part.folder,
                              pattern = paste0("^", ALGO.RESULT.FILE.PREFIX, ".*\\.txt$"))
    nb.result.file = length(result.files)
    
    for (i in 1:nb.result.file) {
        algo.result.filename = result.files[i]
		print(algo.result.filename)
        # load the resulting partition file
        mbrshp <-
            load.external.partition(part.folder, algo.result.filename, algo.name, keep.tmp =
                                        TRUE)
		print(mbrshp)
        mbrshps[[i]] = mbrshp
    }
    
    return(mbrshps)
}




#############################################################################################
# It loads all partitions from the membership files whose the name is structured as "membership<XX>.txt" where XX
#   indicates the partition number (e.g. 0, 1, etc.), and then it stores as a list of vectors. 
#   Note that each partition info, i.e. membership, is stored as a vector.
#   Note that a file "membership0.txt" might be like this (each line is associated with a node and each line contains cluster info):
#     1
#     2
#     1
#     3
#     .
#     .
#
# part.folder: the partitioning result folder
#############################################################################################
load.membership.files = function(part.folder) {
    mbrshps = list()
    #paste0("^membership",".*\\.txt$")
    mbrshp.files = list.files(path = part.folder,
                              pattern = paste0("^", MBRSHP.FILE.PREFIX, ".*\\.txt$"))
    nb.mbrshp.file = length(mbrshp.files)
    
    if (nb.mbrshp.file > 0) {
        for (id in 0:(nb.mbrshp.file - 1)) {
            # load the resulting partition file
            table.file = file.path(part.folder, paste0(MBRSHP.FILE.PREFIX, id, ".txt"))
            mbrshp <-
                as.numeric(as.matrix(read.table(
                    file = table.file, header = FALSE
                )))
            mbrshps[[id + 1]] = mbrshp
        }
    }
    
    return(mbrshps)
}



#############################################################################################
# It removes/deletes the raw (partition) result files
#
# part.folder: the partitioning result folder
#############################################################################################
remove.result.files = function(part.folder) {
    result.files = list.files(path = part.folder,
                              pattern = paste0("^", ALGO.RESULT.FILE.PREFIX, ".*\\.txt$"))
    unlink(x = file.path(part.folder, result.files))
}



#############################################################################################
# It removes/deletes the log files created during ExCC execution. Specifically,
#   - 'logcplex.txt': the log file from ExCC execution containing some infos related to
#       evolution of Gap, Number of branch nodes, Best Ineger, etc.
#   - 'log.txt': it records the number of valid cuts before enumrating all solutions. 
#       Note that finding valid cuts make the algorithm faster
#
# part.folder: the partitioning result folder
#############################################################################################
remove.log.files = function(part.folder) {
    log.files = c("log.txt","logcplex.txt")
    unlink(x = file.path(part.folder,log.files))
}





#############################################################################################
# It writes all partition results into file in 'membership' format
#   Note that a file "membership0.txt" might be like this (each line is associated with a node and each line contains cluster info):
#     1
#     2
#     1
#     3
#     .
#     .
#
# part.folder: the partitioning result folder
# mbrshps: all partition results as a list of vectors
#############################################################################################
write.membership.files = function(part.folder, mbrshps) {
    n = length(mbrshps)
    for (id in 0:(n - 1)) {
        table.file = file.path(part.folder, paste0(MBRSHP.FILE.PREFIX, id, ".txt"))
        write.table(
            x = mbrshps[[id + 1]],
            file = table.file,
            row.names = FALSE,
            col.names = FALSE
        )
    }
}


#############################################################################################
# It checks if the R script has launched the partitioning method for optimal solutions.
# This method should be used by combining with the method 'is.membership.file.created()'
#   If the partitioning method had completed the execution of the method, it means that there are some raw result files. 
#
# algo.name: the name of the partitioning algorithm which produced the results
# part.folder: the partitioning result folder
#############################################################################################
is.raw.result.file.created = function(algo.name, part.folder) {
    files = list.files(path = part.folder,
                       pattern = paste0("^", ALGO.RESULT.FILE.PREFIX, "*"))
    nb.raw.result.file = length(files)
    return(nb.raw.result.file != 0)
}



#############################################################################################
# It checks if the R script has processed the raw partition result files:
#   the process of the raw partition result files consists of converting from 'solXX.txt' into'membershipXX.txt'.
#   If 'solXX.txt' files exist, it means that the R script has not processed the raw partition result files.
#
# part.folder: the partitioning result folder
#############################################################################################
is.membership.file.created = function(part.folder) {
    files = list.files(path = part.folder, pattern = paste0("^", MBRSHP.FILE.PREFIX, "*"))
    nb.mbrshp.file = length(files)
    return(nb.mbrshp.file != 0)
}


#############################################################################################
# It executes the partitioning method through command line, since the method is written in Java.
# The whole workflow is described as follows:
# - it launches the bash script which tracks the used memory during Cplex. This script will be stopped automatically
#     when it can not find any Cplex processes
# - it executes the partitioning method
# - it writes the total execution time into file
# - it updates the name of the resulting file if the aim is not to generate all optimal solutions, but only one
#    The reason is that when the method is run for only for one optimal solution, it will be named as 'sol0.txt by default.
#     To prevent from having the '0' suffix, the renaming is established.
#
# part.folder: 
# in.folder: 
# algo.name: 
# g.name:
# plot.formats:
#############################################################################################
apply.partitioning.algorithm = function(part.folder, in.folder, algo.name, g.name, plot.formats) {
    
    # ----------------------------------
    ## start to record used memory during Cplex in infinite loop until the command has been terminated - asynchrone command
    #system(
    #    command = paste("bash", file.path(LIB.FOLDER, RECORD.MEM.INFO.SCRIPTNAME), part.folder),
    #    wait = FALSE
    #    ) # save the output in the corresponding part folder
    ## ----------------------------------
    

    # apply the correlation clustering algorithm (which are external programs)
    tlog(32, "Applying ", algo.name, " to the network")
    # set the external command and invoke it
    cmd <-
        get.algo.commands(
            algo.names = algo.name,
            input.folder = in.folder,
            out.folder = part.folder,
            graph.name = g.name
        )
    tlog(32, "Command: ", cmd)
    if(cmd != "NONE"){
        start = Sys.time()
        system(command = cmd)
        end = Sys.time()
        exec.time = as.numeric(end) - as.numeric(start)
        # save exec.time: write into file
        write(x = exec.time, file = file.path(part.folder, EXEC.TIME.FILENAME))
        prepare.algo.output.filename(part.folder, algo.name, g.name)
    }
}



#############################################################################################
# It partitions the specified network, using the specified algorithm, and records the result
# as a table. Optionnaly, it can plot the network with the partition information: as many network plots as obtained solutions
#
# g: graph to process.
# algo.name: (normalized) name of the community detection algorithm.
# part.folder: folder in which to write the result (partition) files.
# graph.folder: folder of the processed network (for external tools).
# plot.format: plot format(s), such as PDF, PNG or NA (it means 'no plotting')
# plot.layout: plot layout, such as "kamada.kawai", "fruchterman.reingold", "bn_zheng" or "circle"
# force: indicates whether existing result files should be loaded and used vs. replaced by new ones.
#
# returns:
#############################################################################################
process.partitioning.algorithm <- function(part.folder, in.folder, algo.name, graph.name, keep.algo.log.files, 
                                           plot.format, plot.layout, force)
    {
        #tlog("n=",vcount(g), " m=",ecount(g), " d=",graph.density(g))
        #tlog(" connected=", is.connected(g,mode="weak"))
        
        # check if algo result files already exist: if the result files are not created, process the algo
        process.algo <- !is.raw.result.file.created(algo.name, part.folder)
        process.mbrshp <- !is.membership.file.created(part.folder)
        
        mbrshps = NA
        if ((process.algo && process.mbrshp) || force) {
            tlog(24, "Applying algorithm ", algo.name, " on folder ", part.folder)
            
            apply.partitioning.algorithm(part.folder, in.folder, algo.name, graph.name, plot.formats)
            
            #if(!startsWith(algo.name,ENUM.POPULATE.ExCC)){
            #    mbrshps = load.result.files(part.folder, algo.name)
            #
            #    write.membership.files(part.folder, mbrshps)
            #    remove.result.files(part.folder)
            #    
            #    if(!keep.algo.log.files)
            #        remove.log.files(part.folder)
            #}
		
			mbrshps = load.result.files(part.folder, algo.name)
			write.membership.files(part.folder, mbrshps)
			remove.result.files(part.folder)

            mbrshps = load.membership.files(part.folder)
			print(mbrshps)
        }	
        else if (!process.algo && process.mbrshp) {
            tlog(
                24, "Raw algo result file is already present for algorithm ",
                algo.name, " on folder ", part.folder, 
                " so we just load the existing result, and write it into membership file"
            )
            
            mbrshps = load.result.files(part.folder, algo.name)
            write.membership.files(part.folder, mbrshps)
            remove.result.files(part.folder)

            mbrshps = load.membership.files(part.folder)
        } else{
            tlog(
                24, "Membership file is already present for algorithm ",
                algo.name, " on folder ", part.folder, " so we just load the existing result"
            )
            mbrshps = load.membership.files(part.folder)
        }
        
        
        # # ---------------------------------------------------------------------------
        # # plot: if plot.formats = NA, then no plotting
        # network.path = file.path(in.folder,graph.name)
        # g = read.graph.ils(network.path)
        # 
        # nb.mbrshp.file = length(mbrshps)
        # for(id in seq(0,nb.mbrshp.file-1)){
        # 	mbrshp = mbrshps[[id+1]]
        # 	# print(mbrshp)
        # 	plot.file = file.path(part.folder, paste0(MBRSHP.FILE.PREFIX,id))
        # 	plot.network(g, membership=mbrshp, plot.file, format=plot.format, method=plot.layout)
        # }
        # # ---------------------------------------------------------------------------
        
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
partition.network = function(n, l0, d, prop.mispl, prop.neg, network.no,
                             cor.clu.heur.algos, heur.reps, cor.clu.exact.algos, exact.reps, 
                             keep.algo.log.files, plot.format, plot.layout, force)
{
    net.folder = get.input.network.folder.path(n, l0, d, prop.mispl, prop.neg, network.no)
    tlog(16, "start to partition networks with exact algorithms")

    for (graph.desc.name in c(SIGNED.UNWEIGHTED.FILE)) {
        tlog(20, "partitioning networks => graph.desc.name: ", graph.desc.name)

        tlog(20, "start to partition networks with exact algorithms")	
        for(algo.name in cor.clu.exact.algos){
            tlog(24, "start to partition networks with the exact algorithm: ", algo.name)

            for(rep.no in exact.reps){	
                part.folder = get.part.folder.path(n, l0, d, prop.mispl,
                                prop.neg, network.no, algo.name, graph.desc.name, rep.no)
                if (dir.exists(net.folder)) {
                    #if (dir.exists(part.folder)){
                    #    dir.create(paste0(part.folder,"/../BACKUP"))
                    #    file.copy(part.folder, paste0(part.folder,"/../BACKUP"), recursive=TRUE)
                    #}
                    if (!dir.exists(part.folder))
                        dir.create(path = part.folder, showWarnings = FALSE, recursive = TRUE)
                    
                    #if(!file.exists(file.path(part.folder,"..","..",get.ExCC.code(enum.all=FALSE),"signed-unweighted","strengthedModelAfterRootRelaxation.lp"))){
                    #    print(file.path(part.folder,"..","..",get.ExCC.code(enum.all=FALSE),"signed-unweighted","strengthedModelAfterRootRelaxation.lp"))
                    #    print(part.folder)
                    #    print("!!!!!!!!!!!!!!!!!")
                    #    unlink(part.folder, recursive=TRUE)
                    #}

                    
                    tlog(28, "partitioning networks in ", part.folder)
                    graph.name = paste0(graph.desc.name, ".G")
                    process.partitioning.algorithm(part.folder, net.folder, algo.name, graph.name, keep.algo.log.files,
                                                   plot.format, plot.layout, force)
                }
            
            }
            
        }
        
        
        tlog(20, "start to partition networks with heuristic algorithms")	
        for(algo.name in cor.clu.heur.algos){
            tlog(24, "start to partition networks with the heuristic algorithm: ", algo.name)
            
            for(rep.no in heur.reps){	
                part.folder = get.part.folder.path(n, l0, d, prop.mispl,
                                                   prop.neg, network.no, algo.name, graph.desc.name, rep.no)

                if (dir.exists(net.folder)) {
                    if (!dir.exists(part.folder))
                        dir.create(path = part.folder, showWarnings = FALSE, recursive = TRUE)

                    tlog(28, "partitioning networks in ", part.folder)
                    graph.name = paste0(graph.desc.name, ".G")
                    process.partitioning.algorithm(part.folder, net.folder, algo.name, graph.name, keep.algo.log.files,
                                                   plot.format, plot.layout, force)
                }
                
            }
            
        }
        
        
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
partition.networks = function(graph.sizes, d, l0, prop.mispls, prop.negs, in.rand.net.folders,
                              cor.clu.heur.algos, heur.reps, cor.clu.exact.algos, exact.reps, 
                              keep.algo.log.files, plot.format, plot.layout, force)
{
    tlog("starts partitioning networks")
    for (n in graph.sizes) {
        tlog(4, "partitioning networks => n: ", n)
        
        for (prop.mispl in prop.mispls) {
            tlog(8, "partitioning networks => prop.mispl: ", prop.mispl)
            
            my.prop.negs = prop.negs # if we do not do that, for each n value, prop negs will not be the initial value(s)
            if(is.na(my.prop.negs) && d == 1){
                my.prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
            }
            
            for (prop.neg in my.prop.negs) {
                tlog(12, "partitioning networks => prop.neg: ", prop.neg)
                
                for (network.no in in.rand.net.folders) {
                    tlog(16, "partitioning networks => network.no: ", network.no)
                    
                    partition.network(n, l0, d, prop.mispl, prop.neg, network.no, 
                                      cor.clu.heur.algos, heur.reps, cor.clu.exact.algos, exact.reps, 
                                      keep.algo.log.files, plot.format, plot.layout, force)
                }
                
            }
            
        }
        
    }
}
