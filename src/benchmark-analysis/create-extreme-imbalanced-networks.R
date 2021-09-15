


#################################################################
#
# Retrieves the statistics regarding misplaced negative and positive links from a signed graph and its partition.
#   - prop.mispl.int
#   - imb.int.indxs
#   - prop.mispl.ext
#   - imb.ext.indxs
#
##################################################################
get.prop.imbalance.neg.pos.links.info = function(g, membership){
	el = get.edgelist(g)
	el = cbind(el, E(g)$weight)
	comembership = membership[as.integer(el[,1])+1] == membership[as.integer(el[,2])+1]
	
	neg.ext.indxs = which(!comembership & el[,3]<0)
	pos.ext.indxs = which(!comembership & el[,3]>0)
	nb.pos.ext = sum(as.integer(el[pos.ext.indxs,3]))
	nb.neg.ext = abs(sum(as.integer(el[neg.ext.indxs,3])))
	prop.pos.ext = (nb.pos.ext)/(nb.pos.ext+nb.neg.ext)
	#print(length(pos.ext.indxs))
	#print(prop.pos.ext)
	prop.mispl.ext = prop.pos.ext
	
	neg.int.indxs = which(comembership & el[,3]<0)
	pos.int.indxs = which(comembership & el[,3]>0)
	nb.pos.int = sum(as.integer(el[pos.int.indxs,3]))
	nb.neg.int = abs(sum(as.integer(el[neg.int.indxs,3])))
	prop.neg.int = (nb.neg.int)/(nb.pos.int+nb.neg.int)
	#print(length(neg.int.indxs))
	#print(prop.neg.int)
	prop.mispl.int = prop.neg.int
	
	result = list()
	result[["prop.mispl.int"]] = prop.mispl.int
	result[["imb.int.indxs"]] = neg.int.indxs
	result[["prop.mispl.ext"]] = prop.mispl.ext
	result[["imb.ext.indxs"]] = pos.ext.indxs
	return(result)
}

#################################################################
# It performs the following steps:
#   1) It creates three most perturbed perturbed signed graphs with respect to different link type strategy used in the perturbation process. 
#       - one graph, where only internal links are perturbed. Note that this graph is the most pertubed one with respect to internal links.
#       - one graph, where only external links are perturbed. Note that this graph is the most pertubed one with respect to internal links.
#       - one graph, where both internal and external links are perturbed. Note that this graph is the most pertubed one with respect to internal links.
#   2) Then, for each most perturbed grah descrived in step 1, it creates new signed graphs by decreasingly the amount of perturbation step by step.

#   Remark: the amount of perturbed edges amounts to be the underlying imbalance in the optimal solution, since the perturbation process starts with a perfectly balanced graph.
#
# n: number of vertices
# l0.values: number of modules
# d: graph density
# prop.mispl.int: proportion of misplaced internal links
# prop.mispl.ext: proportion of misplaced extenal links
# prop.neg: proportion of negative edges
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
# plot.formats: plot format(s), such as PDF, PNG or NA (it means 'no plotting')
#
##################################################################
create.extreme.imbalanced.networks = function(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, force, plot.formats)
{
	prop.mispl = 0
	graph.desc.name = "signed-unweighted"
	in.folder = get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg)
	graph.name = paste0(graph.desc.name,".G")
	network.path = file.path(in.folder,graph.name)
	
	print(network.path)
	if(file.exists(network.path)){
		g = read.graph.ils(network.path)
		
		algo.name = COR.CLU.ExCC
		part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, COR.CLU.ExCC)
		if (!dir.exists(part.folder))
			dir.create(path = part.folder, showWarnings = FALSE, recursive = TRUE)
		
		membership.filepath = file.path(part.folder, paste0(MBRSHP.FILE.PREFIX,0,".txt"))
		
		if(!file.exists(membership.filepath) || force){ 
			apply.partitioning.algorithm(part.folder, in.folder, algo.name, graph.name, plot.formats)
			# mbrshp = load.membership.files(part.folder)[[1]] # there will be a unique membership
			mbrshps = load.result.files(part.folder, algo.name)
			write.membership.files(part.folder, mbrshps)
			remove.result.files(part.folder)
		}
		membership = read.table(membership.filepath)$V1
		
		strengthed.lp.model.file.path = file.path(part.folder, "strengthedModelAfterRootRelaxation.lp")
		
		## =================================================
		#
		#strengthed.lp.model.file.path = file.path(part.folder, "strengthedModelAfterRootRelaxation.lp")
		#orig.strengthed.lp.model.file.path = file.path(part.folder, "strengthedModelAfterRootRelaxation_orig.lp")
		#
		#if(file.exists(strengthed.lp.model.file.path) && !file.exists(orig.strengthed.lp.model.file.path)){
		#    print(strengthed.lp.model.file.path)
		#    con <- file(strengthed.lp.model.file.path, "r")
		#    lines <- readLines(con)
		#    close(con)
		#
		#    # rename 
		#    file.rename(from=strengthed.lp.model.file.path, to=orig.strengthed.lp.model.file.path)
		#
		#
		#    n = vcount(g)
		#    A = as_adjacency_matrix(g, type = c("both"), attr="weight")
		#    A = as.matrix(A)
		#
		#    new.lines = perform.filtering.lp.miyauchi(A, n, lines)
		#    write.table(file=strengthed.lp.model.file.path, x=paste0(new.lines, collapse="\n"), row.names=FALSE, col.names=FALSE, quote=FALSE)
		#}
		## =================================================
		
		
		#fractional.graph.filepath = file.path(part.folder, "fractionalGraph.G")
		
		out.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, NA)
		
		
		# ------------------
		# Part 1
		# ------------------
		
		vars.filepath = file.path(out.folder,"vars.txt")
		bounds.filepath = file.path(out.folder,"b.txt")
		constrs.filepath = file.path(out.folder,"S.jld")
		
		if(!file.exists(vars.filepath) || !file.exists(bounds.filepath) || !file.exists(constrs.filepath) || force){
			cmd = "julia"
			cmd = paste(cmd, "src/benchmark-analysis/write-sparse-matrix-for-set-of-ineqs_extended.jl")
			cmd = paste(cmd, in.folder)
			cmd = paste(cmd, out.folder)
			cmd = paste(cmd, network.path)
			cmd = paste(cmd, membership.filepath)
			cmd = paste(cmd, strengthed.lp.model.file.path)
			#cmd = paste(cmd, fractional.graph.filepath)
			print(cmd)
			
			start = Sys.time()
			system(command = cmd)
			end = Sys.time()
			exec.time = as.numeric(end) - as.numeric(start)
			## save exec.time: write into file
			#write(x = exec.time, file = file.path(part.folder, EXEC.TIME.FILENAME))
		}
		
		
		# ------------------
		# Part 2
		# ------------------
		
		cplex.threads = 8
		
		perturbed.internal.graph.filepath  = file.path(out.folder,"perturbed_graph_internal.G")
		perturbed.external.graph.filepath  = file.path(out.folder,"perturbed_graph_external.G")
		perturbed.all.graph.filepath  = file.path(out.folder,"perturbed_graph_all.G")
		
		if(!file.exists(perturbed.internal.graph.filepath) || !file.exists(perturbed.external.graph.filepath) || !file.exists(perturbed.all.graph.filepath) || force){
			start = Sys.time()
			
			# internal
			if(!file.exists(perturbed.internal.graph.filepath) || force){
				cmd = "julia"
				cmd = paste(cmd, "src/benchmark-analysis/perturb-signed-graph.jl")
				#cmd = paste(cmd, "src/benchmark-analysis/perturb-signed-graph-with-resuming.jl") # it might be better to use this one
				                                                                                  #, since it can resume the perturbation process, whis has been already started before
				cmd = paste(cmd, out.folder)
				cmd = paste(cmd, network.path)
				cmd = paste(cmd, membership.filepath)
				cmd = paste(cmd, vars.filepath)
				cmd = paste(cmd, constrs.filepath)
				cmd = paste(cmd, bounds.filepath)
				cmd = paste(cmd, cplex.threads)
				cmd = paste(cmd, "internal")
				#cmd = paste(cmd, fractional.graph.filepath)
				print(cmd)
				system(command = cmd)
			}
			
			# external
			if(!file.exists(perturbed.external.graph.filepath) || force){
				cmd = "julia"
				cmd = paste(cmd, "src/benchmark-analysis/perturb-signed-graph.jl")
				#cmd = paste(cmd, "src/benchmark-analysis/perturb-signed-graph-with-resuming.jl") # it might be better to use this one
				                                                                                  #, since it can resume the perturbation process, whis has been already started before
				cmd = paste(cmd, out.folder)
				cmd = paste(cmd, network.path)
				cmd = paste(cmd, membership.filepath)
				cmd = paste(cmd, vars.filepath)
				cmd = paste(cmd, constrs.filepath)
				cmd = paste(cmd, bounds.filepath)
				cmd = paste(cmd, cplex.threads)
				cmd = paste(cmd, "external")
				#cmd = paste(cmd, fractional.graph.filepath)
				print(cmd)
				system(command = cmd)
			}
			
			# all
			if(!file.exists(perturbed.all.graph.filepath) || force){
				cmd = "julia"
				cmd = paste(cmd, "src/benchmark-analysis/perturb-signed-graph.jl")
				#cmd = paste(cmd, "src/benchmark-analysis/perturb-signed-graph-with-resuming.jl") # it might be better to use this one
				                                                                                  #, since it can resume the perturbation process, whis has been already started before
				cmd = paste(cmd, out.folder)
				cmd = paste(cmd, network.path)
				cmd = paste(cmd, membership.filepath)
				cmd = paste(cmd, vars.filepath)
				cmd = paste(cmd, constrs.filepath)
				cmd = paste(cmd, bounds.filepath)
				cmd = paste(cmd, cplex.threads)
				cmd = paste(cmd, "all")
				#cmd = paste(cmd, fractional.graph.filepath)
				print(cmd)
				system(command = cmd)
			}
			
			end = Sys.time()
			exec.time = as.numeric(end) - as.numeric(start)
		}
		
		
		# ------------------
		# Part 3
		# ------------------
		
		pert.new.graph.filepaths = c()
		#for(pert.graph.filepath in c(perturbed.internal.graph.filepath)){	
		for(pert.graph.filepath in c(perturbed.internal.graph.filepath, perturbed.external.graph.filepath, perturbed.all.graph.filepath)){
			#print(pert.graph.filepath)
			g = read.graph.ils(pert.graph.filepath)
			
			neg.indxs = which(E(g)$weight<0)
			pos.indxs = which(E(g)$weight>0)
			prop.neg = length(neg.indxs)/(length(neg.indxs)+length(pos.indxs))
			
			result = get.prop.imbalance.neg.pos.links.info(g, membership)
			prop.mispl.int = result[["prop.mispl.int"]]
			prop.mispl.ext = result[["prop.mispl.ext"]]
			
			out.folder = get.benchmark.input.network.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg)
			#print(out.folder)
			if (!dir.exists(out.folder))
				dir.create(path = out.folder, showWarnings = FALSE, recursive = TRUE)
			new.filepath = file.path(out.folder,"signed-unweighted.G")
			
			if(!file.exists(new.filepath)){
				file.copy(from=pert.graph.filepath, to=new.filepath, overwrite = TRUE, recursive = FALSE)
				write.graph(file=file.path(out.folder,"signed-unweighted.graphml"),graph=g,"graphml")
			}
			pert.new.graph.filepaths = c(pert.new.graph.filepaths, new.filepath) # for Part 4
		}
		
		print("part4")
		# ------------------
		# Part 4
		# ------------------
		for(pert.new.graph.filepath in pert.new.graph.filepaths){
			print("=============================================")
			cat("Based on the following perturbed graph: ", pert.new.graph.filepath,"\n")
			print(pert.new.graph.filepath)
			g = read.graph.ils(pert.new.graph.filepath)
			
			neg.indxs = which(E(g)$weight<0)
			pos.indxs = which(E(g)$weight>0)
			prop.neg = length(neg.indxs)/(length(neg.indxs)+length(pos.indxs))
			
			result = get.prop.imbalance.neg.pos.links.info(g, membership)
			prop.mispl.int = result[["prop.mispl.int"]]
			raw.imb.int.indxs = result[["imb.int.indxs"]] # all the indexes of the links contributing to the imalance inside the modules
			prop.mispl.ext = result[["prop.mispl.ext"]]
			raw.imb.ext.indxs = result[["imb.ext.indxs"]] # all the indexes of the links contributing to the imalance between the modules
			
			print(raw.imb.int.indxs)
			print(raw.imb.ext.indxs)
			imb.int.indxs = c() # init
			imb.ext.indxs = c() # init
			
			prop.step = 0.05
			int.props = seq(0,prop.mispl.int,prop.step)
			ext.props = seq(0,prop.mispl.ext,prop.step)
			
			# transform 'raw.imb.ext.indxs' into 'imb.ext.indxs' based on 'max.desired.prop.mispl.ext'
			# note that prop.mispl.ext > max.desired.prop.mispl.ext. This is why we construct 'imb.ext.indxs'
			ext.nb.interval = length(ext.props)
			max.desired.prop.mispl.ext = ext.props[ext.nb.interval] # last item
			if(length(raw.imb.ext.indxs)>0){
				rel.prop = (1*max.desired.prop.mispl.ext)/prop.mispl.ext
				raw.target.nb.links = length(raw.imb.ext.indxs)*rel.prop
				target.nb.links = raw.target.nb.links - (raw.target.nb.links %% ext.nb.interval)
				indexs = sample(x=seq(1,length(raw.imb.ext.indxs)),size=target.nb.links)
				imb.ext.indxs = raw.imb.ext.indxs[indexs]
			}
			
			# transform 'raw.imb.int.indxs' into 'imb.int.indxs' based on 'max.desired.prop.mispl.int'
			# note that prop.mispl.int > max.desired.prop.mispl.int. This is why we construct 'imb.int.indxs'
			int.nb.interval = length(int.props)
			max.desired.prop.mispl.int = int.props[int.nb.interval]
			if(length(raw.imb.int.indxs)>0){
				rel.prop = (1*max.desired.prop.mispl.int)/prop.mispl.int
				raw.target.nb.links = length(raw.imb.int.indxs)*rel.prop
				target.nb.links = raw.target.nb.links - (raw.target.nb.links %% int.nb.interval)
				indexs = sample(x=seq(1,length(raw.imb.int.indxs)),size=target.nb.links)
				imb.int.indxs = raw.imb.int.indxs[indexs]
			}
			
			# -------------------------------------------------------------
			# At every iteration, we start the pertrubed graph of the previous iteration and we make the same number of links balanced
			g2 = g
			int.target.nb.link.step = length(imb.int.indxs)/int.nb.interval
			ext.target.nb.link.step = length(imb.ext.indxs)/ext.nb.interval
			for(int.i in 1:int.nb.interval){
				for(ext.i in 1:ext.nb.interval){
					if(length(imb.int.indxs)>0 || length(imb.ext.indxs)>0){
						
						if(length(imb.int.indxs)>0){
							indexs = sample(x=seq(1,length(imb.int.indxs)),size=int.target.nb.link.step)
							target.indexs = imb.int.indxs[indexs]
							E(g2)$weight[target.indexs] = -E(g2)$weight[target.indexs] # inverse the sign
							imb.int.indxs = imb.int.indxs[-indexs]
						}
						
						if(length(imb.ext.indxs)>0){
							indexs = sample(x=seq(1,length(imb.ext.indxs)),size=ext.target.nb.link.step)
							target.indexs = imb.ext.indxs[indexs]
							E(g2)$weight[target.indexs] = -E(g2)$weight[target.indexs] # inverse the sign
							imb.ext.indxs = imb.ext.indxs[-indexs]
						}
						
						info = get.prop.imbalance.neg.pos.links.info(g2, membership)
						neg.indxs = which(E(g2)$weight<0)
						pos.indxs = which(E(g2)$weight>0)
						prop.neg2 = length(neg.indxs)/(length(neg.indxs)+length(pos.indxs))
						cat("int.i: ",int.i, ", ext.i: ",ext.i,", prop.mispl.int: ", info[["prop.mispl.int"]], ", prop.mispl.ext:", info[["prop.mispl.ext"]],"prop.neg: ",prop.neg2,"\n")
						
						
						out.folder = get.benchmark.input.network.folder.path(n, l0, d, info[["prop.mispl.int"]], info[["prop.mispl.ext"]], prop.neg2)
						#print(out.folder)
						if (!dir.exists(out.folder))
							dir.create(path = out.folder, showWarnings = FALSE, recursive = TRUE)
						new.filepath = file.path(out.folder,"signed-unweighted.G")
						new.filepath2 = file.path(out.folder,"signed-unweighted.graphml")
						write.graph.ils(graph=g2, new.filepath)
						write.graph(graph=g2, new.filepath2, "graphml")
						
						
						overall.part.folder = get.benchmark.part.folder.path(n, l0, d, info[["prop.mispl.int"]], info[["prop.mispl.ext"]], prop.neg2, NA)
						#print(overall.part.folder)
						if (!dir.exists(overall.part.folder))
							dir.create(path = overall.part.folder, showWarnings = FALSE, recursive = TRUE)
						file.copy(from=part.folder, to=overall.part.folder, recursive=TRUE)
					}
				}
				
			}
			# -------------------------------------------------------------

			# copy the ExCC folder into the folder associated with the inital perturbed graph 
			overall.part.folder = get.benchmark.part.folder.path(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, NA)
			#print(overall.part.folder)
			if (!dir.exists(overall.part.folder))
				dir.create(path = overall.part.folder, showWarnings = FALSE, recursive = TRUE)
			file.copy(from=part.folder, to=overall.part.folder, recursive=TRUE)
			
		}
		
	}
}



#################################################################
# This is the starter method for perturbing perfectly balanced signed networks.
#
# graph.sizes: number of vertices
# d: graph density
# l0.values: number of modules
# prop.negs: a set of values from the range [0,1] for the proportion of negative edges
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
# plot.formats: plot format(s), such as PDF, PNG or NA (it means 'no plotting')
#
##################################################################
create.all.extreme.imbalanced.networks = function(graph.sizes, d, l0.values, prop.negs, force, plot.formats)
{
	prop.mispl.int = 0
	prop.mispl.ext = 0
	
	tlog("starts fast jump model stats")
	for(n in graph.sizes){
		for(l0 in l0.values){
			tlog(8, "fast jump model stats => n: ", n, ", l0: ",l0)
			
			my.prop.negs = prop.negs # if we do not do that, for each n value, prop negs will not be the initial value(s)
			if(is.na(my.prop.negs) && d == 1){
				prop.mispl = 0
				my.prop.negs = compute.prop.neg(n, d, l0, prop.mispl)
			}
			
			for (prop.neg in my.prop.negs) {
				tlog(12, "fast jump model stats => prop.neg: ", prop.neg)
				
				create.extreme.imbalanced.networks(n, l0, d, prop.mispl.int, prop.mispl.ext, prop.neg, force, plot.formats)
			}
		}
	}
}
