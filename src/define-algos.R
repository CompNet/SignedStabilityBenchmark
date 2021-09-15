#################################################################
# Correlation Clustering (CC) problem
#################################################################


# ===============================================================
# Exact Approach: ExCC
# ===============================================================

NB.THREAD = 5

COR.CLU.ExCC <- "ExCC"
COR.CLU.ExCC.ENUM.ALL <- "ExCC-all" # TODO
ExCC.LIB.FOLDER = file.path(LIB.FOLDER,COR.CLU.ExCC)
ExCC.JAR.PATH = paste(ExCC.LIB.FOLDER,"ExCC.jar",sep="/") # gaia cluster - CERI
#CPLEX.BIN.PATH = "/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/"
CPLEX.BIN.PATH = "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/"
ExCC.MAX.TIME.LIMIT = 3600



# ===============================================================
# Exact Approach: EnumCC
# ===============================================================

ENUMCC = "EnumCC"
ENUMCC.LIB.FOLDER = file.path(LIB.FOLDER,ENUMCC)
ENUMCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0(ENUMCC,".jar"),sep="/") # gaia cluster - CERI
ENUMCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0("EnumCC.jar"),sep="/") # gaia cluster - CERI
#ENUM.POPULATE.CC.JAR.PATH = paste(ENUM.POPULATE.CC.LIB.FOLDER,paste0("MyPopulateCC_CplexFast_NEW.jar"),sep="/") # gaia cluster - CERI

ENUMCC.MAX.NB.SOLS = 50000 #10000 # we know that our method is not efficient for very large number of solutions
ENUMCC.MAX.TIME.LIMIT = 3600*12 # 12 hours


# ===============================================================
# Exact Approach: RNSCC (which is part of EnumCC)
# ===============================================================

CoNSCC = "CoNS"
RNSCC = "RNSCC"
RNSCC.JAR.PATH = paste(ENUMCC.LIB.FOLDER,paste0(RNSCC,".jar"),sep="/")





############################################################################
# It reads a .G graph file, and returns the contents as a data frame object.
#
# network.path: the file path which stores the .G graph file
#
############################################################################
read.graph.ils.file.as.df = function(network.path){
	# skip the first line bc it does not contain graph info
	df = read.table(
			file=network.path, 
			header=FALSE, 
			sep="\t", 
			skip=1, 
			check.names=FALSE
	)
	# df$V1: vertex1
	# df$V2: vertex2
	# df$V3: weight
	return(df)
}



############################################################################
#  It reads a .G graph file, and returns the contents as a igraph graph object.
#  To handle isolated nodes, first we had to find the max vertex id.
#  Then, indicate explicitely vertices ids in graph.data.frame()
#
# network.path: the file path which stores the .G graph file
#
############################################################################
read.graph.ils = function(network.path){
	df = read.graph.ils.file.as.df(network.path)
	
	edg.list = df[,c(1, 2)]
	max.v.id = max(unique(c(edg.list[,1], edg.list[,2])))
	
	g <- graph.data.frame(edg.list, vertices=seq(0,max.v.id), directed=FALSE)
	cat("max id: ",max.v.id, "\n")
	E(g)$weight = df[, 3]
	# V(g)$id = seq(0,max.v.id)
	# V(g)$id = seq(1,max.v.id+1)
	
	return(g)
}


############################################################################
# It writes the graph object into a file
#
# graph: the igraph graph object
# file path: the file path which will store the graph content in the format of .G graph file
############################################################################
write.graph.ils = function(graph, file.path){
# export using a format compatible with pILS
	t <-get.edgelist(graph=graph)
	t =  matrix(as.integer(t), nrow(t), ncol(t))
	if(t[1,1] == 1)
		t <- t - 1	# start numbering nodes at zero
	
	t <- cbind(t,E(graph)$weight)		# add weights as the third column
	
	write.table(data.frame(vcount(graph),ecount(graph)), file=file.path, append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) # write header
	write.table(t, file=file.path, append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE) # write proper graph
}









#############################################################################################
# 
#############################################################################################
get.ExCC.code <- function(enum.all)
{
	result <- COR.CLU.ExCC
	if(enum.all)
		result <- COR.CLU.ExCC.ENUM.ALL
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.ExCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	is.cp = "true" # in any case, use cutting plane approach
	is.enum.all = "false"
	# tilim = 3600 # 1 hour
	#if(algo.name == COR.CLU.ExCC.ENUM.ALL){
	#	is.enum.all = "true"
	#}
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	input.file.for.g = file.path(input.folder, graph.name)
	g = read.graph.ils(input.file.for.g)
	
	initSolutionFilePath = file.path(out.folder,"..","..",COR.CLU.ExCC,GRAPH.FILENAME,"membership0.txt")
	if(!file.exists(initSolutionFilePath))
		initSolutionFilePath="''"
	
	cmd = "NONE"
	if(is.enum.all == "false"){
		# An example:
		# java -Djava.library.path=/users/narinik/Cplex/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux/
		# -DinFile="in/""$name" -DoutDir="out/""$modifiedName" -DenumAll=false -Dcp=true -DMaxTimeForRelaxationImprovement=20
		# -DuserCutInBB=false -DinitSolutionFilePath="$initSolutionFilePath" -DLPFilePath="$LPFilePath"
		# -DonlyFractionalSolution=false -DfractionalSolutionGapPropValue=-1.0 -DnbThread=2 -Dverbose=true -Dtilim=200 -jar exe/ExCC.jar
		
		# -------------------------------------------------------------------------
		# TODO handle this in a better way, for instance, use small values for sparse networks
		# TODO write a function called 'estimate.max.time.for.relaxation.improvement(..)'
		maxTimeForRelaxationImprovement = "600"
		if(vcount(g)>39 && vcount(g)<=50)
			maxTimeForRelaxationImprovement = "4500" # 1h15m
		else if(vcount(g)>50)
			maxTimeForRelaxationImprovement = "10000" # 166 mins 
		# -------------------------------------------------------------------------
		
		cmd = 
			paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DinFile=", input.file),
				paste0("-DoutDir=", out.folder),
				paste0("-Dcp=",is.cp),
				paste0("-DenumAll=",is.enum.all),
				paste0("-Dtilim=",ExCC.MAX.TIME.LIMIT),
				paste0("-DtilimForEnumAll=",-1),
				paste0("-DsolLim=",1),
				paste0("-DMaxTimeForRelaxationImprovement=",maxTimeForRelaxationImprovement),
				"-DlazyInBB=false",
				"-DuserCutInBB=false",
				paste0("-DinitSolutionFilePath=",initSolutionFilePath),
				"-Dverbose=true",
				paste0("-DnbThread=", NB.THREAD),
				"-DLPFilePath=''",
				"-DonlyFractionalSolution=false",
				"-DfractionalSolutionGapPropValue=-1",
				"-jar",
				ExCC.JAR.PATH,
				sep=" "
			)
	} else { # if(is.enum.all == "true")
		# An example:
		# java -DinFile="in/""$name" -DoutDir="out/""$modifiedName" -DenumAll=true
		# -Dcp=false -DinitSolutionFilePath="$initSolutionFilePath" -DLPFilePath="$LPFilePath"
		# -DnbThread=2 -Dverbose=true -Dtilim=-1 -DtilimForEnumAll=60 -DsolLim=100 -jar exe/ExCC.jar
		
#		LP.filepath = file.path(out.folder,"..","..",COR.CLU.ExCC,GRAPH.FILENAME,"strengthedModelAfterRootRelaxation.lp")
#		if(file.exists(LP.filepath)) {
#			cmd = 
#				paste(
#					"java",		
#					paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
#					paste0("-DinFile=", input.file),
#					paste0("-DoutDir=", out.folder),
#					paste0("-Dcp=","false"),
#					paste0("-DenumAll=","true"),
#					paste0("-Dtilim=",-1),
#					paste0("-DtilimForEnumAll=",ExCCAll.MAX.TIME.LIMIT),
#					paste0("-DsolLim=",ExCCAll.MAX.NB.SOLS),
#					paste0("-DMaxTimeForRelaxationImprovement=","-1"), # no specific time limit, use the default one
#					"-DlazyInBB=false",
#					"-DuserCutInBB=false",
#					paste0("-DinitSolutionFilePath=",initSolutionFilePath),
#					"-Dverbose=true",
#					paste0("-DnbThread=", NB.THREAD),
#					paste0("-DLPFilePath='",LP.filepath,"'"),
#					"-DonlyFractionalSolution=false",
#					"-DfractionalSolutionGapPropValue=-1",
#					"-jar",
#					ExCC.JAR.PATH,
#					sep=" "
#				)
#		}
	}

	print(cmd)
	return(cmd)
}


#############################################################################################
# 
#############################################################################################
get.EnumCC.code <- function(maxNbEdit)
{
	result <- paste0(ENUMCC,"-maxNbEdit",maxNbEdit)
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.EnumCC.command <- function(algo.name, input.folder, out.folder, graph.name)
{
	print(algo.name)
	base.algo.name <- strsplit(x=algo.name, split="-", fixed=TRUE)[[1]][1]
	params.str <- gsub(paste0(base.algo.name,"-"),"",algo.name)
	print(params.str)
	
	print(graph.name)
	input.file = paste("'", input.folder, "/", graph.name, "'", sep="")
	
	maxNbEdit = as.integer(gsub("maxNbEdit","",params.str))
	
	cmd = "NONE"
	LP.filepath = file.path(out.folder,"..","..",COR.CLU.ExCC,GRAPH.FILENAME,"strengthedModelAfterRootRelaxation.lp")
	
	if(file.exists(LP.filepath)) {
		cmd = 
			paste(
				"java",		
				paste("-Djava.library.path=", CPLEX.BIN.PATH, sep=""),
				paste0("-DinFile=", input.file),
				paste0("-DoutDir=", out.folder),
				paste0("-DLPFilePath=", LP.filepath),
				paste0("-DinitMembershipFilePath=", file.path(out.folder,"..","..",COR.CLU.ExCC,GRAPH.FILENAME,"membership0.txt")),
				paste0("-DnbThread=",NB.THREAD),
				paste0("-DmaxNbEdit=",maxNbEdit),
				paste0("-DsolLim=",ENUMCC.MAX.NB.SOLS),
				paste0("-Dtilim=",ENUMCC.MAX.TIME.LIMIT),
				paste0("-DJAR_filepath_EnumCC=",paste(ENUMCC.LIB.FOLDER,paste0("EnumPopulateCCOnePass.jar"),sep="/")),
				"-jar",
				ENUMCC.JAR.PATH,
				sep=" "
			)
	}
	
	print(cmd)
	return(cmd)
}


#############################################################################################
# 
#############################################################################################
get.CoNSCC.code <- function(maxNbEdit, pruning.without.MVMO)
{
	result <- paste0(CoNSCC,"-maxNbEdit",maxNbEdit)
	if(pruning.without.MVMO)
		result <- paste0(result,"-pruningWithoutMVMO")
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.RNSCC.code <- function(maxNbEdit, pruning.without.MVMO)
{
	result <- paste0(RNSCC,"-maxNbEdit",maxNbEdit)
	if(pruning.without.MVMO)
		result <- paste0(result,"-pruningWithoutMVMO")
	return(result)
}


#############################################################################################
# 
#############################################################################################
get.RNSCC.command <- function(algo.name, sol.lim, input.folder, out.folder, graph.name)
{
	print(algo.name)
	base.algo.name <- strsplit(x=algo.name, split="-", fixed=TRUE)[[1]][1]
	params.str <- gsub(paste0(base.algo.name,"-"),"",algo.name)
	print(params.str)
	
	pruningWithoutMVMO = "false"
	maxNbEdit = 1
	
	params.str.list <- unlist(strsplit(x=params.str, split="-", fixed=TRUE))
	if(length(params.str.list) == 1)
		maxNbEdit = as.integer(gsub("maxNbEdit","",params.str))
	else #  length(params.str.list) == 2
	{
		maxNbEdit = as.integer(gsub("maxNbEdit","",params.str.list[1]))
		pruningWithoutMVMO = "true"
	}
	
	graph.name = paste0("signed-unweighted", ".G")
	network.path = file.path(input.folder,graph.name)
	
	cmd = 
		paste(
			"java",
			paste0("-DinitMembershipFilePath=", file.path(out.folder,paste0(MBRSHP.FILE.PREFIX,"0.txt"))),
			paste0("-DallPreviousResultsFilePath=", file.path(out.folder,"allResults.txt")), # it should be an empty file at startup
			paste0("-DinputFilePath=", network.path),
			paste0("-DoutDir=", out.folder),
			paste0("-DmaxNbEdit=", maxNbEdit),
			"-Dtilim=-1",
			paste0("-DsolLim=",sol.lim),
			"-DnbThread=1",
			paste0("-DisBruteForce=",pruningWithoutMVMO),
			"-DisIncrementalEditBFS=true",
			"-jar",
			RNSCC.JAR.PATH,
			sep=" "
		)
	
	print(cmd)
	return(cmd)
}


#############################################################################################
#
#############################################################################################
prepare.algo.output.filename = function(part.folder, algo.name, g.name){
    
    if(algo.name == COR.CLU.ExCC)
    {
#        ExCC.output.file <- file.path(part.folder, "ExCC-result.txt")
#        id=0
#        file.rename(from=ExCC.output.file, to=file.path(part.folder, paste0(ALGO.RESULT.FILE.PREFIX,id,".txt")))
    }
	else if(startsWith(algo.name,ENUMCC)){
		# do nothing
	}
	else if(startsWith(algo.name,RNSCC)){
		# do nothing
	}
    else {
        # TODO
        print("!!!!!!!!!! in TODO")
        # do nothing
    }

    
}




#############################################################################################
# Returns the full name based on the normalized (short) name. Note that for parameterized 
# algorithms, this will just return a clean version of the short name, since it contains 
# the parameter values.
#
# algo.names: short names of the considered algorithms.
#
# returns: the corresponding full names, to be used in plots for instance.
#############################################################################################
get.algo.names <- function(algo.names)
{	result <- c()
	
	for(algo.name in algo.names)
	{
		# parameters
		result <- c(result, gsub(pattern="_", replacement=" ", x=algo.name, fixed=TRUE))
	}
	
	return(result)
}



#############################################################################################
# Returns the inline command for the specified algorithm. The "..." parameters are fetched
# to the algorithm-specific function.
#
# algo.name: short code associated to the algorithm, including its parameter values.
#
# returns: the command allowing to invoke the program externally.
#############################################################################################
get.algo.commands <- function(algo.names, ...)
{	result <- c()

    # substring(x, 1, nchar(prefix)) == prefix
    
    for(algo.name in algo.names)
    {	
        if(startsWith(algo.name,COR.CLU.ExCC))
            result <- c(result, get.ExCC.command(algo.name, ...))
		else if(startsWith(algo.name,ENUMCC))
			result <- c(result, get.EnumCC.command(algo.name, ...))
		else if(startsWith(algo.name,RNSCC))
			result <- c(result, get.RNSCC.command(algo.name, sol.lim=-1, ...))
		else if(startsWith(algo.name,CoNSCC))
			result <- c(result, get.RNSCC.command(algo.name, sol.lim=2, ...))
    }
    
    return(result)
}
