
 
palette <- as.list(rainbow(14))
palette[1:4] = rainbow(4, v=0.7)
# ==============
palette[[5]] = "#E3CF57" # banana yellow
palette[[6]] = "maroon1" #
palette[[7]] = "darkviolet" # 
palette[[8]] = "plum1" #
palette[[9]] = "seagreen4" #
palette[[10]] = "black" # 
palette[[11]] = "burlywood4" # 
palette[[12]] = "yellow" #
palette[[13]] = "lightseagreen"
palette[[14]] = "gray87" # light gray
# ==============


#graph.sizes = c(30,40)
#d=0.25
#l0.values=c(3)
#nb.edit=3
#force=FALSE



############################################################################
#
# For each plot, the x-axis represents graph order n, and execution times (in seconds) 
# are shown in the log-scaled y-axis. A plot line is solid (resp. dashed)
# when it corresponds to CoNS(r) with (resp. without) the use of MVMO property.
# Each plot line in these subfigures corresponds to a specific value of l0 and
# is represented with a specific color. Each plot line includes a shaded colored
# region to depict a range of execution times based on the corresponding initial
# signed network and its perturbed versions.
#
#
# graph.sizes: number of vertices
# d: graph density
# l0.values: number of modules
# force: whether or not the existing files are overwritten by a fresh call of all corresponding methods (e.g partitioning method)
#
############################################################################
plot.benchmark.results = function(graph.sizes, d, l0.values, nb.edit, force)
{
	
	csv.path = file.path(OUTPUT.CSV.FOLDER,paste0("benchmark-results_d=",d,".csv"))
	df = read.csv(csv.path, header=1, check.names=F, row.names=1)
	
	
	# source: https://stackoverflow.com/questions/20924705/plot-negative-values-in-logarithmic-scale-with-ggplot-2
	
	library(ggplot2)
	library(scales)
	library(ggallin)
	
	
	csv.n.values.str = sapply(row.names(df), function(x) unlist(strsplit(x,","))[1])
	csv.n.values = as.integer(sapply(csv.n.values.str, function(x) unlist(strsplit(x,"="))[2]))
	csv.l0.values.str = sapply(row.names(df), function(x) unlist(strsplit(x,","))[2])
	csv.l0.values = as.numeric(sapply(csv.l0.values.str, function(x) unlist(strsplit(x,"="))[2]))
	csv.prop.mispl.values.str = sapply(row.names(df), function(x) unlist(strsplit(x,","))[4])
	csv.prop.mispl.values = as.numeric(sapply(csv.prop.mispl.values.str, function(x) unlist(strsplit(x,"="))[2]))
	
	
	
	desc1 = paste0(CoNSCC," ",nb.edit,"edit")
	desc2 = paste0(CoNSCC," ",nb.edit,"edit without MVMO pruning")
	df.out = c()
	for(n in graph.sizes){
	    for(l0 in l0.values){
	        indxs = which(csv.n.values == n & csv.l0.values == l0)
	        subdf = df[indxs,]
	        without = subdf[,desc2]
	        min.without = min(without, na.rm=TRUE)
	        max.without = max(without, na.rm=TRUE)
	        avg.without = (min.without+max.without)/2
	        with = subdf[,desc1]
	        min.with = min(with, na.rm=TRUE)
	        max.with = max(with, na.rm=TRUE)
	        avg.with = (min.with+max.with)/2
	        df.out = rbind(df.out, c(as.character(n),as.character(l0),min.with,max.with,avg.with,min.without,max.without,avg.without))
	    }
	}
	colnames(df.out) = c("n","l0","min.with","max.with","avg.with","min.without","max.without","avg.without")
	data = data.frame(df.out)
	data[,"n"] = as.integer(data[,"n"])
	data[,"min.with"] = as.numeric(data[,"min.with"])
	data[,"max.with"] = as.numeric(data[,"max.with"])
	data[,"avg.with"] = as.numeric(data[,"avg.with"])
	data[,"min.without"] = as.numeric(data[,"min.without"])
	data[,"max.without"] = as.numeric(data[,"max.without"])
	data[,"avg.without"] = as.numeric(data[,"avg.without"])
	
	print(data)
	
	ggplot(data, aes(x = n))+
	     # ylim(0, 8000)+
	     geom_line(aes(y = avg.with, color = l0), size=1.5)+
	     geom_line(aes(y = avg.without, color = l0), size=1.5, linetype = "dashed")+
	     #scale_color_manual(values = c("darkred", "steelblue", "#E3CF57"))+
	     scale_y_continuous(trans = pseudolog10_trans,breaks=c(-5000,-1000,-100,-10,-1,0,1,10,100,1000,5000), limits=c(0,10000))+
	     geom_ribbon(data=data,aes(ymin=min.with,ymax=max.with,group=l0, fill=l0),alpha=0.2)+
	     geom_ribbon(data=data,aes(ymin=min.without,ymax=max.without,group=l0, fill=l0),alpha=0.2)+
	     theme_light(base_size = 21)
 
 	ggsave(filename=file.path(OUTPUT.CSV.FOLDER,paste0("benchmark_d=",d,"_nbEdit",nb.edit,".pdf")),
		 width=16)
	
	
	#sp = sp + geom_hline(yintercept=0, linetype="dashed", color="black")
	#sp = sp + geom_vline(xintercept=20, linetype="dashed", color="black")
	#sp = sp + geom_vline(xintercept=40, linetype="dashed", color="black")
	#sp = sp + geom_vline(xintercept=60, linetype="dashed", color="black")

}


