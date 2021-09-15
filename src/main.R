
# =============================================================================
# VARIABLES
#   ==> REmark 1: Do not forget to set CPLEX.BIN.PATH correclty in src/define-algos.R
#   ==> Remark 2: PROP.NEGS should be set to 'NA' when DENSITY=1.0
#   ==> Remark 3: CORE.PART.THRESHOLD should be set to 1. Because, when it is less than 1,
#                  there might be multiple way of building core part.
#   ==> Remark 4: It is the responsability of the user who will ensure if the RAM requirement
#                  of his/her system is ok for large graphs, because Cplex may require 
#                  lots of RAM for graphs whose size is larger than 28.
# =============================================================================


## libraries for parallel processing
#library(foreach)
#library(doParallel)

source("src/define-imports.R")


#######################################################################
# STEP 0: Set parameter values
#######################################################################=

BENCHMARK.GRAPH.SIZES = c(20,30) #c(30,40,50)
BENCHMARK.DENSITY = 0.25
BENCHMARK.L0.VALUES = c(3,4)
BENCHMARK.PROP.NEGS = c(0.3,0.5) #c(0.3,0.5,0.7)
BENCHMARK.FORCE = FALSE
BENCHMARK.PLOT.FORMAT = JUST.PLOT


plot.format <- c( # ==========> it is not taken into account everywhere !! TODO
		#PLOT.AS.PDF
		#PLOT.AS.JPEG
		#PLOT.AS.PNG
		JUST.PLOT
)

FORCE = FALSE

keep.algo.log.files = TRUE




#######################################################################
# STEP 1: Generate perfectly balanced signed networks and place them into 'in' folder
#######################################################################

#generate.all.perfectly.balanced.networks(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.PROP.NEGS)




#######################################################################
# STEP 2: Perturb in many ways the perfectly balanced signed networks
#			and obtain new signed networks
#######################################################################

create.all.extreme.imbalanced.networks(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.PROP.NEGS, BENCHMARK.FORCE, BENCHMARK.PLOT.FORMAT)
perform.all.benchmark(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.FORCE)
collect.all.benchmark.results(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, BENCHMARK.FORCE)

plot.benchmark.results(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, 3, BENCHMARK.FORCE)
plot.benchmark.results(BENCHMARK.GRAPH.SIZES, BENCHMARK.DENSITY, BENCHMARK.L0.VALUES, 4, BENCHMARK.FORCE)


