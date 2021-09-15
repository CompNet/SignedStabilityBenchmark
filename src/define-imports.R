#############################################################################################
# Just loads the scripts necessary to process the dataset. This script is needed when
# using foreach (parallel processing), since each worker must load all these dependencies.
# 
#############################################################################################

source("src/define-constants.R")
source("src/define-algos.R")
source("src/common.R")
source("src/define-random-generator.R")

source("src/define-paths.R")
source("src/benchmark-analysis/create-extreme-imbalanced-networks.R")
source("src/benchmark-analysis/perform-benchmark.R")
source("src/benchmark-analysis/collect-benchmark-results.R")
source("src/benchmark-analysis/plot-benchmark.R")


#source("src/filter-lp-miyauchi/filter-lp-miyauchi.R")

source("src/generate-input-networks/generate-perfectly-balanced-networks.R")


source("src/partition-networks/partition-networks.R")
source("src/partition-networks/load-membership.R")

#source("src/generate-input-networks/assess.R")
#source("src/generate-input-networks/common.R")
#source("src/generate-input-networks/generate.R")
#source("src/generate-input-networks/plot.R")

