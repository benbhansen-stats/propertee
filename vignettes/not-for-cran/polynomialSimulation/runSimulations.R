source('simulationFunctions.r')
devtools::load_all()

nreps <- commandArgs(TRUE)
if (length(nreps) > 1) {
  stop("Only expecting one variable corresponding to `NREPS`")
}

# Initialization. Note that `nreps=0` corresponds to no simulations,
# just print results from previously saved simulations.
# In order to re-run the simulations, the `nreps`
# variable should have been set to a positive integer before initiating this script.
# 
# To run the simulations in parallel, using the `parallel` package in
# `R`,
# register a cluster, called `cl` with the desired number of nodes, with
# code similar to the following:
# 
if (length(nreps) == 0) nreps <- 0 else nreps <- as.numeric(nreps)
nreps
if (nreps) {
  library(parallel)
  cl <- makeCluster(5)
  
  library('robustbase')
  # library('rdd')
  # library('RItools')
  library('sandwich')
  # library('nnet')

  clust <- FALSE
  if(require('parallel')) if(exists('cl')) if(inherits(cl,"cluster")) clust <- TRUE

  if(clust){
    clusterEvalQ(cl,{

      library('robustbase')
      # library('rdd')
      # library('RItools')
      library('sandwich')
      # library('nnet')
      devtools::load_all("../..")
      source('simulationFunctions.r')
    })
  } else cl <- NULL
  
  ## Run the simulation
  set.seed(201609)
  st2 <- system.time(totalPoly <- totalPolySim(nreps,
                                               n = 500,
                                               degs = seq_len(5),
                                               tau = 0,
                                               cluster = cl))
  save(totalPoly,file="./totalPolySim.RData")
  cat(paste0(date(), ', nreps=', nreps, '\n'),
      paste(c(names(st2),'\n', collapse=T)),
      st2,
      file='totalPolySim-runtime.txt', append=TRUE)
}
