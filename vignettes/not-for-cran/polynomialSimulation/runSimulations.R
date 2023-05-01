source('simulationFunctions.r')
devtools::load_all()



library('robustbase')
#library('rdd')
#library('RItools')
library('sandwich')
#library('nnet')


if(clust){
  clusterEvalQ(cl,{

    library('robustbase')
    #library('rdd')
    #library('RItools')
    library('sandwich')
    #library('nnet')
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
