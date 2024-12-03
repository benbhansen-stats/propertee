resTab <- function(run,full=FALSE) {
  tabFun <- function(est_type, full) {
    if (est_type == "ll") se_search <- "\\.p" else se_search <- "\\.se"

    runEst <- run[grepl(est_type, rownames(run)) & grepl("est", rownames(run)),]
    runSE <- run[grepl(est_type, rownames(run)) & grepl(se_search, rownames(run)),]
    trueVar <- apply(runEst,1,var,na.rm=TRUE)
    
    if(full){
      ciL=runEst-1.96*runSE
      ciH=runEst+1.96*runSE
      cover=sign(ciL*ciH)<0
    }
    
    tab <- rbind(
      
      VarRelBias=
        runSE^2|>
        sweep(1,trueVar)|>
        rowMeans(na.rm=TRUE)|>
        (\(x) x/trueVar)(),
      RMSE=apply(runEst,1,function(x) sqrt(mean(x^2,na.rm=TRUE))))
    
    if(full) tab <- rbind(tab,
                          bias=rowMeans(runEst,na.rm=TRUE),
                          sd=sqrt(trueVar),
                          ciCover=rowMeans(cover,na.rm=TRUE))
    tab
  }
  list(tabFLEX=tabFun("flex", full=full)
       , tabOLS=tabFun("ols", full=full)
       # , tabLLR=tabFun("ll", full=full)
       )
}

prntTab <- function(totalPoly,maxDeg=4,full=TRUE) {
  tab <- NULL
  
  # ncompare_rows <- nrow(totalPoly[[1]])-2 # this includes LLR
  ncompare_rows <- nrow(totalPoly[[1]])
  if(maxDeg>ncompare_rows/4){
    maxDeg <- ncompare_rows/4
    if(!full)
      warning(
        paste('There don\'t seem to be that many degrees in the simulation.',
              "Setting maxDeg=",maxDeg)
      )
  }
  
  ctab <- function(runname,full){
    run <- totalPoly[[runname]]
    res <- resTab(run,full=full)
    
    cbind(res[['tabFLEX']][,1:maxDeg]
          , res[['tabOLS']][,1:maxDeg]
          # , res[['tabLLR']]
          )
  }
  
  for (dgm in c('lin','antiSym','wass')) {
    tab <- rbind(tab,ctab(paste0(dgm,'_t'),full))
    if(full) {
      if(paste0(dgm,'_norm')%in%names(totalPoly)) {
        tab <- rbind(tab,ctab(paste0(dgm,'_norm'),full))
      }
    }
  }
  
  return(tab)
}
