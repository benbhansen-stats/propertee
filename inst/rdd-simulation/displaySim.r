

###########################################################################################
### Table 4 simulation (polynomial)
###########################################################################################


#' @export
## dgms <- function(tp){

##     ### the DGMs
##     curve(0.5*x,from=-1,to=1,ylim=c(-2,2),xlab='r',ylab='$\\EE[Y|R=r]$',main='Linear')
##     abline(v=0,lty=2)
##     curve(ifelse(abs(x)>0.5, 3*x+sign(x)*(0.5-3)*0.5,0.5*x),from=-1,to=1,ylim=c(-2,2),xlab='r',ylab='$\\EE[Y|R=r]$',main='Anti-Symmetric')
##     abline(v=0,lty=2)
##   #curve(ifelse(x>0.5,3*x+(0.5-3)*0.5,0.5*x),from=-1,to=1,ylim=c(-2,2),xlab='r',ylab='$\\EE[Y|R=r]$',main='One-Sided')
##   curve(mu4,from=-1,to=1,ylim=c(-2,2),xlab='r',ylab='$\\EE[Y|R=r]$',main='One-Sided')
##     abline(v=0,lty=2)
## }

dgms <- function(){
  lin <- function(x) 0.5*x
  as <- function(x) ifelse(abs(x)>0.5,3*x+sign(x)*(0.5-3)*0.5,lin(x))
  mu4 <- function(x) sin(3*x)
  p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))+xlim(-1,1)+
    theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      #axis.title.x='$r$',
      #axis.title.y='$\\EE[Y|R=r]$',
      legend.position="none",
      ## panel.background=element_blank(),
      ## panel.border=element_blank(),
      ## panel.grid.major=element_blank(),
      ## panel.grid.minor=element_blank(),
      ## plot.background=element_blank()
    )+xlab('$r$')+ylab('$\\EE[Y|R=r]$')

  gridExtra::grid.arrange(p+stat_function(fun=lin)+ggtitle('Linear'),
    p+stat_function(fun=as)+ggtitle('Anti-Symmetric'),
    p+stat_function(fun=mu4)+ggtitle('Sine'),nrow=1)
}

resTab <- function(run,full=FALSE){
    llr <- run[c(nrow(run)-1,nrow(run)),]
    run <- run[-c(nrow(run)-1,nrow(run)),]

    SEs <- run[seq(1,nrow(run)-1,2),]
    ests <- run[seq(2,nrow(run),2),]

    runSEflex <- SEs[seq(1,nrow(SEs),2),]
    runSEols <- SEs[seq(2,nrow(SEs),2),]

    runEstflex <- ests[seq(1,nrow(ests),2),]
    runEstols <- ests[seq(2,nrow(ests),2),]


    tabFun <- function(runSE,runEst,full){
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
    list(tabFLEX=tabFun(runSEflex,runEstflex,full=full),
         tabOLS=tabFun(runSEols,runEstols,full=full),
         tabLLR=tabFun(rbind(llr[1,]),rbind(llr[2,]),full=full))
}


## resTab <- function(run,full=FALSE){
##   run <- run[,apply(run,2,function(x) !any(is.na(x)))]
##   tabFun <- function(ests)
##       rbind(
##         bias=rowMeans(rbind(ests)),
##         RMSE=sqrt(rowMeans(rbind(ests^2))))

##   tabFunP <- function(ps)
##     apply(rbind(ps),1,function(x) mean(x<0.05,na.rm=TRUE))


##   out <- sapply(c('sh','ik','ll'),function(mm) tabFun(run[grep(paste0(mm,'.est'),rownames(run)),]),
##     simplify=FALSE)

##   if(full) out <- sapply(c('sh','ik','ll'),
##     function(mm) rbind(out[[mm]],level= tabFunP(run[grep(paste0(mm,'.p'),rownames(run)),])),
##     simplify=FALSE)

##   out
## }

#' @export
prntTab <- function(totalPoly,maxDeg=4,full=TRUE,md=FALSE){
    tab <- NULL

    if(maxDeg>(nrow(totalPoly[[1]])-2)/4){
        maxDeg <- (nrow(totalPoly[[1]])-2)/4
        if(!full)
            warning(paste('There don\'t seem to be that many degrees in the simulation.\n Setting maxDeg=',maxDeg))
    }

    ctab <- function(runname,full){
        run <- totalPoly[[runname]]
        res <- resTab(run,full=full)

        cbind(res[['tabFLEX']][,1:maxDeg],
              res[['tabOLS']][,1:maxDeg],
              res[['tabLLR']])
    }
    for(dgm in c('lin','antiSym','wass')){
        tab <- rbind(tab,ctab(paste0("reg=",dgm,'_t'),full))
        if(full) if(paste0(dgm,'_norm')%in%names(totalPoly)) tab <- rbind(tab,ctab(paste0(dgm,'_norm'),full))
    }
    if(md){
        colnames(tab) <- c(paste('LRD, deg=',1:maxDeg),paste('OLS, deg=',1:maxDeg),'Loc.Lin')
        rownames(tab) <- paste(rep(c('lin','antiSym','sine'),
          each=nrow(tab)/sum(rownames(tab)=='bias')),#,times=2),
          #rep(c('t err','norm err'),each=nrow(tab)/2),
          rownames(tab))
    }
    return(tab)
}

#' @export
polyLatex <- function(tab,full,caption='',label='tab:poly',braks='ht'){
    if(NCOL(tab)!=9) stop('This only works with polynomial degree=1,...,4')
    cat('
        \\begin{table}[ht]
\\centering
\\begin{tabular}{cr|llll|llll|l',ifelse(full,'|llll|llll|l}','}'),'
  \\hline \n')
    if(full) cat('&&\\multicolumn{9}{c|}{$t_3$ Errors} &\\multicolumn{9}{c|}{$\\mathcal{N}(0,1)$ Errors} \\\\ \n')

 cat('&& \\multicolumn{4}{c|}{Limitless} &  \\multicolumn{4}{c|}{OLS} &\\makecell[c]{Local\\\\Linear}',
        ifelse(full,'\\multicolumn{4}{c|}{Limitless} &  \\multicolumn{4}{c|}{OLS} &\\makecell[c]{Local\\\\Linear}',''),'\\\\
 \\multicolumn{2}{r|}{\\makecell[r]{Polynomial\\\\Degree}}&1&2&3&4&1&2&3&4&',ifelse(full,'&1&2&3&4&1&2&3&4&n/a','n/a'),' \\\\
')
    for(rr in 1:nrow(tab)){
        if(rr==1) cat('\\hline\n\\hline\n\\multirow{',ifelse(full,4,2),'}{*}{',ifelse(full,'\\begin{sideways}Linear\\end{sideways}','Linear'),'}')
        if(rr==3) cat('\\hline\n\\hline\n\\multirow{',ifelse(full,4,2),'}{*}{',ifelse(full,'\\begin{sideways}Anti-Sym\\end{sideways}','Anti-Sym'),'}')
        if(rr==5) cat('\\hline\n\\hline\n\\multirow{',ifelse(full,4,2),'}{*}{',ifelse(full,'\\begin{sideways}One-Side\\end{sideways}','One-Side'),'}')
        cat('&',rownames(tab)[rr],'&')
        cat(paste(sprintf("%.1f", round(tab[rr,],1)),collapse='&'))
        cat('\\\\ \n')
    }
    cat('
 \\hline
\\end{tabular}
\\caption{',caption,'}
\\label{',label,'}
\\end{table}\n',sep='')

}

#' @export
polyLatex5 <- function(tab,full,caption='',label='tab:poly'){
  if(NCOL(tab)!=11) stop('This only works with polynomial degree=1,...,5')
  tab2 <- tab
  for(i in 1:nrow(tab)) for(j in 1:ncol(tab))
                          tab2[i,j] <- ifelse(tab[i,j]>10,
                            sprintf("%i",as.integer(round(tab[i,j]))),
                            sprintf("%.1f",round(tab[i,j],1)))

    cat('
        \\begin{table}[ht]
\\centering
\\begin{tabular}{ll|ccccc|ccccc|c',ifelse(full,'|llll|llll|l}','}'),'
  \\hline \n')
    if(full) cat('&&\\multicolumn{11}{c|}{$t_3$ Errors} &\\multicolumn{11}{c|}{$\\mathcal{N}(0,1)$ Errors} \\\\ \n')

 cat('&& \\multicolumn{5}{c|}{Limitless} &  \\multicolumn{5}{c|}{OLS} &Local',
        ifelse(full,'\\multicolumn{5}{c|}{Limitless} &  \\multicolumn{5}{c|}{OLS} &Local',''),'\\\\
 && \\multicolumn{5}{c|}{Polynomial Degree}&\\multicolumn{5}{c|}{Polynomial Degree}&Linear',
 ifelse(full,'\\multicolumn{5}{c|}{Polynomial Degree}&\\multicolumn{5}{c|}{Polynomial Degree}&Linear',''),'\\\\
 DGM&Measure&1&2&3&4&5&1&2&3&4&5&',ifelse(full,'&1&2&3&4&1&2&3&4&5&',''),' \\\\
')
    for(rr in 1:nrow(tab)){
        if(rr==1) cat('\\hline\n\\hline\n\\multirow{',ifelse(full,4,2),'}{*}{',ifelse(full,'\\begin{sideways}Linear\\end{sideways}','Linear'),'}',sep='')
        if(rr==3) cat('\\hline\n\\hline\n\\multirow{',ifelse(full,4,2),'}{*}{',ifelse(full,'\\begin{sideways}Anti-Sym\\end{sideways}','\\makecell[c]{Anti-\\\\Sym}'),'}',sep='')
        if(rr==5) cat('\\hline\n\\hline\n\\multirow{',ifelse(full,4,2),'}{*}{',ifelse(full,'\\begin{sideways}One-Side\\end{sideways}','Sine'),'}',sep='')
        cat('&',rownames(tab)[rr],'&')
        cat(paste(tab2[rr,],collapse='&'))
        cat('\\\\ \n')
    }
    cat('
 \\hline
\\end{tabular}
\\caption{',caption,'}
\\label{',label,'}
\\end{table}\n',sep='')

}


simpleSummaryOneEstimator=function(est,eff='TOT'){
  err <- est['est',]-est[eff,]
  V <- var(err,na.rm=TRUE)
  c(
    `bias/V`=mean(err,na.rm=TRUE)/V,
    RMSE=sqrt(mean(err^2,na.rm=TRUE)),
    samplingVarRelativeBias=(mean(est['se',]^2,na.rm=TRUE)-V)/V,
    ciCoverage=mean(est['ciL',]<est[eff,] & est['ciH',]>est[eff,],na.rm=TRUE) 
  )
}

simpleSummary <- function(run,eff='TOT'){
  run=run[apply(run,1,function(x) mean(is.na(x)))<1,]

  effs=run[c('ATE','TOT'),]
  run=run[-which(rownames(run)%in%c('ATE','TOT')),]
  rownamesplit <- do.call('rbind',strsplit(rownames(run),'\\.'))
  
  summ <- NULL
  for(estimator in unique(rownamesplit[,1])){
    est <- run[rownamesplit[,1]==estimator,]
    rownames(est) <- rownamesplit[rownamesplit[,1]==estimator,2]
    est <- rbind(est,effs)
    summ1=simpleSummaryOneEstimator(est,eff=eff)
    names(summ1) <- paste(estimator,names(summ1),sep='.')
    summ <- c(summ,summ1)
  }
  summ
}
