tryNA <- function(expr,num=1){
    out <- try(expr,silent=TRUE)
    if(inherits(out,'try-error')) return(rep(NA,num))
    out
}


###############
### table 4 simulation (polynomial)
###############
mu3 <- function(x){
    ifelse(x<0,1.27*x-.5*7.18*x^2+0.7*20.21*x^3+1.1*21.54*x^4+1.5*7.33*x^5,
    .84*x-0.1*3.00*x^2-0.3*7.99*x^3-0.1*9.01*x^4+3.56*x^5)
}

### from Wasserman "all of nonparametric stats" p.99 (just the sine function really)
mu4 <- function(x)
  sin(3*x)

makeDataShapes <- function(n,shape,tdist=FALSE,tau=0,plt=FALSE){
    curve <- 3
    R <- runif(n,-1,1)

    yc <- .5*R
    if(shape=='antiSym')
        yc <- ifelse(abs(R)>0.5,3*R+sign(R)*(0.5-3)*0.5,yc)
    if(shape=='oneSide')
        yc <- yc+ ifelse(R> 0.5,3*R+(0.5-3)*0.5,yc)
    if(shape=='sym')
        yc <- ifelse(abs(R)> 0.5,sign(R)*3*R+(sign(R)*0.5-3)*0.5,yc)
    if(shape=='cct')
        yc <- mu3(R)
    if(shape=='poly3')
      yc <- 1.75*R^3
    if(shape=='wass')
      yc <- mu4(R)

    if(plt) plot(R,yc)

    yc <- yc+if(tdist) rt(n,3) else rnorm(n)

    Z <- R>0

    Y <- yc+Z*tau

    id <- 1:n

    data.frame(R=R,Z=Z,Y=Y,id=id)
}



olsPoly <- function(dat,deg){
    mod <- lm(Y~poly(R,deg,raw=TRUE)*Z,data=dat)
    Z.pos <- which(grepl('Z',names(coef(mod))) & !grepl(':',names(coef(mod))))
    setNames(c(se=summary(mod)$coef[Z.pos,'Std. Error'],est=coef(mod)[Z.pos]),
             paste0(c('ols.se.','ols.est.'),deg))
}

flexPoly <- function(dat,deg){
  des <- rd_design(Z ~ forcing(R) + unitid(id), data=dat)

### we can use poly(.) once issue #61 is resolved
### covForm <- Y~Z+poly(R,deg,raw=TRUE)
### until then
  covForm <- "Y~R+Z"
  if(deg>1)
    covForm <- paste0(covForm,'+',paste(paste0("I(R^",2:deg,")"),collapse='+'))
  covForm <-as.formula(covForm)

  covMod <- lmrob(covForm,data=dat,method='MM')

  res <- lmitt(Y~assigned(),design=des,offset=cov_adj(covMod),data=dat)|>
    summary()|>
    getElement('coefficients')


  setNames(c(res[2,'Std. Error'],res[2,'Estimate']),
            paste0(c('flex.se.','flex.est.'),deg))
}

llPoly <- function(dat){
    mod <- RDestimate(Y~R,data=dat)
    setNames(c(mod$p[1],mod$est[1]),c('ll.p','ll.est'))
}


polySim <- function(n,degs=1:5,shape='lin',tdist=TRUE,tau=0){
    dat <- makeDataShapes( n=n,shape=shape,tdist=tdist,tau=tau)

    func <- function(deg){
        c(tryNA(flexPoly(dat,deg),2),tryNA(olsPoly(dat,deg),2))
    }

    c(do.call('c',lapply(degs,func)),tryNA(llPoly(dat),2))
}


polyDisp <- function(sim){
    if(nrow(sim)<ncol(sim)) sim <- t(sim)


    simp <- sim[,grep('p',colnames(sim))]
    simEst <- sim[,grep('est',colnames(sim))]
    print('level of test')
    print(apply(simp,2,function(x) mean(x<0.05,na.rm=TRUE)))
    print('mean of est')
    print(apply(simEst,2,mean,na.rm=TRUE))
    print('SD of est')
    print(apply(simEst,2,sd,na.rm=TRUE))
    print('RMSE')
    print(apply(simEst,2,function(x) mean(x^2,na.rm=TRUE)))

}


#' Run the polynomial simulation from Table 4
#'
#' @import robustbase
#'
#' @param nreps Number of simulation replications
#'
#' @return list of output for each simulation run
#' @export
#'
totalPolySim <- function(nreps=5000,cluster=NULL){
  appFunc <- if(is.null(cluster)) sapply else function(X,FUN) parSapply(cl=cluster,X=X,FUN=FUN)
  res <- list()
    #B=5000
    n=500
    tau=0
    degs <- 1:5

    for(shape in c('lin','antiSym','wass')){
        for(tdist in c(TRUE)){#,FALSE)){
            res[[paste0(shape,'_',ifelse(tdist,'t','norm'))]] <-
                appFunc(1:nreps,function(i) polySim(n,degs,shape,tdist,tau))
        }
    }
    res
}
