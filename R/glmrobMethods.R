##' @title Extract empirical estimating equations from a \code{glmbrob} model
##'   fit
##' @param x a fitted \code{glmrob} object
##' @param ... arguments passed to methods
##' @return matrix, estimating functions evaluated at data points and fitted
##'   parameters
##' @importFrom stats family model.matrix
##' @importFrom utils getFromNamespace
##' @rdname glmrob_methods
##' @exportS3Method
estfun.glmrob <- function(x, ...) {

  ## numbers refer to which lines in
  ## https://github.com/cran/robustbase/blob/335b69f2310bd21ca4cdfc17a2a99ebbcad84017/R/glmrobMqle.R
  ## the expressions were derived or copied from

  ## this was copied verbatim from lines 121-155
  switch(family(x)$family,
         "binomial" = {
           Epsi.init <- utils::getFromNamespace("EpsiBin.init", "robustbase")
           Epsi <- utils::getFromNamespace("EpsiBin", "robustbase")
           EpsiS <- utils::getFromNamespace("EpsiSBin", "robustbase")
           Epsi2 <- utils::getFromNamespace("Epsi2Bin", "robustbase")
           phiEst <- phiEst.cl <- 1
         },
         "poisson" = {
           Epsi.init <- utils::getFromNamespace("EpsiPois.init", "robustbase")
           Epsi <- utils::getFromNamespace("EpsiPois", "robustbase")
           EpsiS <- utils::getFromNamespace("EpsiSPois", "robustbase")
           Epsi2 <- utils::getFromNamespace("Epsi2Pois", "robustbase")
           phiEst <- phiEst.cl <- expression({1})
         },
         "gaussian" = {
           Epsi.init <- utils::getFromNamespace("EpsiGaussian.init", "robustbase")
           Epsi <- utils::getFromNamespace("EpsiGaussian", "robustbase")
           EpsiS <- utils::getFromNamespace("EpsiSGaussian", "robustbase")
           Epsi2 <- utils::getFromNamespace("Epsi2Gaussian", "robustbase")
           phiEst.cl <- utils::getFromNamespace("phiGaussianEst.cl", "robustbase")
           phiEst <- utils::getFromNamespace("phiGaussianEst", "robustbase")
         },
         "Gamma" = { ## added by ARu
           # Epsi.init <- utils::getFromNamespace("EpsiGamma.init", "robustbase")
           #### The `Gmn` function inside `EpsiGamma.init` was not being
           #### properly found. Bringing it in here manually and using
           #### `getFromNamespace` to locate it. JE 3/2024.
           Epsi.init <- expression({
             nu <- 1/phi
             snu <- 1/sqrt(phi)
             pPtc <- pgamma(snu + c(-tcc, tcc), shape = nu, rate = snu)
             pMtc <- pPtc[1]
             pPtc <- pPtc[2]
             aux2 <- tcc * snu
             GLtcc <- getFromNamespace("Gmn", "robustbase")(-tcc, nu)
             GUtcc <- getFromNamespace("Gmn", "robustbase")(tcc, nu)
           })
           Epsi <- utils::getFromNamespace("EpsiGamma", "robustbase")
           EpsiS <- utils::getFromNamespace("EpsiSGamma", "robustbase")
           Epsi2 <- utils::getFromNamespace("Epsi2Gamma", "robustbase")
           phiEst.cl <- utils::getFromNamespace("phiGammaEst.cl", "robustbase")
           phiEst <- utils::getFromNamespace("phiGammaEst", "robustbase")
           },
         ## else
         stop(gettextf("family '%s' not yet implemented", family$family),
              domain=NA)
         )

  ### change names of returned objects to what they were in original robustbase
  ### code
  ni <- x$ni
  tcc <- x$tcc
  eta <- x$linear.predictors
  phi <- x$dispersion
  mu <- x$fitted.values

  ### get V(mu)
  rr <- x$y-mu
  sni <- sqrt(ni)
  sVF <- sni*rr/x$residuals ## 164 (inverted)
  sV <- sVF*sqrt(x$dispersion) #168
  Vmu <- sVF^2 ## 163 (inverted)


 ### robustbase code for calculating cpsi=psi_c - E[psi_c]
  ##177-78
  K <- floor(mu*ni + tcc* sni*sV)
  H <- floor(mu*ni - tcc* sni*sV)

  eval(Epsi.init)
  residPS <- x$residual/sqrt(phi) # scaled Pearson residuals (169)
  cpsi <- pmax.int(-tcc,pmin.int(residPS,tcc))-eval(Epsi) ## 210

  ### dPsi/dBeta=dPsi/d eta d eta/d Beta=dPsi/d eta X
  dmu.deta <- family(x)$mu.eta(eta)
  X <- model.matrix(x)

  cpsi*x$w.x*sni/sV*dmu.deta*X
}

##' @title Extract bread matrix from an \code{lmrob()} fit
##' @param x a fitted \code{lmrob} object
##' @param ... arguments passed to methods
##' @return matrix, inverse Hessian of loss as evaluated at fitted parameters
##' @rdname glmrob_methods
##' @exportS3Method
bread.glmrob <- function(x, ...) {
  return(solve(x$matM))
}
