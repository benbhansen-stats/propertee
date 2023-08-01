estfun.glmrob <- function(mod,...){

  switch(family(mod)$family,
         "binomial" = {
           Epsi.init <- robustbase:::EpsiBin.init
           Epsi <- robustbase:::EpsiBin
           EpsiS <- robustbase:::EpsiSBin
           Epsi2 <- robustbase:::Epsi2Bin
           phiEst <- phiEst.cl <- 1
         },
         "poisson" = {
           Epsi.init <- robustbase:::EpsiPois.init
           Epsi <- robustbase:::EpsiPois
           EpsiS <- robustbase:::EpsiSPois
           Epsi2 <- robustbase:::Epsi2Pois
           phiEst <- phiEst.cl <- expression({1})
         },
         "gaussian" = {
           Epsi.init <- robustbase:::EpsiGaussian.init
           Epsi <- robustbase:::EpsiGaussian
           EpsiS <- robustbase:::EpsiSGaussian
           Epsi2 <- robustbase:::Epsi2Gaussian
           phiEst.cl <- robustbase:::phiGaussianEst.cl
           phiEst <- robustbase:::phiGaussianEst
         },
         "Gamma" = { ## added by ARu
           Epsi.init <- robustbase:::EpsiGamma.init
           Epsi <- robustbase:::EpsiGamma
           EpsiS <- robustbase:::EpsiSGamma
           Epsi2 <- robustbase:::Epsi2Gamma
           phiEst.cl <- robustbase:::phiGammaEst.cl
           phiEst <- robustbase:::phiGammaEst
           },
         ## else
         stop(gettextf("family '%s' not yet implemented", family$family),
              domain=NA)
         )

  ### change names of returned objects to what they were in original robustbase code
  ni <- mod$ni
  tcc <- mod$tcc
  eta <- mod$linear.predictors
  phi <- mod$dispersion
  mu <- mod$fitted.values

  ### get V(mu)
  rr <- mod$y-mu
  sni <- sqrt(ni)
  sVF <- sni*rr/mod$residuals ## 164 (inverted)
  sV <- sVF*sqrt(mod$dispersion) #168
  Vmu <- sVF^2 ## 163 (inverted)


  ### robustbase code for calculating cpsi=psi_c - E[psi_c]
  ##177-78
  K <- floor(mu*ni + tcc* sni*sV)
  H <- floor(mu*ni - tcc* sni*sV)

  eval(Epsi.init)
  residPS <- mod$residual/sqrt(phi) # scaled Pearson residuals
  cpsi <- pmax.int(-tcc,pmin.int(residPS,tcc))-eval(Epsi) ## 210

  ### dPsi/dBeta=dPsi/d eta d eta/d Beta=dPsi/d eta X
  dmu.deta <- family(mod)$mu.eta(mod$linear.predictors)
  X <- model.matrix(mod)

  cpsi*w.x*sni/sV*dmu.deta*X
}

bread.glmrob <- function(mod,...)
  return(solve(mod$matM))
