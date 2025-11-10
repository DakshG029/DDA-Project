set.seed(101)
Gamma.VaR.CVaR=function(alphaX=2, betaX=2, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qgamma(pu,alphaX,betaX)
  X <- rgamma(n,alphaX,betaX)
  print(head(X))
  Xu <- X[X>u]-u
  MXu=max(Xu)
  m <- length(Xu)
  
  # Real value
  VaRT <- qgamma(pk,alphaX,betaX)
  CVaRT <- 1/((1-pk)*gamma(alphaX)*betaX)*gamma_inc(alphaX+1,betaX*VaRT)
  
  ## Empiric
  VaRE <- quantile(X,pk)[[1]]
  CVaRE <- mean(X[X>=VaRE])
  
  ## Parametric
  library(gsl)
  GammaNewton <- function(x,alpha0)
  {
    m <- mean(x)
    mlog <- mean(log(x))
    alpha0*(1+(digamma(alpha0)+log(m/alpha0)-mlog)/(1-alpha0*trigamma(alpha0)))
  }
  emvalphaGamma <- function(x)
  {
    m <- mean(x)
    v <- var(x)
    alpha0 <- m^2/v
    error <- 1
    i <- 0
    while(error>1E-4)
    {
      i <- i+1
      alpha1 <- GammaNewton(x,alpha0)
      error <- abs(alpha1-alpha0); alpha0 <- alpha1
    }
    alpha1
  }
  alphaX1 <- emvalphaGamma(X); betaX1 <- alphaX1/mean(X)
  VaRP <- qgamma(pk,alphaX1,betaX1)
  CVaRP <- 1/((1-pk)*gamma(alphaX1)*betaX1)*gamma_inc(alphaX1+1,betaX1*VaRP)
  
  ## Bootstrap
  MR <- function(x){
    rX <- sample(X,n,rep=T)
    VaRB <- quantile(rX,pk)[[1]]
    CVaRB <- mean(rX[rX>=VaRB])
    c(VaRB,CVaRB)
  }
  EstR <- apply(replicate(B,MR(sample(X,n,rep=T))),1,mean)
  VaRB <- 2*VaRE - EstR[1]
  CVaRB <- 2*CVaRE - EstR[2]
  
  ## MCMC
  pcola <- 1-(1-pk)/(1-pu)
  nburn <- 3000
  nthin <- 50
  ndraw <- 10000
  
  xiMHi <- 0.1; sd_xiMH <- 0.01
  sigmaMHi <- ifelse(m>=2,var(Xu)/mean(Xu),1/betaX); sd_sigmaMH <- 1
  alphaBDi <- mean(X)^2/var(X); sd_alphaBD <- 0.5
  betaBDi <- mean(X)/var(X); sd_betaBD <- 0.5
  xiIBDi <- -0.032+0.014/alphaBDi
  sigmaIBDi <- 1/betaBDi*(0.5+0.5*sqrt(alphaBDi))
  for(i in 1:nburn){
    #MH
    xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
    sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
    mXu <- min(1+xiMHcan*Xu/sigmaMHcan)
    while(mXu<0||sigmaMHcan<0||xiMHcan<(-0.5))
    {
      xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
      sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
      mXu <- min(1+xiMHcan*Xu/sigmaMHcan)
    }
    pMH <- -(m+1)*log(sigmaMHcan) -
      (1+xiMHcan)/xiMHcan*sum(log(1+xiMHcan*Xu/sigmaMHcan)) - 
      log(1+xiMHcan) - 0.5*log(1+2*xiMHcan) + 
      (m+1)*log(sigmaMHi) + 
      (1+xiMHi)/xiMHi*sum(log(1+xiMHi*Xu/sigmaMHi)) + 
      log(1+xiMHi) + 0.5*log(1+2*xiMHi)
    if(log(runif(1))<pMH)
    {
      xiMHi <- xiMHcan
      sigmaMHi <- sigmaMHcan
    }
    VaRMH <- sigmaMHi/xiMHi*((1-pcola)^(-xiMHi)-1) + u
    CVaRMH <- sigmaMHi/xiMHi*(((1-pcola)^(-xiMHi))/(1-xiMHi)-1) + u
    
    #BD
    alphaBDcan <- rnorm(1,alphaBDi,sd_alphaBD)
    betaBDcan <- rnorm(1,betaBDi,sd_betaBD)
    while(alphaBDcan<0||betaBDcan<0)
    {
      alphaBDcan <- rnorm(1,alphaBDi,sd_alphaBD)
      betaBDcan <- rnorm(1,betaBDi,sd_betaBD)
    }
    pBD <- -n*log(gamma(alphaBDcan)) + 
      (n*alphaBDcan-1)*log(betaBDcan) + 
      (alphaBDcan-1)*sum(log(X)) - betaBDcan*sum(X) + 
      0.5*log(alphaBDcan*trigamma(alphaBDcan)-1) +
      n*log(gamma(alphaBDi)) - 
      (n*alphaBDi-1)*log(betaBDi) -
      (alphaBDi-1)*sum(log(X)) + betaBDi*sum(X) - 
      0.5*log(alphaBDi*trigamma(alphaBDi)-1) 
    if(log(runif(1))<pBD)
    {
      alphaBDi <- alphaBDcan 
      betaBDi <- betaBDcan
    }
    VaRBD <- qgamma(pk,alphaBDi,betaBDi)
    CVaRBD <- 1/((1-pk)*gamma(alphaBDi)*betaBDi)*gamma_inc(alphaBDi+1,betaBDi*VaRBD)
    
    #IBD
    xiIBDcan <- -0.032+0.014/alphaBDi
    sigmaIBDcan <- 1/betaBDi * (0.5+0.5*sqrt(alphaBDi))
    pIBD <- -m*log(sigmaIBDcan) -
      (1+1/xiIBDcan)*sum(log(1+xiIBDcan*Xu/sigmaIBDcan)) +
      m*log(sigmaIBDi) + 
      (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
    if(log(runif(1))<pIBD)
    {
      xiIBDi <- xiIBDcan
      sigmaIBDi <- sigmaIBDcan
    }
    VaRIBD <- sigmaIBDi/xiIBDi*((1-pcola)^(-xiIBDi)-1) + u
    CVaRIBD <- sigmaIBDi/xiIBDi*(((1-pcola)^(-xiIBDi))/(1-xiIBDi)-1) + u
  }
  draw <- c(xiMHi, sigmaMHi, VaRMH, CVaRMH,
            alphaBDi, betaBDi, VaRBD, CVaRBD, 
            xiIBDi, sigmaIBDi, VaRIBD, CVaRIBD)
  for(i in 1:ndraw)
  {
    for(j in 1:nthin)
    {
      #MH
      xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
      sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
      mXu <- min(1+xiMHcan*Xu/sigmaMHcan)
      while(mXu<0||sigmaMHcan<0||xiMHcan<(-0.5))
      {
        xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
        sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
        mXu <- min(1+xiMHcan*Xu/sigmaMHcan)
      }
      pMH <- -(m+1)*log(sigmaMHcan) -
        (1+xiMHcan)/xiMHcan*sum(log(1+xiMHcan*Xu/sigmaMHcan)) - 
        log(1+xiMHcan) - 0.5*log(1+2*xiMHcan) + 
        (m+1)*log(sigmaMHi) + 
        (1+xiMHi)/xiMHi*sum(log(1+xiMHi*Xu/sigmaMHi)) + 
        log(1+xiMHi) + 0.5*log(1+2*xiMHi)
      if(log(runif(1))<pMH)
      {
        xiMHi <- xiMHcan
        sigmaMHi <- sigmaMHcan
      }
      VaRMH <- sigmaMHi/xiMHi*((1-pcola)^(-xiMHi)-1) + u
      CVaRMH <- sigmaMHi/xiMHi*(((1-pcola)^(-xiMHi))/(1-xiMHi)-1) + u
      
      #BD
      alphaBDcan <- rnorm(1,alphaBDi,sd_alphaBD)
      betaBDcan <- rnorm(1,betaBDi,sd_betaBD)
      while(alphaBDcan<0||betaBDcan<0)
      {
        alphaBDcan <- rnorm(1,alphaBDi,sd_alphaBD)
        betaBDcan <- rnorm(1,betaBDi,sd_betaBD)
      }
      pBD <- -n*log(gamma(alphaBDcan)) + 
        (n*alphaBDcan-1)*log(betaBDcan) + 
        (alphaBDcan-1)*sum(log(X)) - betaBDcan*sum(X) + 
        0.5*log(alphaBDcan*trigamma(alphaBDcan)-1) +
        n*log(gamma(alphaBDi)) - 
        (n*alphaBDi-1)*log(betaBDi) -
        (alphaBDi-1)*sum(log(X)) + betaBDi*sum(X) - 
        0.5*log(alphaBDi*trigamma(alphaBDi)-1) 
      if(log(runif(1))<pBD)
      {
        alphaBDi <- alphaBDcan
        betaBDi <- betaBDcan
      }
      VaRBD <- qgamma(pk,alphaBDi,betaBDi)
      CVaRBD <- 1/((1-pk)*gamma(alphaBDi)*betaBDi)*gamma_inc(alphaBDi+1,betaBDi*VaRBD)
      
      #IBD
      xiIBDcan <- -0.032+0.014/alphaBDi
      sigmaIBDcan <- 1/betaBDi * (0.5+0.5*sqrt(alphaBDi))
      pIBD <- -m*log(sigmaIBDcan) -
        (1+1/xiIBDcan)*sum(log(1+xiIBDcan*Xu/sigmaIBDcan)) +
        m*log(sigmaIBDi) +
        (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
      if(log(runif(1))<pIBD)
      {
        xiIBDi <- xiIBDcan
        sigmaIBDi <- sigmaIBDcan
      }
      VaRIBD <- sigmaIBDi/xiIBDi*((1-pcola)^(-xiIBDi)-1) + u
      CVaRIBD <- sigmaIBDi/xiIBDi*(((1-pcola)^(-xiIBDi))/(1-xiIBDi)-1) + u
    }
    draw <- c(draw, c(xiMHi, sigmaMHi, VaRMH, CVaRMH,
                      alphaBDi, betaBDi, VaRBD, CVaRBD,
                      xiIBDi, sigmaIBDi, VaRIBD, CVaRIBD))
  }
  draw <- matrix(draw,ncol=12,byrow = TRUE)
  list(alpha=alphaX, beta=betaX, length=n, threshold=pu, p=pk, 
       real.value.VaR=VaRT, real.value.CVaR=CVaRT,
       empirical.estimation.VaR=VaRE, empirical.estimation.CVaR=CVaRE,
       bootstrap.VaR=VaRB, bootstrap.CVaR=CVaRB,
       parametric.VaR=VaRP, parametric.CVaR=CVaRP,
       EVB.VaR=mean(draw[,3]), EVB.CVaR=mean(draw[,4]),
       parametric.Bayesian.VaR=mean(draw[,7]), parametric.Bayesian.CVaR=mean(draw[,8]),
       IPB.VaR=mean(draw[,11]), IPB.CVaR=mean(draw[,12]))
}