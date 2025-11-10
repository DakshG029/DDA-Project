set.seed(101)
Normal.VaR.CVaR=function(sdX=1, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qnorm(pu,sd=sdX)
  X <- rnorm(n,0,sdX)
  print(head(X))
  Xu <- X[X>u]-u
  MXu=max(Xu)
  m <- length(Xu)
  
  # Real value
  VaRT <- sdX*qnorm(pk)
  CVaRT <- sdX*dnorm(qnorm(pk))/(1-pk)
  
  ## Empiric
  VaRE <- quantile(X,pk)[[1]]
  CVaRE <- mean(X[X>=VaRE])
  
  ## Parametric
  VaRP <- mean(X)+qnorm(pk)*sd(X)
  CVaRP <- -mean(X)+sd(X)*dnorm(qnorm(1-pk))/(1-pk)
  
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
  
  kMHi <- 0.15 
  deltaMHi <- sdX/kMHi
  a1 <- 0.1; b1 <- 0.1; sdk <- 0.1
  a2 <- 0.1; b2 <- 0.1; sddelta <- 0.1
  a3 <- 1; b3 <- 1
  mu_xiIBD <- -0.7+0.61*pu
  sd_xiIBD <- 0.03
  mu_sigmaIBD <- (0.34+3.18*(1-pu)-12.4*(1-pu)^2)
  sd_sigmaIBD <- exp(-41.58+83.55*pu-46.24*pu^2)
  xiIBDi <- rnorm(1,mu_xiIBD,sd_xiIBD); sdxi <- 0.03
  sigmaIBDi <- rnorm(1,mu_sigmaIBD*sd(X),sd_sigmaIBD)
  while(MXu>(-sigmaIBDi/xiIBDi)||sigmaIBDi<0||xiIBDi>0)
  {
    xiIBDi <- rnorm(1,mu_xiIBD,sd_xiIBD)
    sigmaIBDi <- rnorm(1,mu_sigmaIBD*sd(X),sd_sigmaIBD)
  }
  sdsigma <- 1
  for(i in 1:nburn)
  {
    #MH
    kMHcan <- rnorm(1,kMHi,sdk)
    while(kMHcan<0) kMHcan <- rnorm(1,kMHi,sdk)
    deltaMHcan <- rnorm(1,deltaMHi,sddelta)
    while(deltaMHcan<MXu) deltaMHcan <- rnorm(1,deltaMHi,sddelta)
    pkMH <- (a1-m-1)*log(kMHcan/kMHi)+b1*(kMHi-kMHcan) + 
      (1/kMHcan-1/kMHi)*sum(log(1-Xu/deltaMHi))
    pdeltaMH <- (a2-m-1)*log(deltaMHcan/deltaMHi) +           
      b2*(deltaMHi-deltaMHcan) + 
      (1/kMHi-1)*sum(log(1-Xu/deltaMHcan)) - 
      (1/kMHi-1)*sum(log(1-Xu/deltaMHi))
    v <- log(runif(2))
    if(v[1]<pkMH) kMHi <- kMHcan
    if(v[2]<pdeltaMH) deltaMHi <- deltaMHcan
    xiMHi <- -kMHi
    sigmaMHi <- kMHi*deltaMHi
    VaRMH <- sigmaMHi/xiMHi*((1-pcola)^(-xiMHi)-1) + u
    CVaRMH <- sigmaMHi/xiMHi*(((1-pcola)^(-xiMHi))/(1-xiMHi)-1) + u
    
    #BD
    sdBDi <- sqrt(1/rgamma(1,a3+n/2,b2+sum(X^2)/2))
    
    VaRBD <- qnorm(pk)*sdBDi
    CVaRBD <- sdBDi*dnorm(qnorm(pk))/(1-pk)
    
    #IBD
    xiIBDcan <- rnorm(1,xiIBDi,sdxi)
    sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
    while(MXu>(-sigmaIBDcan/xiIBDcan)||sigmaIBDcan<0||xiIBDcan>0)
    {
      xiIBDcan <- rnorm(1,xiIBDi,sdxi)
      sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
    }
    mu_sigma <- mu_sigmaIBD*sdBDi
    pIBD <- m*log(sigmaIBDi/sigmaIBDcan) -
      0.5/sd_xiIBD^2 * ((xiIBDcan-mu_xiIBD)^2 - (xiIBDi-mu_xiIBD)^2) -
      0.5/sd_sigmaIBD^2 * ((sigmaIBDcan-mu_sigma)^2 - (sigmaIBDi - mu_sigma)^2) -
      (1+1/xiIBDcan)*sum(log(1+xiIBDcan*Xu/sigmaIBDcan)) +
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
            sdBDi, VaRBD, CVaRBD,
            xiIBDi, sigmaIBDi, VaRIBD, CVaRIBD)
  
  for(i in 1:ndraw){
    for(j in 1:nthin)
    {
      #MH
      kMHcan <- rnorm(1,kMHi,sdk)
      while(kMHcan<0) kMHcan <- rnorm(1,kMHi,sdk)
      deltaMHcan <- rnorm(1,deltaMHi,sddelta)
      while(deltaMHcan<MXu) deltaMHcan <- rnorm(1,deltaMHi,sddelta)
      pkMH <- (a1-m-1)*log(kMHcan/kMHi)+b1*(kMHi-kMHcan) + 
        (1/kMHcan-1/kMHi)*sum(log(1-Xu/deltaMHi))
      pdeltaMH <- (a2-m-1)*log(deltaMHcan/deltaMHi)+b2*(deltaMHi-deltaMHcan) + 
        (1/kMHi-1)*sum(log(1-Xu/deltaMHcan))-(1/kMHi-1)*sum(log(1-Xu/deltaMHi))
      v <- log(runif(2))
      if(v[1]<pkMH) kMHi <- kMHcan
      if(v[2]<pdeltaMH) deltaMHi <- deltaMHcan
      xiMHi <- -kMHi; sigmaMHi <- kMHi*deltaMHi
      VaRMH <- sigmaMHi/xiMHi*((1-pcola)^(-xiMHi)-1) + u
      CVaRMH <- sigmaMHi/xiMHi*(((1-pcola)^(-xiMHi))/(1-xiMHi)-1) + u
      
        #BD
        sdBDi <- sqrt(1/rgamma(1,a3+n/2,b2+sum(X^2)/2))
        VaRBD <- qnorm(pk)*sdBDi
        CVaRBD <- sdBDi*dnorm(qnorm(pk))/(1-pk)
      
      #IBD
      xiIBDcan <- rnorm(1,xiIBDi,sdxi)
      sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
      while(MXu>=(-sigmaIBDcan/xiIBDcan)||sigmaIBDcan<=0||xiIBDcan>=0){
        xiIBDcan <- rnorm(1,xiIBDi,sdxi)
        sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
      }
      mu_sigma <- mu_sigmaIBD*sdBDi
      pIBD <- m*log(sigmaIBDi/sigmaIBDcan) -
        0.5/sd_xiIBD^2 * ((xiIBDcan-mu_xiIBD)^2 - (xiIBDi-mu_xiIBD)^2) -
        0.5/sd_sigmaIBD^2 * ((sigmaIBDcan-mu_sigma)^2 - (sigmaIBDi - mu_sigma)^2) -
        (1+1/xiIBDcan)*sum(log(1+xiIBDcan*Xu/sigmaIBDcan)) +
        (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
      if(log(runif(1))<pIBD){
        xiIBDi <- xiIBDcan
        sigmaIBDi <- sigmaIBDcan
      }
      VaRIBD <- sigmaIBDi/xiIBDi*(-1+(1-pcola)^(-xiIBDi)) + u 
      CVaRIBD <- sigmaIBDi/xiIBDi*(((1-pcola)^(-xiIBDi))/(1-xiIBDi)-1) + u  
    }
    draw <- c(draw,
              c(xiMHi, sigmaMHi, VaRMH, CVaRMH,
                sdBDi, VaRBD, CVaRBD,
                xiIBDi, sigmaIBDi, VaRIBD, CVaRIBD))
  }
  
  draw <- matrix(draw, ncol=11, byrow = TRUE)
  list(standar.deviation=sdX, length=n, threshold=pu, p=pk, 
       real.value.VaR=VaRT, real.value.CVaR=CVaRT,
       empirical.estimation.VaR=VaRE, empirical.estimation.CVaR=CVaRE,
       bootstrap.VaR=VaRB, bootstrap.CVaR=CVaRB,
       parametric.VaR=VaRP, parametric.CVaR=CVaRP,
       EVB.VaR=mean(draw[,3]), EVB.CVaR=mean(draw[,4]),
       parametric.Bayesian.VaR=mean(draw[,6]), parametric.Bayesian.CVaR=mean(draw[,7]),
       IPB.VaR=mean(draw[,10]), IPB.CVaR=mean(draw[,11]))
}

