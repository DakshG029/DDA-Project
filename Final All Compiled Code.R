set.seed(101010)


#Normal
Normal.VaR.CVaR=function(sdX=1, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qnorm(pu,sd=sdX)
  X <- rnorm(n,0,sdX)
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

Normal.VaR.CVaR()

#Exponential
Exp.VaR.CVaR=function(lambdaX=1, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qexp(pu,lambdaX)
  X <- rexp(n,lambdaX)
  Xu <- X[X>u]-u
  MXu=max(Xu)
  m <- length(Xu)
  
  # Real value
  VaRT <- -log(1-pk)/lambdaX
  CVaRT <- (1-log(1-pk))/lambdaX
  
  ## Empiric
  VaRE <- quantile(X,pk)[[1]]
  CVaRE <- mean(X[X>=VaRE])
  
  ## Parametric
  VaRP <- -log(1-pk)*mean(X)
  CVaRP <- mean(X)*(1-log(1-pk))
  
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
  
  a1 <- 1; b1 <- 1 
  a2 <- 1; b2 <- 1
  sigmaIBDi <- mean(X)
  for(i in 1:nburn)
  {
    #MH
    sigmaMHi <- 1/rgamma(1,a1+m,b1+sum(Xu))
    VaRMH <- -log(1-pcola)*sigmaMHi + u 
    CVaRMH <- (1-log(1-pcola))*sigmaMHi + u 
    
    #BD
    lambdaBDi <- rgamma(1,a2+n,b2+sum(X))
    VaRBD <- -log(1-pk)/lambdaBDi
    CVaRBD <- (1-log(1-pk))/lambdaBDi
    
    #IBD
    sigmaIBDcan <- 1/lambdaBDi
    pIBD <- m*log(sigmaIBDi/sigmaIBDcan) + sum(Xu)*(1/sigmaIBDi-1/sigmaIBDcan)
    if(log(runif(1))<pIBD) sigmaIBDi <- sigmaIBDcan
    VaRIBD <- -log(1-pcola)*sigmaIBDi + u
    CVaRIBD <- (1-log(1-pcola))*sigmaIBDi + u 
  }
  draw <- c(sigmaMHi, VaRMH, CVaRMH,
            lambdaBDi, VaRBD, CVaRBD,
            sigmaIBDi, VaRIBD, CVaRIBD)
  
  for(i in 1:ndraw)
  {
    for(j in 1:nthin)
    {
      #MH
      sigmaMHi <- 1/rgamma(1,a1+m,b1+sum(Xu))
      VaRMH <- -log(1-pcola)*sigmaMHi + u
      CVaRMH <- (1-log(1-pcola))*sigmaMHi + u 
      #BD
      lambdaBDi <- rgamma(1,a2+n,b2+sum(X))
      VaRBD <- -log(1-pk)/lambdaBDi
      CVaRBD <- (1-log(1-pk))/lambdaBDi
      
      #IBD
      sigmaIBDcan <- 1/lambdaBDi
      pIBD <- m*log(sigmaIBDi/sigmaIBDcan) +sum(Xu)*(1/sigmaIBDi-1/sigmaIBDcan)
      if(log(runif(1))<pIBD) sigmaIBDi <- sigmaIBDcan
      VaRIBD <- -log(1-pcola)*sigmaIBDi + u
      CVaRIBD <- (1-log(1-pcola))*sigmaIBDi + u 
    }
    draw <- c(draw,c(sigmaMHi, VaRMH, CVaRMH,
                     lambdaBDi, VaRBD, CVaRBD,
                     sigmaIBDi, VaRIBD, CVaRIBD))
  }
  draw <- matrix(draw,ncol=9,byrow = TRUE)
  
  list(lambda=lambdaX, length=n, threshold=pu, p=pk, 
       real.value.VaR=VaRT, real.value.CVaR=CVaRT,
       empirical.estimation.VaR=VaRE, empirical.estimation.CVaR=CVaRE,
       bootstrap.VaR=VaRB, bootstrap.CVaR=CVaRB,
       parametric.VaR=VaRP, parametric.CVaR=CVaRP,
       EVB.VaR=mean(draw[,2]), EVB.CVaR=mean(draw[,3]),
       parametric.Bayesian.VaR=mean(draw[,5]), parametric.Bayesian.CVaR=mean(draw[,6]),
       IPB.VaR=mean(draw[,8]), IPB.CVaR=mean(draw[,9]))
}

Exp.VaR.CVaR()

#Cauchy
Cauchy.VaR.CVaR=function(deltaX=1, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qcauchy(pu,0,deltaX)
  X <- rcauchy(n,0,deltaX)
  Xu <- X[X>u]-u
  MXu=max(Xu)
  m <- length(Xu)
  
  # Real value
  VaRT <- deltaX*tan(pi*(pk-0.5))
  
  ## Empiric
  VaRE <- quantile(X,pk)[[1]]
  
  ## Parametric
  VaRP <- IQR(X)*tan(pi*(pk-0.5))/2
  
  ## Bootstrap
  B <- 1000
  MR <- function(x)
  {
    rX <- sample(X,n,rep=T)
    quantile(rX,pk)[[1]]
  }
  EstR <- mean(replicate(B,MR(sample(X,n,rep=T))))
  VaRB <- 2*VaRE - EstR
  
  ## MCMC
  pcola <- 1-(1-pk)/(1-pu)
  nburn <- 3000
  nthin <- 50
  ndraw <- 10000
  
  xiMHi <- 1
  a1 <- 10^(-6); b1 <- 10^(-3)
  sd_xiMH <- 0.01
  a2 <- 1; b2 <- 1
  sigmaMHi <- 1/rgamma(1,a2,b2); sd_sigmaMH <- 1
  deltaBDi <- as.numeric((quantile(X,3/4)-quantile(X,1/4))/2)
  a3 <- 1; b3 <- 1; sd_deltaBD <- 1
  mu_xiIBD <- 1; sd_xiIBD <- 0.065
  mu_sigmaIBD <- 1/(pi*(1-pu))
  sd_sigmaIBD <- exp(266.13-588.51*pu+323.57*pu^2)
  xiIBDi <- rnorm(1,mu_xiIBD,sd_xiIBD); 
  sdxi <- 0.1
  while(xiIBDi<0) xiIBDi <- rnorm(1,mu_xiIBD,sd_xiIBD)
  sigmaIBDi <- rnorm(1,mu_sigmaIBD*deltaBDi,sd_sigmaIBD)
  while(sigmaIBDi<0) sigmaIBDi <- rnorm(1,mu_sigmaIBD*deltaBDi,sd_sigmaIBD)
  sdsigma <- 0.1
  for(i in 1:nburn)
  {
    #MH
    xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
    while(xiMHcan<b1) xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
    sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
    while(sigmaMHcan<0) sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH) 
    pxiMH <- (a1+1)*log(xiMHi/xiMHcan) + (1+1/xiMHi)*sum(log(1+xiMHi*Xu/sigmaMHi)) - 
      (1+1/xiMHcan)*sum(log(1+xiMHcan*Xu/sigmaMHi))
    psigmaMH <- (m+a2+1)*log(sigmaMHi/sigmaMHcan) + 
      b2*(1/sigmaMHi - 1/sigmaMHcan) + 
      (1+1/xiMHi)*sum(log(1+xiMHi*Xu/sigmaMHi)) -
      (1+1/xiMHi)*sum(log(1+xiMHi*Xu/sigmaMHcan))
    v <- log(runif(2))
    if(v[1]<pxiMH) xiMHi <- xiMHcan
    if(v[2]<psigmaMH) sigmaMHi <- sigmaMHcan
    VaRMH <- sigmaMHi/xiMHi*((1-pcola)^(-xiMHi)-1) + u
    
    #BD
    deltaBDcan <- rnorm(1,deltaBDi,sd_deltaBD)
    while(deltaBDcan<0) deltaBDcan <- rnorm(1,deltaBDi,sd_deltaBD)
    pBD <- (a3-n-1)*log(deltaBDcan/deltaBDi) - b3*(deltaBDcan-deltaBDi) - 
      sum(log(1+X^2/deltaBDcan^2)) + sum(log(1+X^2/deltaBDi^2))
    if(log(runif(1))<pBD) deltaBDi <- deltaBDcan
    VaRBD <- deltaBDi*tan(pi*(pk-0.5))
    
    #IBD
    xiIBDcan <- rnorm(1,xiIBDi,sdxi)
    while(xiIBDcan<0) xiIBDcan <- rnorm(1,xiIBDi,sdxi)
    sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
    while(sigmaIBDcan<0) sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
    mu_sigma <- mu_sigmaIBD*deltaBDi
    pxiIBD <- - 0.5/sd_xiIBD^2 * ((xiIBDcan-mu_xiIBD)^2 - (xiIBDi-mu_xiIBD)^2) -
      (1+1/xiIBDcan)*sum(log(1+xiIBDcan*Xu/sigmaIBDi)) + 
      (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
    psigmaIBD <- m*log(sigmaIBDi/sigmaIBDcan) - 
      0.5/sd_sigmaIBD^2 * ((sigmaIBDcan-mu_sigma)^2 - (sigmaIBDi-mu_sigma)^2) -
      (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDcan)) +
      (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
    v <- log(runif(2))
    if(v[1]<pxiIBD) xiIBDi <- xiIBDcan
    if(v[2]<psigmaIBD) sigmaIBDi <- sigmaIBDcan
    VaRIBD <- sigmaIBDi/xiIBDi*((1-pcola)^(-xiIBDi)-1) + u
  }
  draw <- c(xiMHi, sigmaMHi, VaRMH,
            deltaBDi, VaRBD,
            xiIBDi, sigmaIBDi, VaRIBD)
  for(i in 1:ndraw)
  {
    for(j in 1:nthin)
    {
      #MH
      xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
      while(xiMHcan<b1) xiMHcan <- rnorm(1,xiMHi,sd_xiMH)
      sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
      while(sigmaMHcan<0) sigmaMHcan <- rnorm(1,sigmaMHi,sd_sigmaMH)
      pxiMH <- (a1+1)*log(xiMHi/xiMHcan) + 
        (1+1/xiMHi)*sum(log(1+xiMHi*Xu/sigmaMHi)) - 
        (1+1/xiMHcan)*sum(log(1+xiMHcan*Xu/sigmaMHi))
      psigmaMH <- (m+a2+1)*log(sigmaMHi/sigmaMHcan) + 
        b2*(1/sigmaMHi - 1/sigmaMHcan) + 
        (1+1/xiMHi)*sum(log(1+xiMHi*Xu/sigmaMHi)) -
        (1+1/xiMHi)*sum(log(1+xiMHi*Xu/sigmaMHcan))
      v <- log(runif(2))
      if(v[1]<pxiMH) xiMHi <- xiMHcan
      if(v[2]<psigmaMH) sigmaMHi <- sigmaMHcan
      VaRMH <- sigmaMHi/xiMHi*((1-pcola)^(-xiMHi)-1) + u
      #BD
      deltaBDcan <- rnorm(1,deltaBDi,sd_deltaBD)
      while(deltaBDcan<0) deltaBDcan <- rnorm(1,deltaBDi,sd_deltaBD)
      pBD <- (a3-n-1)*log(deltaBDcan/deltaBDi) -
        b3*(deltaBDcan-deltaBDi) - 
        sum(log(1+X^2/deltaBDcan^2)) + 
        sum(log(1+X^2/deltaBDi^2))
      if(log(runif(1))<pBD) deltaBDi <- deltaBDcan
      VaRBD <- deltaBDi*tan(pi*(pk-0.5))
      
      #IBD
      xiIBDcan <- rnorm(1,xiIBDi,sdxi)
      while(xiIBDcan<0) xiIBDcan <- rnorm(1,xiIBDi,sdxi)
      sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
      while(sigmaIBDcan<0) sigmaIBDcan <- rnorm(1,sigmaIBDi,sdsigma)
      mu_sigma <- mu_sigmaIBD*deltaBDi
      pxiIBD <- - 0.5/sd_xiIBD^2 * ((xiIBDcan-mu_xiIBD)^2 - (xiIBDi-mu_xiIBD)^2) -
        (1+1/xiIBDcan)*sum(log(1+xiIBDcan*Xu/sigmaIBDi)) +
        (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
      psigmaIBD <- m*log(sigmaIBDi/sigmaIBDcan) - 
        0.5/sd_sigmaIBD^2 * ((sigmaIBDcan-mu_sigma)^2 - (sigmaIBDi-mu_sigma)^2) -
        (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDcan)) +
        (1+1/xiIBDi)*sum(log(1+xiIBDi*Xu/sigmaIBDi))
      v <- log(runif(2))
      if(v[1]<pxiIBD) xiIBDi <- xiIBDcan
      if(v[2]<psigmaIBD) sigmaIBDi <- sigmaIBDcan
      VaRIBD <- sigmaIBDi/xiIBDi*((1-pcola)^(-xiIBDi)-1) + u
    }
    draw <- c(draw, c(xiMHi, sigmaMHi, VaRMH,
                      deltaBDi, VaRBD,
                      xiIBDi, sigmaIBDi, VaRIBD))
  }
  draw <- matrix(draw,ncol=8,byrow = TRUE)
  
  list(delta=deltaX, length=n, threshold=pu, p=pk, 
       real.value.VaR=VaRT, 
       empirical.estimation.VaR=VaRE, 
       bootstrap.VaR=VaRB, 
       parametric.VaR=VaRP, 
       EVB.VaR=mean(draw[,3]), 
       parametric.Bayesian.VaR=mean(draw[,5]), 
       IPB.VaR=mean(draw[,8]))
}

Cauchy.VaR.CVaR()

#Gamma
Gamma.VaR.CVaR=function(alphaX=2, betaX=2, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qgamma(pu,alphaX,betaX)
  X <- rgamma(n,alphaX,betaX)
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

library(gsl)
Gamma.VaR.CVaR()

set.seed(101)

# ------------------------------
# 1) Load libraries
# ------------------------------
library(quantmod)
library(fitdistrplus)
library(evir)

# ------------------------------
# 2) Prepare data (NIFTY50)
# ------------------------------
ticker <- "^NSEI"
start_date <- as.Date("2010-01-01")
end_date <- Sys.Date()
dataset_xts <- getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
dataset <- data.frame(Date = index(dataset_xts), coredata(dataset_xts))
dataset <- na.omit(dataset)

prices <- dataset$NSEI.Close
returns <- diff(log(prices))
losses <- -as.numeric(returns)    # working in LOSSES
losses <- losses[!is.na(losses)]

pk <- 0.95      # confidence level for VaR/CVaR
pu <- 0.90      # tail threshold
pu_thresh <- quantile(losses, pu)
exceedances <- losses[losses > pu_thresh] - pu_thresh

# ------------------------------
# 3) Helper: GPD log-likelihood + VaR formula
# ------------------------------
gpd_loglik <- function(xi, beta, y) {
  if (beta <= 0) return(-Inf)
  n <- length(y)
  if (abs(xi) < 1e-8) {
    ll <- -n * log(beta) - sum(y)/beta
  } else {
    if (any(1 + xi * y / beta <= 0)) return(-Inf)
    ll <- -n * log(beta) - (1/xi + 1) * sum(log(1 + xi * y / beta))
  }
  return(ll)
}

gpd_VaR_from_params <- function(threshold, xi, beta, pu, pk) {
  q <- (pk - pu)/(1 - pu)
  if (abs(xi) < 1e-8) {
    return(threshold - beta * log(1 - q))
  } else {
    return(threshold + beta/xi * ((1 - q)^(-xi) - 1))
  }
}

# ------------------------------
# 4) Method 1: Empirical
# ------------------------------
VaR_emp <- quantile(losses, pk)
CVaR_emp <- mean(losses[losses >= VaR_emp])

# ------------------------------
# 5) Method 2: Bootstrap
# ------------------------------
B <- 500
VaR_boot <- numeric(B)
CVaR_boot <- numeric(B)
for (i in 1:B) {
  sample_X <- sample(losses, length(losses), replace = TRUE)
  VaR_boot[i] <- quantile(sample_X, pk)
  CVaR_boot[i] <- mean(sample_X[sample_X >= VaR_boot[i]])
}
VaR_bootstrap <- 2 * VaR_emp - mean(VaR_boot)
CVaR_bootstrap <- 2 * CVaR_emp - mean(CVaR_boot)

# ------------------------------
# 6) Method 3: Parametric (t-distribution) [FIXED]
# ------------------------------
library(MASS)

# Fit a Student-t distribution using MASS::fitdistr
fit_t <- fitdistr(losses, densfun = "t", start = list(m = mean(losses), s = sd(losses), df = 5))

m_t <- fit_t$estimate["m"]
s_t <- fit_t$estimate["s"]
df_t <- fit_t$estimate["df"]

VaR_par <- m_t + s_t * qt(pk, df = df_t)
CVaR_par <- m_t + s_t * (dt(qt(pk, df_t), df_t) / (1 - pk)) * ((df_t + qt(pk, df_t)^2) / (df_t - 1))


# ------------------------------
# 7) Method 4: Parametric Bayesian (t-dist priors)
# ------------------------------
n_iter <- 5000
burn_in <- 1000
pk <- 0.95  # confidence level

# Initialize chains
mu_chain <- numeric(n_iter)
sigma_chain <- numeric(n_iter)
df_chain <- numeric(n_iter)

# Starting values
mu_chain[1] <- mean(losses)
sigma_chain[1] <- sd(losses)
df_chain[1] <- 5

# MCMC sampling
for (i in 2:n_iter) {
  # --- Proposal step ---
  mu_prop <- rnorm(1, mu_chain[i-1], 0.02)
  sigma_prop <- abs(rnorm(1, sigma_chain[i-1], 0.02))
  df_prop <- abs(rnorm(1, df_chain[i-1], 0.5))
  
  # --- Log-likelihood ---
  ll_curr <- sum(log(dt((losses - mu_chain[i-1]) / sigma_chain[i-1],
                        df_chain[i-1]) / sigma_chain[i-1]))
  ll_prop <- sum(log(dt((losses - mu_prop) / sigma_prop,
                        df_prop) / sigma_prop))
  
  # --- Log-prior ---
  prior_curr <- dnorm(mu_chain[i-1], mean(losses), sd(losses), log = TRUE) +
    dgamma(sigma_chain[i-1], 2, 2, log = TRUE)
  prior_prop <- dnorm(mu_prop, mean(losses), sd(losses), log = TRUE) +
    dgamma(sigma_prop, 2, 2, log = TRUE)
  
  # --- MH acceptance ratio ---
  log_acc <- (ll_prop + prior_prop) - (ll_curr + prior_curr)
  
  if (log(runif(1)) < log_acc) {
    mu_chain[i] <- mu_prop
    sigma_chain[i] <- sigma_prop
    df_chain[i] <- df_prop
  } else {
    mu_chain[i] <- mu_chain[i-1]
    sigma_chain[i] <- sigma_chain[i-1]
    df_chain[i] <- df_chain[i-1]
  }
}

# -------------------------------
# Burn-in and posterior estimates
# -------------------------------
mu_chain <- mu_chain[(burn_in + 1):n_iter]
sigma_chain <- sigma_chain[(burn_in + 1):n_iter]
df_chain <- df_chain[(burn_in + 1):n_iter]

mu_post <- mean(mu_chain)
sigma_post <- mean(sigma_chain)
df_post <- mean(df_chain)

# -------------------------------
# VaR and CVaR computation
# -------------------------------
VaR_parB <- mu_post + sigma_post * qt(pk, df_post)

CVaR_parB <- mu_post + sigma_post *
  (dt(qt(pk, df_post), df_post) / (1 - pk)) *
  ((df_post + qt(pk, df_post)^2) / (df_post - 1))
# ------------------------------
# 8) Method 5: Extreme Value Bayesian (EVB)
# ------------------------------
gpd_fit <- gpd(losses, pu_thresh)
xi_evb <- gpd_fit$par.ests["xi"]
beta_evb <- gpd_fit$par.ests["beta"]

VaR_evb <- gpd_VaR_from_params(pu_thresh, xi_evb, beta_evb, pu, pk)
q <- (pk - pu)/(1 - pu)
CVaR_evb <- pu_thresh + (beta_evb / (1 - xi_evb)) * ((1 - q)^(-xi_evb) - 1)

# ------------------------------
# 9) Method 6: Informative Prior Bayesian (IPB)
# ------------------------------
mu_xi_prior <- 0.1; sd_xi_prior <- 0.06
mean_beta_prior <- max(mean(exceedances, na.rm = TRUE), 1e-6)
sd_beta_prior <- 0.5 * mean_beta_prior

n_mcmc2 <- 8000; burn_mcmc2 <- 3000
xi_chain2 <- numeric(n_mcmc2); beta_chain2 <- numeric(n_mcmc2)
xi_chain2[1] <- mu_xi_prior
beta_chain2[1] <- mean_beta_prior
sd_xi_prop2 <- 0.03
sd_beta_prop2 <- 0.1 * beta_chain2[1]
n_sim_tail <- 5000

for (i in 2:n_mcmc2) {
  xi_prop <- rnorm(1, xi_chain2[i-1], sd_xi_prop2)
  beta_prop <- rnorm(1, beta_chain2[i-1], sd_beta_prop2)
  if (beta_prop <= 0) { xi_chain2[i] <- xi_chain2[i-1]; beta_chain2[i] <- beta_chain2[i-1]; next }
  
  logprior_curr <- dnorm(xi_chain2[i-1], mu_xi_prior, sd_xi_prior, log = TRUE) -
    dgamma(beta_chain2[i-1], shape = (mean_beta_prior/sd_beta_prior)^2, rate = mean_beta_prior/(sd_beta_prior^2), log = TRUE)
  logprior_prop <- dnorm(xi_prop, mu_xi_prior, sd_xi_prior, log = TRUE) -
    dgamma(beta_prop, shape = (mean_beta_prior/sd_beta_prior)^2, rate = mean_beta_prior/(sd_beta_prior^2), log = TRUE)
  
  ll_curr <- gpd_loglik(xi_chain2[i-1], beta_chain2[i-1], exceedances)
  ll_prop <- gpd_loglik(xi_prop, beta_prop, exceedances)
  
  logacc <- (ll_prop + logprior_prop) - (ll_curr + logprior_curr)
  if (log(runif(1)) < logacc) {
    xi_chain2[i] <- xi_prop; beta_chain2[i] <- beta_prop
  } else {
    xi_chain2[i] <- xi_chain2[i-1]; beta_chain2[i] <- beta_chain2[i-1]
  }
}

xi_post2 <- xi_chain2[(burn_mcmc2+1):n_mcmc2]
beta_post2 <- beta_chain2[(burn_mcmc2+1):n_mcmc2]
n_post2 <- length(xi_post2)

VaR_ibd_draws <- numeric(n_post2)
CVaR_ibd_draws <- numeric(n_post2)
for (i in 1:n_post2) {
  xi_i <- xi_post2[i]; beta_i <- beta_post2[i]
  VaR_i <- gpd_VaR_from_params(pu_thresh, xi_i, beta_i, pu, pk)
  U <- runif(n_sim_tail)
  if (abs(xi_i) < 1e-8) {
    Ysim <- -beta_i * log(1 - U)
  } else {
    Ysim <- beta_i/xi_i * ((1 - U)^(-xi_i) - 1)
  }
  Xsim <- pu_thresh + Ysim
  CVaR_i <- mean(Xsim[Xsim >= VaR_i])
  VaR_ibd_draws[i] <- VaR_i
  CVaR_ibd_draws[i] <- CVaR_i
}
VaR_ipb <- mean(VaR_ibd_draws, na.rm = TRUE)
CVaR_ipb <- mean(CVaR_ibd_draws, na.rm = TRUE)

# ------------------------------
# 10) Compile all results
# ------------------------------
results <- data.frame(
  Method = c("Empirical", "Bootstrap", "Parametric (t)", 
             "Parametric Bayesian", "Extreme Value Bayesian", "Informative Prior Bayesian"),
  VaR = c(VaR_emp, VaR_bootstrap, VaR_par, VaR_parB, VaR_evb, VaR_ipb),
  CVaR = c(CVaR_emp, CVaR_bootstrap, CVaR_par, CVaR_parB, CVaR_evb, CVaR_ipb)
)
print("NIFTY50")
print(results)


# ==============================
# GARCH vs IPB for VaR and CVaR estimation
# ==============================

library(quantmod)
library(rugarch)
library(ggplot2)

set.seed(101)

# ------------------------------
# 1) Load data
# ------------------------------
ticker <- "^NSEI"
start_date <- as.Date("2010-01-01")
end_date <- Sys.Date()

dataset_xts <- getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
dataset <- data.frame(Date = index(dataset_xts), coredata(dataset_xts))
dataset <- na.omit(dataset)

prices <- dataset$NSEI.Close
returns <- diff(log(prices))
losses <- -returns

n <- length(losses)
train_end <- floor(0.8 * n)
train_losses <- losses[1:train_end]
test_losses <- losses[(train_end + 1):n]

# ------------------------------
# 2) GARCH(1,1) model
# ------------------------------
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std"
)

fit_garch <- ugarchfit(spec, train_losses)

# Forecast volatility for test period
forecast <- ugarchforecast(fit_garch, n.ahead = length(test_losses))
sigma_forecast <- sigma(forecast)
shape <- coef(fit_garch)["shape"]

# VaR and CVaR from forecast
pk <- 0.95
q_pk <- qt(pk, df = shape)
VaR_GARCH <- q_pk * sigma_forecast

u_grid <- seq(pk, 1, length.out = 1000)
CVaR_GARCH <- sigma_forecast * mean(qt(u_grid, df = shape))

# ------------------------------
# 3) IPB on training data
# ------------------------------
pu <- 0.90
pu_thresh <- quantile(train_losses, pu)
exceedances <- train_losses[train_losses > pu_thresh] - pu_thresh
n_sim_tail <- 5000

# GPD log-likelihood
gpd_loglik <- function(xi, beta, y) {
  if (beta <= 0) return(-Inf)
  n <- length(y)
  if (abs(xi) < 1e-8) {
    ll <- -n * log(beta) - sum(y) / beta
  } else {
    if (any(1 + xi * y / beta <= 0)) return(-Inf)
    ll <- -n * log(beta) - (1 / xi + 1) * sum(log(1 + xi * y / beta))
  }
  return(ll)
}

# VaR from params
gpd_VaR_from_params <- function(threshold, xi, beta, pu, pk) {
  q <- (pk - pu) / (1 - pu)
  if (abs(xi) < 1e-8) {
    return(threshold - beta * log(1 - q))
  } else {
    return(threshold + beta / xi * ((1 - q)^(-xi) - 1))
  }
}

# Informative priors
mu_xi_prior <- 0.1
sd_xi_prior <- 0.06
mean_beta_prior <- max(mean(exceedances, na.rm = TRUE), 1e-6)
sd_beta_prior <- 0.5 * mean_beta_prior

# MCMC
n_mcmc <- 5000
burn <- 2000
xi_chain <- numeric(n_mcmc)
beta_chain <- numeric(n_mcmc)
xi_chain[1] <- mu_xi_prior
beta_chain[1] <- mean_beta_prior

sd_xi_prop <- 0.03
sd_beta_prop <- 0.1 * beta_chain[1]

for (i in 2:n_mcmc) {
  xi_prop <- rnorm(1, xi_chain[i - 1], sd_xi_prop)
  beta_prop <- rnorm(1, beta_chain[i - 1], sd_beta_prop)
  
  if (beta_prop <= 0) {
    xi_chain[i] <- xi_chain[i - 1]
    beta_chain[i] <- beta_chain[i - 1]
    next
  }
  
  logprior_curr <- dnorm(xi_chain[i - 1], mu_xi_prior, sd_xi_prior, log = TRUE) -
    dgamma(beta_chain[i - 1],
           shape = (mean_beta_prior / sd_beta_prior)^2,
           rate = mean_beta_prior / (sd_beta_prior^2),
           log = TRUE)
  
  logprior_prop <- dnorm(xi_prop, mu_xi_prior, sd_xi_prior, log = TRUE) -
    dgamma(beta_prop,
           shape = (mean_beta_prior / sd_beta_prior)^2,
           rate = mean_beta_prior / (sd_beta_prior^2),
           log = TRUE)
  
  ll_curr <- gpd_loglik(xi_chain[i - 1], beta_chain[i - 1], exceedances)
  ll_prop <- gpd_loglik(xi_prop, beta_prop, exceedances)
  
  logacc <- (ll_prop + logprior_prop) - (ll_curr + logprior_curr)
  
  if (log(runif(1)) < logacc) {
    xi_chain[i] <- xi_prop
    beta_chain[i] <- beta_prop
  } else {
    xi_chain[i] <- xi_chain[i - 1]
    beta_chain[i] <- beta_chain[i - 1]
  }
}

xi_post <- xi_chain[(burn + 1):n_mcmc]
beta_post <- beta_chain[(burn + 1):n_mcmc]

# Compute VaR / CVaR
VaR_IPB_draws <- numeric(length(xi_post))
CVaR_IPB_draws <- numeric(length(xi_post))

for (i in 1:length(xi_post)) {
  xi_i <- xi_post[i]
  beta_i <- beta_post[i]
  VaR_i <- gpd_VaR_from_params(pu_thresh, xi_i, beta_i, pu, pk)
  
  U <- runif(n_sim_tail)
  if (abs(xi_i) < 1e-8) {
    Ysim <- -beta_i * log(1 - U)
  } else {
    Ysim <- beta_i / xi_i * ((1 - U)^(-xi_i) - 1)
  }
  Xsim <- pu_thresh + Ysim
  CVaR_i <- mean(Xsim[Xsim >= VaR_i])
  
  VaR_IPB_draws[i] <- VaR_i
  CVaR_IPB_draws[i] <- CVaR_i
}

VaR_IPB <- mean(VaR_IPB_draws)
CVaR_IPB <- mean(CVaR_IPB_draws)

# ------------------------------
# 4) Backtesting
# ------------------------------
violations_GARCH <- mean(test_losses > VaR_GARCH)
violations_IPB <- mean(test_losses > VaR_IPB)

ES_error_GARCH <- mean((test_losses[test_losses > VaR_GARCH] - CVaR_GARCH)^2)
ES_error_IPB <- mean((test_losses[test_losses > VaR_IPB] - CVaR_IPB)^2)

# ------------------------------
# 5) Print results
# ------------------------------
cat("GARCH VaR:", mean(VaR_GARCH),
    " Violations:", violations_GARCH,
    " CVaR Error:", ES_error_GARCH, "\n")

cat("IPB VaR:", VaR_IPB,
    " Violations:", violations_IPB,
    " CVaR Error:", ES_error_IPB, "\n")

# ------------------------------
# 6) Plots
# ------------------------------
true_VaR <- quantile(test_losses, pk)
true_CVaR <- mean(test_losses[test_losses >= true_VaR])

VaR_GARCH_daily <- rep(mean(VaR_GARCH), length(test_losses))
CVaR_GARCH_daily <- rep(mean(CVaR_GARCH), length(test_losses))
VaR_IPB_daily <- rep(VaR_IPB, length(test_losses))
CVaR_IPB_daily <- rep(CVaR_IPB, length(test_losses))

test_dates <- dataset$Date[(train_end + 2):(train_end + 1 + length(test_losses))]

# ---- VaR Plot ----
plot_df_VaR <- data.frame(
  Date = test_dates,
  Daily_Losses = test_losses,
  VaR_GARCH = VaR_GARCH_daily,
  VaR_IPB = VaR_IPB_daily
)

ggplot(plot_df_VaR, aes(x = Date)) +
  geom_line(aes(y = Daily_Losses, color = "Losses"), size = 0.8) +
  geom_line(aes(y = VaR_GARCH, color = "VaR GARCH"), size = 0.8) +
  geom_line(aes(y = VaR_IPB, color = "VaR IPB"), size = 0.8) +
  geom_hline(yintercept = true_VaR, linetype = "dashed", color = "green", size = 1, alpha = 0.7) +
  scale_color_manual(values = c("Losses" = "black", "VaR GARCH" = "red", "VaR IPB" = "blue")) +
  labs(title = "Daily Losses with Predicted VaR", x = "Date", y = "Losses / Returns (%)", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ---- CVaR Plot ----
plot_df_CVaR <- data.frame(
  Date = test_dates,
  Daily_Losses = test_losses,
  CVaR_GARCH = CVaR_GARCH_daily,
  CVaR_IPB = CVaR_IPB_daily
)

ggplot(plot_df_CVaR, aes(x = Date)) +
  geom_line(aes(y = Daily_Losses, color = "Losses"), size = 0.8) +
  geom_line(aes(y = CVaR_GARCH, color = "CVaR GARCH"), size = 0.8) +
  geom_line(aes(y = CVaR_IPB, color = "CVaR IPB"), size = 0.8) +
  geom_hline(yintercept = true_CVaR, linetype = "dashed", color = "green", size = 1, alpha = 0.7) +
  scale_color_manual(values = c("Losses" = "black", "CVaR GARCH" = "darkred", "CVaR IPB" = "darkblue")) +
  labs(title = "Daily Losses with Predicted CVaR", x = "Date", y = "Losses / Returns (%)", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "bottom")
