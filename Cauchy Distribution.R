set.seed(101)
Cauchy.VaR.CVaR=function(deltaX=1, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qcauchy(pu,0,deltaX)
  X <- rcauchy(n,0,deltaX)
  print(head(X))
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