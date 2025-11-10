set.seed(101)
Exp.VaR.CVaR=function(lambdaX=1, n=10^3, B=10^3, pu=0.9, pk=0.95)
{
  u <- qexp(pu,lambdaX)
  X <- rexp(n,lambdaX)
  print(head(X))
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