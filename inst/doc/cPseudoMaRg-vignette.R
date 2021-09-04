## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- collapse = TRUE---------------------------------------------------------
realTheta1 <- .2 + .3
realTheta2 <- .2
realParams <- c(realTheta1, realTheta2)
numObs <- 10
realX <- rnorm(numObs, mean = 0, sd = sqrt(realTheta2))
realY <- rnorm(numObs, mean = realX, sd = sqrt(realTheta1 - realTheta2))

## ---- collapse=TRUE-----------------------------------------------------------
library(cPseudoMaRg)
numImportanceSamps <- 1000
numMCMCIters <- 1000
randomWalkScale <- 1.5
recordEveryTh <- 1

# create the function that performs sampling
sampler <- makeCPMSampler(
  paramKernSamp = function(params){
    params + rnorm(2)*randomWalkScale
  },
  logParamKernEval = function(oldTheta, newTheta){
    dnorm(newTheta[1], oldTheta[1], sd = randomWalkScale, log = TRUE)
    + dnorm(newTheta[2], oldTheta[2], sd = randomWalkScale, log = TRUE)
  },
  logPriorEval = function(theta){
    if( (50 > theta[1]) & (theta[1] > theta[2]) & (theta[2] > 0) ){
      -7.130899 # - log of 50^2/2
    }else{
      -Inf
    }
  },
  logLikeApproxEval = function(y, thetaProposal, uProposal){
    if(  (50 > thetaProposal[1]) & (thetaProposal[1] > thetaProposal[2]) & (thetaProposal[2] > 0)  ){
      xSamps <- uProposal*sqrt(thetaProposal[2])
      logCondLikes <- sapply(xSamps,
                             function(xsamp) {
                               sum(dnorm(y,
                                         xsamp, 
                                         sqrt(thetaProposal[1] - thetaProposal[2]),
                                         log = TRUE)) })
      m <- max(logCondLikes)
      log(sum(exp(logCondLikes - m))) + m - log(length(y))
    }else{
      -Inf
    }
  },
  realY, 
  numImportanceSamps, 
  numMCMCIters, 
  .99, # change to 0 for original pseudo-marginal method
  recordEveryTh)

## ---- collapse=TRUE, out.width='80%'------------------------------------------
res <- sampler(realParams)
print(res)
plot(res)

## ---- collapse=TRUE, out.width='80%'------------------------------------------
samplerExact <- makeCPMSampler(
  paramKernSamp = function(params){
    return(params + rnorm(2)*randomWalkScale)
  },
  logParamKernEval = function(oldTheta, newTheta){
    dnorm(newTheta[1], oldTheta[1], sd = randomWalkScale, log = TRUE)
    + dnorm(newTheta[2], oldTheta[2], sd = randomWalkScale, log = TRUE)
  },
  logPriorEval = function(theta){
    if( (50 > theta[1]) & (theta[1] > theta[2]) & (theta[2] > 0) ){
      -7.130899 # - log of 50^2/2
    }else{
      -Inf
    }
  },
  logLikeApproxEval = function(y, thetaProposal, uProposal){
    # this is exact now!
    if( (50 > thetaProposal[1]) & (thetaProposal[1] > thetaProposal[2]) & (thetaProposal[2] > 0) ){
      sum(dnorm(y, mean = 0, sd = sqrt(thetaProposal[1]), log = TRUE))
    }else{
      -Inf
    }
  },
  realY, 
  numImportanceSamps, # doesn't this matter because Us are not used
  numMCMCIters, 
  .99, # doesn't this matter because Us are not used
  recordEveryTh)
res2 <- samplerExact(realParams)
print(res2)
plot(res2)

## ---- collapse=TRUE, out.width='80%', eval = FALSE----------------------------
#  library(cPseudoMaRg)
#  devtools::install_github("tbrown122387/pfexamplesinr@e4e2a80")
#  library(pfexamplesinr)
#  
#  returnsData <- read.csv("data/return_data.csv", header=F)[,1]
#  numParticles <- 500 # THIS MUST MATCH "#define NP 500" in src/likelihoods.cpp
#  numMCMCIters <- 500
#  randomWalkScale <- .1
#  recordEveryTh <- 1
#  numUs <- length(returnsData)*(numParticles+1)
#  
#  # some helper functions
#  transformParams <- function(untrans){
#    p <- vector(mode = "numeric", length = 3)
#    p[1] <- boot::logit(.5*(untrans[1] + 1))
#    p[2] <- untrans[2]
#    p[3] <- log(untrans[3])
#    return(p)
#  }
#  revTransformParams <- function(trans){
#    p <- vector(mode = "numeric", length = 3)
#    p[1] <- 2*boot::inv.logit( trans[1] )-1
#    p[2] <- trans[2]
#    p[3] <- exp(trans[3])
#    return(p)
#  }
#  
#  sampler <- makeCPMSampler(
#    paramKernSamp = function(params){
#      revTransformParams(transformParams(params) + rnorm(3)*randomWalkScale)
#    },
#    logParamKernEval = function(oldTheta, newTheta){
#      unconstrainedNew <- transformParams(newTheta)
#      unconstrainedOld <- transformParams(oldTheta)
#      dnorm(unconstrainedNew[1], unconstrainedOld[1], sd = randomWalkScale, log = TRUE) #phi
#      + 0.6931472 + unconstrainedNew[1] - 2*log(1 + exp(unconstrainedNew[1])) # phi jacobian
#      + dnorm(unconstrainedNew[2], unconstrainedOld[2], sd = randomWalkScale, log = TRUE) #beta
#      + dnorm(unconstrainedNew[3], unconstrainedOld[3], sd = randomWalkScale, log = TRUE) #sigmaSquared
#      - unconstrainedNew[3] # jacobian
#    },
#    logPriorEval = function(theta){
#      if( (abs(theta[1]) >= 1.0) || theta[3] <= 0.0 ){
#        -Inf
#      }else{
#         log(.5) +
#          dnorm(theta[2], mean = 0, sd = 10, log = T) +
#          dgamma(x = 1/theta[3], shape = 1.3, rate = .3, log = T)
#      }
#    },
#    logLikeApproxEval = svolApproxLL, # c++ function from pfexamplesinr
#    returnsData, numUs, numMCMCIters, .99, recordEveryTh, FALSE
#  )

## ---- collapse=TRUE, out.width='80%', eval = FALSE----------------------------
#  svolSampleResults <- sampler( c(.9, 1, .1))
#  mean(svolSampleResults)
#  print(svolSampleResults)

