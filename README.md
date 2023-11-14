# HighDimRD
This repository contains the R-Package HighDimRD for easy inclusion of high-dimensional covariates in regression discontinuity design problems. The main function for this purpose is HighDim_rd. Details of its functionality are contained in the paper "Inference in Regression Discontinuity Designs with High-Dimensional Covariates" which I have written with Christoph Rothe. The paper has appeared in The Econometrics Journal and can be found <a href="https://doi.org/10.1093/ectj/utac029">here</a>. A pre-print is available on <a href="https://arxiv.org/abs/2110.13725">arxiv</a>. The basic functionality is that in a first step the algorithm selects relevant covariates by means of a Lasso procedure. In a second step the treatment effect is estimated by using these selected covariates. In addition, there are some functions which create interactions, cross-interactions and Fourier expansions in order to generate high-dimensional covariates.

For the package to work, the packages glmnet, Matrix, rdrobust and RDHonest are required. The package and all its dependencies can be installed like this:
```
library(remotes)
install_github("kolesarm/RDHonest")
library(RDHonest)
install_github("akreiss/HighDimRD")
library(HighDimRD)

```

Here is an Example for how the package can be used: We first generate some data.
```
## Load Libraries
library(mvtnorm)

## Set seed
set.seed(12345)

## Set variables
n <- 1000
sigma_y <- 0.1295
sigma_z <- 0.1353

## Running Variable
x <- 2*rbeta(n,2,4)-1

## Covariates
Z <- rmvnorm(n,rep(0,5),sigma=diag(sigma_z^2,5))

## Noise
epsilon <- rnorm(n,mean=0,sd=sigma_y)

## Treatment Indicator
T <- as.numeric(x>=0)

## Outcome
Y <- (1-T)*(0.36+0.96*x+5.47*x^2+15.28*x^3+15.87*x^4+5.14*x^5+0.22*Z[,1])+
  T *(0.38+0.62*x-2.84*x^2+ 8.42*x^3-10.24*x^4+4.31*x^5+0.28*Z[,1])+epsilon
```
In the above example the coavariates `Z` are five-dimensional and only its first component is used. We add some Fourier expansions and interaction terms to the covariates:
```
## High-Dimensional Covariates
Z1 <- fourier_basis(Z,4)
Z_HighDim <- cbind(Z,interaction_terms(Z),Z1,interaction_terms(Z1))
```
Using these high-dimensional covariates, we perform the estimation with different methods of tuning parameter choice and different rd packages:
```
## Compute estimate with model selection using rdrobust
cov_rd_CV_robust  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV" ,rd="robust")
cov_rd_BCH_robust <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH",rd="robust")
cov_rd_LV_robust  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="robust")

## Compute estimate with model selection using RDHonest
cov_rd_CV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV" ,rd="honest",C=30)
cov_rd_BCH_honest <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH",rd="honest",C=30)
cov_rd_LV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="honest",C=30)
```
The selected covariates can be obtained via, e.g., `cov_rd_LV_honest$covs`. The estimated value for the RD parameter is contained in the `$rd` part of, e.g., `cov_rd_CV_honest`. The formatting of the this output is the same as in `rdhonest` or `rdrobust` (depending on the choice of `rd` in the call of `HighDim_rd`). Please see the references of the corresponding package.

In the above example the covariates are very sparse. In the following we reproduce a single run from the example of Section 4.2 of the paper where the covariate influence is actually not exactly sparse, but the covariates are still high-dimensional. The data generating process is as follows:
```
## Set seed
set.seed(12345)

## Set variables
n <- 1000 # Number of Obserbvations
p <- 200  # Number of Covariates

## Create Sigma
sigma_y <- 0.1295
sigma_z <- 0.1353
v <- 0.8*sqrt(6)*sigma_y^2/pi/(1:p)

Sigma <- sigma_z^2*diag(p+1)
Sigma[1,1] <- sigma_y^2
Sigma[1,2:(p+1)] <- v
Sigma[2:(p+1),1] <- v

## Running Variable
x <- 2*rbeta(n,2,4)-1

## Covariates and noise
randomness <- rmvnorm(n,rep(0,p+1),sigma=Sigma)
epsilon <- randomness[,1]
Z <- randomness[,2:(p+1)]

## Treatment Indicator
T <- as.numeric(x>=0)

## Outcome
alpha <- matrix(2/(1:p)^2,ncol=1)
Y <- as.numeric((1-T)*(0.36+0.96*x+5.47*x^2+15.28*x^3+15.87*x^4+5.14*x^5+0.22*Z%*%alpha)+
  T *(0.38+0.62*x-2.84*x^2+ 8.42*x^3-10.24*x^4+4.31*x^5+0.28*Z%*%alpha))+epsilon
```
In this case the covariates are already fairly high-dimensional. Therefore we do not add addititional Fourier expansions or interaction terms. The estimation proceeds as before:
```
## Compute estimate with model selection using rdrobust
cov_rd_CV_robust  <- HighDim_rd(Y,x,Z,tpc="CV" ,rd="robust")
cov_rd_BCH_robust <- HighDim_rd(Y,x,Z,tpc="BCH",rd="robust")
cov_rd_LV_robust  <- HighDim_rd(Y,x,Z,tpc="LV" ,rd="robust")

## Compute estimate with model selection using RDHonest
cov_rd_CV_honest  <- HighDim_rd(Y,x,Z,tpc="CV" ,rd="honest",C=30)
cov_rd_BCH_honest <- HighDim_rd(Y,x,Z,tpc="BCH",rd="honest",C=30)
cov_rd_LV_honest  <- HighDim_rd(Y,x,Z,tpc="LV" ,rd="honest",C=30)
```
The selected covariates and the estimates of their impact can be obtained as follows:
```
## Selected Covariates
cov_rd_CV_honest$covs
## Estimated Values
cov_rd_CV_honest$Zpars[cov_rd_CV_honest$covs]
```
