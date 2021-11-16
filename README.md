# HighDimRD
This repository contains the R-Package HighDimRD for easy inclusion of high-dimensional covariates in regression discontinuity design problems. The main function for this purpose is HighDim_rd. Details of its functionality are contained in the paper "Inference in Regression Discontinuity Designs with High-Dimensional Covariates" which I have written with Christoph Rothe. A pre-print is available on <a href="https://arxiv.org/abs/2110.13725">arxiv</a>. The basic functionality is that in a first step the algorithm selects relevant covariates by means of a Lasso procedure. In a second step the treatment effect is estimated by using these selected covariates. In addition, there are some functions which create interactions, cross-interactions and Fourier expansions in order to generate high-dimensional covariates.

For the package to work, the packages glmnet, Matrix, rdrobust and RDHonest are required. The package and all its dependencies can be installed like this:
```
library(remotes)
install_github("kolesarm/RDHonest")
library(RDHonest)
install_github("akreiss/HighDimRD")
library(HighDimRD)

```

Here is an Example for how the package can be used:
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

## High-Dimensional Covariates
Z1 <- fourier_basis(Z,4)
Z_HighDim <- cbind(Z,interaction_terms(Z),Z1,interaction_terms(Z1))
## Compute estimate with model selection using rdrobust
cov_rd_CV_robust  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV" ,rd="robust")
cov_rd_BCH_robust <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH",rd="robust")
cov_rd_LV_robust  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="robust")

## Compute estimate with model selection using RDHonest
cov_rd_CV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV" ,rd="honest",C=30)
cov_rd_BCH_honest <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH",rd="honest",C=30)
cov_rd_LV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="honest",C=30)
```
