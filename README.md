# HighDimRD
This repository provides R-Code for easy inclusion of high-dimensional covariates in regression discontinuity design problems.

The package can be installed like this:
<code>
library(devtools)
install_github("akreiss/HighDimRD")
library(HighDimRD)
</code>

Here is an Example for how the package can be used:
<code>
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

## Compute estimate with model selection
cov_rd_CV  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV")
cov_rd_BCH <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH")
cov_rd_LV  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV")
</code>