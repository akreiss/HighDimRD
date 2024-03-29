% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions.R
\name{HighDim_rd}
\alias{HighDim_rd}
\title{High-Dimensional Covariates in RDD Estimation}
\usage{
HighDim_rd(
  Y,
  X,
  Z,
  c = 0,
  rd = "robust",
  level = 0.95,
  b = NULL,
  bfactor = 1,
  h = NULL,
  tpc,
  kernel = "triangular",
  alpha = 0.05,
  M = NULL,
  L = 100,
  OPC = 50,
  C,
  sclass = "H",
  se.initial = "SilvermanNN"
)
}
\arguments{
\item{Y}{Vector of n responses.}

\item{X}{Vector of n running variables.}

\item{Z}{n x p - Matrix of covariates (each column corresponds to one
covariate).}

\item{c}{Cutoff point for treatment, individuals with X>=c are treated,
default is c=0.}

\item{rd}{Sets which RD functions should be used to find bandwidths,
estimators and confidence sets. If "honest", then the RDHonest package is
used. In this case the user has to specify C and should be aware of sclass.
If "robust", then the rdrobust package is used.}

\item{level}{Level of the computed confidence set, between 0 and 1. Default
is 0.95.}

\item{b}{Reference bandwidth used for doing the model selection step. If
b=NULL, the default, bandwidth selection from rdrobust without covariates
is used to find b. The final bandwidth used for model selection is given by
bfactor*b. b can be provided as the output of rdbwselect. It has to be a
list with at least the element bws which contains the bandwidth.}

\item{bfactor}{The bandwidth which is used for model selection is given by
bfactor*b (see above), the default is bfactor=1.}

\item{h}{Bandwidth for doing the actual RD estimation. If h=NULL, the
default, bandwidth selection from rdrobust with the selected covariate set
is used. See b for the format.}

\item{tpc}{Method for choosing the tuning the parameter of the Lasso.
Possible values are: "CV" for cross-validation. "LV" for the
bootstrap-procedure based on Lederer and Vogt (2020), this requires the
user specified parameters alpha, M and L. "BCH" for the procedure adapted
from Belloni et al. (2013) and "OPC" (observations per covariate) where the
number of selected covariates is no larger than n_effective/OPC, where
n_effective is the number of observations which receive a positive kernel
weight and OPC is specified by the user.}

\item{kernel}{Specifies which kernel to be used for both model selection and
estimation, possible choices are triangular (the default), epanechnikov,
uniform.}

\item{alpha, M, L}{Parameters for "LV" method for tuning parameter choice,
ignored if tpc!="LV".}

\item{OPC}{Parameter for "OPC" method for tuning parameter choice, ignored if
tpc!="OPC".}

\item{C, sclass, se.initial}{Parameters which are used for RDHonest: C is the
bound in the smoothness class specified in sclass. The default sclass is H,
see RDHonest for details. C must be specified by the user. See RDHonest for
a description of se.initial.}
}
\value{
List with five elements: \tabular{ll}{ \code{rd} \tab Output from
rdrobust or RDHonest (depending on the choice of rd) using the selected
covariates, see rdrobust or RDHonest for a description of the output, the
estimator can be found in rd$Estimate (for rdrobust) or in rd$estimate (for
RDHonest). \cr \code{covs} \tab Vector of selected covariate indices. \cr
\code{b} \tab Bandwidth used for model selection \cr \code{h} \tab
Bandwidth used for the estimation \cr \code{Z} \tab pars Vector of length p
containing the estimated parameters of the covariates.}
}
\description{
\code{HighDim_rd} function computes the treatment effect estimator in an RD
setting when a large list of covariates is supplied by the user.
}
\details{
For the exact functionality of \code{HighDim_rd}, we refer to Kreiss and
Rothe (2023). The vectors \code{Y},\code{X} contain the responses and the
running variable, respectively, that is, an individual is exposed to the
treatment if \code{X>=c}. The matrix \code{Z} contains the covariates (each
column contains one covariate). The lengths of \code{Y} and \code{X} must
equal the number of rows in \code{Z} which is the number of observations.
With \code{b} and \code{h} the bandwidths which are used for the model
selection and the actual estimation can be specified. Note that the bandwidth
which is used to do the bias correction in \code{rdrobust} will always be the
same as the same as the bandwidth which is used for the estimation. The
option \code{tpc} specifies which specific algorithm shall be used to find
the tuning parameter in the Lasso, that is, the model selection, step.
Details are provided in the parameter description above. The computationally
quickest method is OPC. \code{bfactor} allows to deliberately use smaller or
larger bandwidths for the model selection.
}
\section{Authors}{
 Alexander Kreiss, Leipzig University, Germany.
\email{alexander.kreiss@math.uni-leipzig.de}

Christoph Rothe, University of Mannheim, Germany
}

\section{References}{
 Belloni, A., Chernozhukov, V. and Hansen, C. (2013):
Inference on treatment effects after selection among high-dimensional
controls. The Review of Economic Studies, 81, 608-650

Kreiss, A. and Rothe, C. (2023): Inference in regression discontinuity
designs with high-dimensional covariates. The Econometrics Journal, Volume 26, 
Issue 2, May 2023, Pages 105–123

Lederer, J. and Vogt, M. (2020): Estimating the Lasso's Effective Noise.
arxiv 2004.11554
}

\examples{
\dontrun{library(mvtnorm)

set.seed(12345)

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
cov_rd_LV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="honest",C=30)}

}
\seealso{
\code{\link{fourier_basis}}, \code{\link{interaction_terms}},
\code{\link{cross_interactions}}
}
