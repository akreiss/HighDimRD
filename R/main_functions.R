#' High-Dimensional Covariates in RDD Estimation
#'
#'
#' \code{HighDim_rd} function computes the treatment effect estimator in an RD
#' setting when a large list of covariates is supplied by the user.
#'
#'
#' For the exact functionality of \code{HighDim_rd}, we refer to Kreiss and
#' Rothe (2021). The vectors \code{Y},\code{X} contain the responses and the
#' running variable, respectively, that is, an individual is exposed to the
#' treatment if \code{X>=c}. The matrix \code{Z} contains the covariates (each
#' column contains one covariate). The lengths of \code{Y} and \code{X} must
#' equal the number of rows in \code{Z} which is the number of observations.
#' With \code{b} and \code{h} the bandwidths which are used for the model
#' selection and the actual estimation can be specified. Note that the bandwidth
#' which is used to do the bias correction in \code{rdrobust} will always be the
#' same as the same as the bandwidth which is used for the estimation. The
#' option \code{tpc} specifies which specific algorithm shall be used to find
#' the tuning parameter in the Lasso, that is, the model selection, step.
#' Details are provided in the parameter description above. The computationally
#' quickest method is OPC. \code{bfactor} allows to deliberately use smaller or
#' larger bandwidths for the model selection.
#'
#' @param Y Vector of n responses.
#' @param X Vector of n running variables.
#' @param Z n x p - Matrix of covariates (each column corresponds to one
#'   covariate).
#' @param c Cutoff point for treatment, individuals with X>=c are treated,
#'   default is c=0.
#' @param rd Sets which RD functions should be used to find bandwidths,
#'   estimators and confidence sets. If "honest", then the RDHonest package is
#'   used. In this case the user has to specify C and should be aware of sclass.
#'   If "robust", then the rdrobust package is used.
#' @param niveau Niveau of the computed confidence set, between 0 and 1. Default
#'   is 0.95.
#' @param b Reference bandwidth used for doing the model selection step. If
#'   b=NULL, the default, bandwidth selection from rdrobust without covariates
#'   is used to find b. The final bandwidth used for model selection is given by
#'   bfactor*b. b can be provided as the output of rdbwselect. It has to be a
#'   list with at least the element bws which contains the bandwidth.
#' @param bfactor The bandwidth which is used for model selection is given by
#'   bfactor*b (see above), the default is bfactor=1.
#' @param h Bandwidth for doing the actual RD estimation. If h=NULL, the
#'   default, bandwidth selection from rdrobust with the selected covariate set
#'   is used. See b for the format.
#' @param tpc Method for choosing the tuning the parameter of the Lasso.
#'   Possible values are: "CV" for cross-validation. "LV" for the
#'   bootstrap-procedure based on Lederer and Vogt (2020), this requires the
#'   user specified parameters alpha, M and L. "BCH" for the procedure adapted
#'   from Belloni et al. (2013) and "OPC" (observations per covariate) where the
#'   number of selected covariates is no larger than n_effective/OPC, where
#'   n_effective is the number of observations which receive a positive kernel
#'   weight and OPC is specified by the user.
#' @param kernel Specifies which kernel to be used for both model selection and
#'   estimation, possible choices are triangular (the default), epanechnikov,
#'   uniform.
#' @param alpha,M,L Parameters for "LV" method for tuning parameter choice,
#'   ignored if tpc!="LV".
#' @param OPC Parameter for "OPC" method for tuning parameter choice, ignored if
#'   tpc!="OPC".
#' @param C,sclass Parameters which are used for RDHonest: C is the bound in the
#'   smoothness class specified in sclass. The default sclass is T, see RDHonest
#'   for details. C must be specified by the user.
#'
#' @return List with five elements: \tabular{ll}{ \code{rd} \tab Output from
#'   rdrobust or RDHonest (depending on the choice of rd) using the selected
#'   covariates, see rdrobust or RDHonest for a description of the output, the
#'   estimator can be found in rd$Estimate (for rdrobust) or in rd$estimate (for
#'   RDHonest). \cr \code{covs} \tab Vector of selected covariate indices. \cr
#'   \code{b} \tab Bandwidth used for model selection \cr \code{h} \tab
#'   Bandwidth used for the estimation \cr \code{Z} \tab pars Vector of length p
#'   containing the estimated parameters of the covariates.}
#'
#'
#' @section Authors: Alexander Kreiss, LSE, United Kingdom.
#'   \email{a.kreiss@@lse.ac.uk}
#'
#'   Christoph Rothe, University of Mannheim, Germany
#'
#' @section References: Belloni, A., Chernozhukov, V. and Hansen, C. (2013):
#'   Inference on treatment effects after selection among high-dimensional
#'   controls. The Review of Economic Studies, 81, 608-650
#'
#'   Kreiss, A. and Rothe, C. (2021): Inference in Regression Discontinuity
#'   Designs with High-Dimensional Covariates. arxiv 2110.13725
#'
#'   Lederer, J. and Vogt, M. (2020): Estimating the Lasso's Effective Noise.
#'   arxiv 2004.11554
#'
#' @examples \dontrun{library(mvtnorm)
#'
#' set.seed(12345)
#'
#' n <- 1000
#' sigma_y <- 0.1295
#' sigma_z <- 0.1353
#'
#' ## Running Variable
#' x <- 2*rbeta(n,2,4)-1
#'
#' ## Covariates
#' Z <- rmvnorm(n,rep(0,5),sigma=diag(sigma_z^2,5))
#'
#' ## Noise
#' epsilon <- rnorm(n,mean=0,sd=sigma_y)
#'
#' ## Treatment Indicator
#' T <- as.numeric(x>=0)
#'
#' ## Outcome
#' Y <- (1-T)*(0.36+0.96*x+5.47*x^2+15.28*x^3+15.87*x^4+5.14*x^5+0.22*Z[,1])+
#' T *(0.38+0.62*x-2.84*x^2+ 8.42*x^3-10.24*x^4+4.31*x^5+0.28*Z[,1])+epsilon
#'
#' ## High-Dimensional Covariates
#' Z1 <- fourier_basis(Z,4)
#' Z_HighDim <- cbind(Z,interaction_terms(Z),Z1,interaction_terms(Z1))
#'
#' ## Compute estimate with model selection using rdrobust
#' cov_rd_CV_robust  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV" ,rd="robust")
#' cov_rd_BCH_robust <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH",rd="robust")
#' cov_rd_LV_robust  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="robust")
#'
#' ## Compute estimate with model selection using RDHonest
#' cov_rd_CV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="CV" ,rd="honest",C=30)
#' cov_rd_BCH_honest <- HighDim_rd(Y,x,Z_HighDim,tpc="BCH",rd="honest",C=30)
#' cov_rd_LV_honest  <- HighDim_rd(Y,x,Z_HighDim,tpc="LV" ,rd="honest",C=30)}
#'
#' @seealso \code{\link{fourier_basis}}, \code{\link{interaction_terms}},
#'   \code{\link{cross_interactions}}
#' @export
HighDim_rd <- function(Y,X,Z,c=0,rd="robust",niveau=0.95,b=NULL,bfactor=1,h=NULL,tpc,kernel="triangular",alpha=0.05,M=NULL,L=100,OPC=50,C,sclass="T") {
  p <- dim(Z)[2]
  n <- length(Y)

  ## Set Treatment Status
  T <- as.numeric(X>=c)

  ## Chose the correct kernel function
  if(kernel=="epanechnikov") {
    kernel_fcn <- epanechnikov
  } else if(kernel=="triangular") {
    kernel_fcn <- triangular
  } else if(kernel=="uniform") {
    kernel_fcn <- uniform
  }

  ## Step 1: Find Bandwidth
  ## If not specified, find bandwidth in RD without covariates
  if(is.null(b)==TRUE) {
    if(rd=="robust") {
      bout <- rdrobust::rdbwselect(Y,X,c=c,bwselect="mserd",kernel=kernel)
      b <- bout$bws[1]
    } else {
      bout <- RDHonest::RDOptBW(Y~X,cutoff=c,M=C,kern=kernel,opt.criterion="FLCI",bw.equal=TRUE,alpha=1-niveau,sclass=sclass,order=1)
      b <- bout$h[1]
    }
  }

  ## Step 2: Model Selection
  ## Compute kernel weights
  kernel_factor <- kernel_fcn((X-c)/(bfactor*b))/(bfactor*b)
  relevant_indices <- which(kernel_factor>0)

  ## Transform Data
  Ymod  <- Y[relevant_indices]
  Tmod  <- T[relevant_indices]
  Xmod  <- X[relevant_indices]-c
  TXmod <- T[relevant_indices]*(X[relevant_indices]-c)
  Zmod  <- Z[relevant_indices,]

  ## Compute Centering and weights
  mu <- 1/n*colSums(Z[relevant_indices,]*matrix(rep(kernel_factor[relevant_indices],p),ncol=p))
  w <- bfactor*b/n*colSums((Z[relevant_indices,]*matrix(rep(kernel_factor[relevant_indices],p),ncol=p)-matrix(rep(mu,length(relevant_indices)),ncol=p,byrow=TRUE))^2)

  ## Regress Covariates on Outcome weighted by the kernel
  if(tpc=="CV") {
    ## Do Cross-Validation
    mod <- glmnet::cv.glmnet(cbind(Tmod,Xmod,TXmod,Zmod),Ymod,alpha=1,weights=kernel_factor[relevant_indices],penalty.factor=c(rep(0,3),sqrt(w)),standardize=FALSE)
    est <- coef(mod)
    Z_pars <- est[5:(4+p)] # Recall that there is the intercept
    theta_tilde <- est[1:4]
  } else if(tpc=="LV") {
    ## Bootstrap parameter choice according to Lederer and Vogt
    V <- cbind(1,Tmod,Xmod,TXmod)
    loadings <- sqrt(w)
    weights <- kernel_factor[relevant_indices]
    out <- bstpc(Ymod,V,Zmod,weights,loadings,M=M,L=L,alpha=alpha)
    Z_pars <- out$par
    theta_tilde <- out$theta_tilde
  } else if(tpc=="BCH") {
    ## Tuning Parameter chocice based on Chernozhukov et al.
    V <- cbind(1,Tmod,Xmod,TXmod)
    weights <- kernel_factor[relevant_indices]
    out <- BCHtpc(Ymod,V,Zmod,weights,n,bfactor*b)
    Z_pars <- out$par
    theta_tilde <- out$theta_tilde
  } else if(tpc=="OPC") {
    mod <-    glmnet::glmnet(cbind(Tmod,Xmod,TXmod,Zmod),Ymod,alpha=1,weights=kernel_factor[relevant_indices],penalty.factor=c(rep(0,3),sqrt(w)),standardize=FALSE)
    d <- max(which(mod$df<=4+length(relevant_indices)/OPC))
    est <- coef(mod,s=mod$lambda[d])
    Z_pars <- est[5:(4+p)]
    theta_tilde <- est[1:4]
  } else {
    stop("No tuning parameter choice method specified")
  }
  sig_cov <- which(abs(Z_pars)>0)

  ## Step 3: Do RD with chosen covariates
  ## Compute bandwidth if necessary
  if(is.null(h)==TRUE) {
    if(length(sig_cov)==0) {
      h <- b
    } else {
      if(rd=="robust") {
        hout <- rdrobust::rdbwselect(Y,X,covs=Z[,sig_cov],c=c,bwselect="mserd",kernel=kernel)
        h <- hout$bws[1]
      } else {
        ## Adjust for Covariate Influence
        cov_fit <- lm(Z[,sig_cov]~1+X+T+X*T,weights=kernel_factor)
        v <- cov_fit$residuals
        Sigma22 <- t(v)%*%(v*kernel_factor)/n
        Sigma21 <- colSums(as.matrix(v*(kernel_factor*Y)))/n
        Sigma <- solve(Sigma22)%*%Sigma21
        Ytilde <- Y-Z[,sig_cov]%*%Sigma

        hout <- RDHonest::RDOptBW(Ytilde~X,cutoff=c,M=C,kern=kernel,opt.criterion="FLCI",bw.equal=TRUE,alpha=1-niveau,sclass=sclass,order=1)
        h <- hout$h[1]
      }
    }
  }

  ## Estimation
  if(length(sig_cov)==0) {
    if(rd=="robust") {
      RDfit <- rdrobust::rdrobust(Y,X,c=c,h=h,b=h,kernel=kernel,level=niveau*100)
    } else {
      RDfit <- RDHonest::RDHonest(Y~X,cutoff=c,M=C,kern=kernel,opt.criterion="FLCI",h=h,alpha=1-niveau,sclass=sclass,order=1)
    }
  } else {
    if(rd=="robust") {
      RDfit <- rdrobust::rdrobust(Y,X,c=c,h=h,b=h,covs=Z[,sig_cov],kernel=kernel,level=niveau*100)
    } else {
      RDfit <- RDHonest::RDHonest(Ytilde~X,cutoff=c,M=C,kern=kernel,opt.criterion="FLCI",h=h,alpha=1-niveau,sclass=sclass,order=1)
    }
  }

  return(list(rd=RDfit,covs=sig_cov,b=b,h=h,Zpars=Z_pars))
}

#' Compute Interaction Terms
#'
#'
#' \code{interaction_terms} computes all interaction terms of the covariates,
#' that is, the products of columns, provided in z.
#'
#'
#' @param z n x p matrix, each column correpsonds to one covariate
#' @return Matrix with n rows and each columns, corresponds to one interaction,
#'   i.e., z\[,i\]*z\[,j\] for i,j=1,...,p and i<j. The column name gives the exact
#'   variables which were interacted.
#'
#'
#' @section Authors:
#' Alexander Kreiss, LSE, United Kingdom. \email{a.kreiss@@lse.ac.uk}
#'
#' Christoph Rothe, University of Mannheim, Germany
#'
#' @seealso \code{\link{fourier_basis}}, \code{\link{HighDim_rd}},
#'   \code{\link{cross_interactions}}
#'
#' @export
interaction_terms <- function(z) {
  p <- dim(z)[2]
  if(p==1) {stop("Input is one-dimensional, interactions make no sense.\n")}

  out <- matrix(NA,nrow=dim(z)[1],ncol=p*(p-1)/2)
  colnames(out) <- rep("a",p*(p-1)/2)
  col_count <- 1
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      out[,col_count] <- z[,i]*z[,j]
      colnames(out)[col_count] <- sprintf("IT %d*%d",i,j)
      col_count <- col_count+1
    }
  }
  return(out)
}

#' Compute Cross-Interactions
#'
#'
#' \code{cross_interactions} computes all cross interactions, that is, products
#' of columns, between the matrices z1 nad z2.
#'
#'
#' @param z1 n x p1 dimensional matrix, each column corresponds to one covariate.
#' @param z2 n x p2 dimensional maitrx, each column correpsonds to one covariate.
#'
#'
#' @return n x (p1*p2) Matrix which contains the cross interactions
#'   z1\[,i\]*z2\[,j\] for i=1,...,p1 and j=1,...,p2. The column name gives involved
#'   covariates.
#'
#' @section Authors:
#' Alexander Kreiss, LSE, United Kingdom. \email{a.kreiss@@lse.ac.uk}
#'
#' Christoph Rothe, University of Mannheim, Germany
#'
#'
#' @seealso \code{\link{fourier_basis}}, \code{\link{interaction_terms}},
#'   \code{\link{HighDim_rd}}
#'
#' @export
cross_interactions <- function(z1,z2,ident="") {
  p1 <- dim(z1)[2]
  p2 <- dim(z2)[2]
  n <- dim(z1)[1]

  out <- matrix(NA,nrow=n,ncol=p1*p2)
  colnames(out) <- rep("a",p1*p2)
  col_count <- 1
  for(i in 1:p1) {
    for(j in 1:p2) {
      out[,col_count] <- z1[,i]*z2[,j]
      colnames(out)[col_count] <- sprintf("CI %s: %d * %d",ident,i,j)
      col_count <- col_count+1
    }
  }

  return(out)
}


#' Compute Fourier Basis Expansions
#'
#'
#' \code{fourier_basis} computes Fourier Bases expansions for the covariates in z up to a given order.
#'
#'
#' @param z n x p Matrix, each column corresponds to one covariate.
#' @param order The order until which the Fourier expansion shall be computed.
#'
#'
#' @return n x (2 x p x order) Matrix, each column contains an expression of the
#'   form sin(2*pi*r*z[,i]) or cos(2*pi*r*z[,i]) for r=1,...,order and
#'   j=1,...,p. The column names are of the form "FB j sin r" or "FB j cos r".
#'
#' @section Authors:
#' Alexander Kreiss, LSE, United Kingdom. \email{a.kreiss@@lse.ac.uk}
#'
#' Christoph Rothe, University of Mannheim, Germany
#'
#'
#' @seealso \code{\link{HighDim_rd}}, \code{\link{interaction_terms}},
#'   \code{\link{cross_interactions}}
#'
#'
#' @export
fourier_basis <- function(z,order) {
  p <- dim(z)[2]

  out <- matrix(NA,nrow=dim(z)[1],ncol=p*2*order)

  colnames(out) <- rep("a",p*2*order)
  col_count <- 1
  for(i in 1:p) {
    for(j in 1:order) {
      out[,col_count] <- sin(2*pi*j*z[,i])
      colnames(out)[col_count] <- sprintf("FB %d sin %d",i,j)
      col_count <- col_count+1
      out[,col_count] <- cos(2*pi*j*z[,i])
      colnames(out)[col_count] <- sprintf("FB %d cos %d",i,j)
      col_count <- col_count+1
    }
  }
  return(out)
}
