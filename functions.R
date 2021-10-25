## Required Libraries
library(glmnet)
library(rdrobust)


## This function computes the treatment effect estimator in an RD setting when
## a large list of covariates is supplied by the user.
## Input Variables:
## Y       - Vector of n responses
## X       - Vector of n running variables
## Z       - n x p - Matrix of covariates (each column corresponds to one
##           covariate)
## c       - Cutoff point for treatment, individuals with X>=c are treated,
##           default is c=0
## b       - Reference bandwidth used for doing the model selection step. If
##           b=NULL, the default, bandwidth selection from rdrobust without
##           covariates is used to find b. The final bandwidth used for model
##           selection is given by bfactor*b. b can be provided as the output
##           of rdbwselect. It has to be a list with at least the element bws
##           which contains the bandwidth.
## bfactor - The bandwidth which is used for model selection is given by
##           bfactor*b (see above), the default is bfactor=1
## h       - Bandwidth for doing the actual RD estimation. If h=NULL, the
##           default, bandwidth selection from rdrobust with the selected
##           covariate set is used. See b for the format.
## tpc     - Method for choosing the tuning the parameter of the Lasso. Possible
##           values are: "CV" for cross-validation. "LV" for the bootstrap-
##           procedure based on Lederer and Vogt (2020), this requires the user
##           specified parameters alpha, M and L. "BCH" for the procedure
##           adapted from Belloni et al. (2013) and "OPC" (observations per
##           covariate) where the number of selected covariates is no larger
##           than n_effective/OPC, where n_effective is the number of
##           observations which receive a positive kernel weight and OPC is
##           specified by the user.
## kernel  - Specifies which kernel to be used for both model selection and
##           estimation, possible choices are triangular (the default),
##           epanechnikov, uniform.
## alpha, 
## M, L    - Parameters for "LV" method for tuning parameter choice, ignored if
##           tpc!="LV"
## OPC     - Parameter for "OPC" method for tuning parameter choice, ignored if
##           tpc!="OPC"
## Output: List with four elements
## rd    - Output from rdrobust using the selected covariates, see rdrobust for
##         a description of the output, the estimator can be found in rd$Estimate
## covs  - Vector of selected covariate indices
## b     - Bandwidth used for model selection (as returned from rdbwselect
##         without covariates if b=NULL)
## h     - Bandwidth used for the estimation (as returned from rdbwselect with
##         the selected covariates if h=NULL)
## Zpars - Vector of length p containing the estimated parameters of the
##         covariates.
HighDim_rd <- function(Y,X,Z,c=0,b=NULL,bfactor=1,h=NULL,tpc,kernel="triangular",alpha=0.05,M=NULL,L=100,OPC=50) {
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
    kernel_fcn <- triangular
  }
  
  ## Step 1: Find Bandwidth
  ## If not specified, find bandwidth in RD without covariates
  if(is.null(b)==TRUE) {
    b <- rdbwselect(Y,X,c=c,bwselect="mserd",kernel=kernel)
  }
  
  ## Step 2: Model Selection
  ## Compute kernel weights
  kernel_factor <- kernel_fcn((X-c)/(bfactor*b$bws[1]))/(bfactor*b$bws[1])
  relevant_indices <- which(kernel_factor>0)
  
  ## Transform Data
  T     <- X>=c
  Ymod  <- Y[relevant_indices]
  Tmod  <- T[relevant_indices]
  Xmod  <- X[relevant_indices]-c
  TXmod <- T[relevant_indices]*(X[relevant_indices]-c)
  Zmod  <- Z[relevant_indices,]
  
  ## Compute Centering and weights
  mu <- 1/n*colSums(Z[relevant_indices,]*matrix(rep(kernel_factor[relevant_indices],p),ncol=p))
  w <- bfactor*b$bws[1]/n*colSums((Z[relevant_indices,]*matrix(rep(kernel_factor[relevant_indices],p),ncol=p)-matrix(rep(mu,length(relevant_indices)),ncol=p,byrow=TRUE))^2)
  
  ## Regress Covariates on Outcome weighted by the kernel
  if(tpc=="CV") {
    ## Do Cross-Validation
    mod <- cv.glmnet(cbind(Tmod,Xmod,TXmod,Zmod),Ymod,alpha=1,weights=kernel_factor[relevant_indices],penalty.factor=c(rep(0,3),sqrt(w)),standardize=FALSE)
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
    out <- BCHtpc(Ymod,V,Zmod,weights,n,bfactor*b$bws[1])
    Z_pars <- out$par
    theta_tilde <- out$theta_tilde
  } else if(tpc=="OPC") {
    mod <-    glmnet(cbind(Tmod,Xmod,TXmod,Zmod),Ymod,alpha=1,weights=kernel_factor[relevant_indices],penalty.factor=c(rep(0,3),sqrt(w)),standardize=FALSE)
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
      h <- rdbwselect(Y,X,covs=Z[,sig_cov],c=c,bwselect="mserd",kernel=kernel)
    }
  }
  
  ## Estimation
  if(length(sig_cov)==0) {
    RDfit <- rdrobust(Y,X,c=c,h=h$bws[1],b=h$bws[1],kernel=kernel)
  } else {
    RDfit <- rdrobust(Y,X,c=c,h=h$bws[1],b=h$bws[1],covs=Z[,sig_cov],kernel=kernel)
  }
  
  ## Compute estimate for variance
#  V <- cbind(1,T,X-c,T*(X-c))
#  r_hat <- Y-V%*%theta_tilde-Z%*%Z_pars
#  f_hat <- density(X,bw="nrd",n=1,from=c,to=c)$y
#  S_sq_hat <- sum(1/f_hat^2*1/(n*h)*kernel_fcn((X-c)/h)^2*(V%*%kappa_inv[2,])^2*r_hat^2)
#  sd_hat <- sqrt(S_sq_hat/n/h)
  
  return(list(rd=RDfit,covs=sig_cov,b=b,h=h,Zpars=Z_pars))
}

## Creates all interaction terms of the covariates provided in z.
## Input:
## z - n x p matrix, each column correpsonds to one covariate
## Output: Matrix with n rows and each columns, corresponds to one interaction,
##         i.e., z[,i]*z[,j] for i,j=1,...,p and i<j. The column name gives the
##         exact variables which were interacted.
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

## Creates all cross interactions between the matrices z1 nad z2.
## Input:
## z1 - n x p1 dimensional matrix, each column corresponds to one covariate.
## z2 - n x p2 dimensional maitrx, each column correpsonds to one covariate.
## Output: n x (p1*p2) Matrix which contains the cross interactions
##         z1[,i]*z2[,j] for i=1,...,p1 and j=1,...,p2. The column name gives
##         the involved covariates.
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


## Computes Fourier Bases expansions for the covariates in z up to a given order.
## Input:
## z     - n x p Matrix, each column corresponds to one covariate.
## order - The order until which the Fourier expansion shall be computed.
## Output: n x (2*p*order) Matrix, each column contains an expression of the form
##         sin(2*pi*r*z[,i]) or cos(2*pi*r*z[,i]) for r=1,...,order and j=1,...,p.
##         The column names are of the form "FB j sin r" or "FB j cos r".
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






















################################################################################
## Internal functions which are normally not called by the user. ###############
################################################################################

## Kernel Functions
epanechnikov <- function(x) {
  return(3/4*(1-x^2)*(abs(x)<=1))
}
triangular <- function(x) {
  return((1-abs(x))*(abs(x)<=1))
}
uniform <- function(x) {
  return(0.5*(abs(x)<=1))
}



## Find LASSO tuning parameter via the method from Vogt and Lederer
bstpc <- function(res,V,covs,weights,loadings,M=NULL,L=100,alpha=0.05) {
  #### Compute the LASSO path. To this end a re-scaling is required (or convenient)
  N <- length(res)
  p <- dim(covs)[2]
  
  ## Re-scale the weights
  s <- sum(weights)/N
  weights_tilde <- weights/s
  
  ## Remove the non-penalized parts
  Vbar  <- 1/N*t(V)%*%diag(weights)%*%V
  VZbar <- 1/N*t(V)%*%diag(weights)%*%covs
  alpha_zero <- solve(Vbar)%*%(1/N*t(V)%*%(weights*res))
  
   res_tilde <- res -V%*%alpha_zero
  covs_tilde <- covs-V%*%solve(Vbar)%*%VZbar
  
  ## Centralize
  mu_res <- sum(weights*res_tilde)/N
  mu_covs <- matrix(rep(colSums(covs_tilde*matrix(rep(weights,p),nrow=N,ncol=p,byrow=FALSE))/N,N),ncol=p,nrow=N,byrow=TRUE)
  
  ## Standardize
  sd_res     <- sqrt(sum(weights_tilde*(res_tilde-1/N*sum(weights_tilde*res_tilde))^2)/N)
  Sigma_covs <- diag(sqrt(colSums((covs_tilde-mu_covs/s)^2*matrix(rep(weights_tilde,p),nrow=N,ncol=p,byrow=FALSE))/N))
  
  mres  <- (res_tilde -mu_res/s )/sd_res
  mcovs <- (covs_tilde-mu_covs/s)%*%solve(Sigma_covs)
  
  ## Scale the loadings
  loadings_tilde <- loadings/diag(Sigma_covs)
  sl <- 1/p*sum(loadings_tilde)
  loadings_bar <- loadings_tilde/sl
  
  ## Compute LASSO path with re-scaled data
  LP <- glmnet(mcovs,mres,weights=weights_tilde,penalty.factor=loadings_bar,intercept=FALSE,standardize=FALSE)
  
  ## Compute lambda sequence to be used
  if(is.null(M)==TRUE) {
    M <- 5*p
  }
  lambda_max <- 2/N*max(abs(t(covs_tilde-mu_covs/s)%*%(weights*(res_tilde-mu_res/s)))/loadings)
  lambda_seq <- seq(from=0,to=lambda_max,length.out=M)
  
  ## Compute Estimators for each element in the lambda sequence
  betas           <- Matrix(0,nrow=p+4,ncol=M)
  betas[5:(p+4),] <- sd_res*solve(Sigma_covs)%*%coef(LP,s=lambda_seq*sl/2/sd_res/s)[2:(p+1),] # without intercept
  betas[1:4,]     <- alpha_zero-solve(Vbar)%*%VZbar%*%betas[5:(p+4),]
  
  ## Compute residuals
  residuals <- matrix(rep(res,M),nrow=N,ncol=M)-cbind(V,covs)%*%betas
  
  ## Compute Bootstrap quantiles
  bs_erros <- rep(rnorm(L*N),M)
  r_tilde <- c(matrix(rep(t(residuals),L),ncol=M,byrow=TRUE))
  r_star <- matrix(r_tilde*bs_erros,ncol=N,byrow=TRUE)
  
  A <- 2*abs(t(t(r_star%*%(weights*covs))/loadings)/N)
  q <- matrix(apply(A,1,max),ncol=L,byrow=TRUE)
  Q <- apply(q,1,quantile,probs=1-alpha)
  
  ## Select tuning parameter
  h <- cumsum(as.numeric(Q>lambda_seq))
  m_hat <- min(which(h==max(h)))
  
  lambda_hat <- lambda_seq[m_hat]
  par <- betas[5:(p+4),m_hat]
  theta_tilde <- betas[1:4,m_hat]
  
  return(list(Q=Q,lambad_seq=lambda_seq,m_hat=m_hat,lambda_hat=lambda_hat,par=par,theta_tilde=theta_tilde))
}

## Find LASSO tuning parameter via plug in rule from Belloni Chernozhukov, Hansen (2014)
BCHtpc <- function(Y,V,covs,weights,n,b,gamma=0.05,c=1.1,nu=0.00001,K=10) {
  #### Extract Information
  N <- length(Y)
  p <- dim(covs)[2]
  
  ## Perform Initial least squares without covariates
  init_model <- lm(Y~1+V[,2]+V[,3]+V[,4],weights=weights)
  residuals <- init_model$residuals
  
  ## Compute Initial Loadings
  loadings <- sqrt(colSums(covs^2*residuals^2*weights^2*b)/n)
  
  ## Scale of the input
  sigmaY <- sd(Y)
  alpha_K <- mean(weights)
  
  ## Tuning Parameter
  lambda <- 2*c*sqrt(n*b)*qnorm(1-gamma/(2*p))
  
  ## Update the loadings
  continue_updating <- TRUE
  count <- 1
  while(continue_updating) {
    ## Scale of the loadings
    alpha_L <- 1/(p+3)*sum(loadings)
    
    ## Compute New Model Fit
    mf <- glmnet(cbind(V[,2:4],covs),Y/sigmaY,weights=weights/alpha_K,standardize=FALSE,penalty.factor=c(0,0,0,loadings/alpha_L))
    a <- coef(mf,s=lambda*alpha_L/2/N/b/sigmaY/alpha_K)
    
    ## New residuals
    residuals <- sigmaY*(Y/sigmaY-cbind(V,covs)%*%a)
    new_loadings <- sqrt(colSums(covs^2*as.numeric(residuals)^2*as.numeric(weights)^2*b)/n)*sqrt(n*b/(n*b-sum(a!=0)+4))
    
    ## Check if we have to continue
    if(count>=K | max(abs(new_loadings-loadings))<nu) {
      continue_updating <- FALSE
    }
    count <- count+1
    loadings <- new_loadings
  }
  
  return(list(theta_tilde=a[1:4]*sigmaY,par=a[5:(p+4)]*sigmaY))
}


