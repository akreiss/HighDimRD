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
  LP <- glmnet::glmnet(mcovs,mres,weights=weights_tilde,penalty.factor=loadings_bar,intercept=FALSE,standardize=FALSE)
  
  ## Compute lambda sequence to be used
  if(is.null(M)==TRUE) {
    M <- 5*p
  }
  lambda_max <- 2/N*max(abs(t(covs_tilde-mu_covs/s)%*%(weights*(res_tilde-mu_res/s)))/loadings)
  lambda_seq <- seq(from=0,to=lambda_max,length.out=M)
  
  ## Compute Estimators for each element in the lambda sequence
  betas           <- Matrix::Matrix(0,nrow=p+4,ncol=M)
  betas[5:(p+4),] <- sd_res*solve(Sigma_covs)%*%coef(LP,s=lambda_seq*sl/2/sd_res/s)[2:(p+1),] # without intercept
  betas[1:4,]     <- alpha_zero-solve(Vbar)%*%VZbar%*%betas[5:(p+4),]
  
  ## Compute residuals
  residuals <- matrix(rep(res,M),nrow=N,ncol=M)-cbind(V,covs)%*%betas
  
  ## Compute Bootstrap quantiles
  bs_erros <- rep(rnorm(L*N),M)
  r_tilde <- c(matrix(rep(Matrix::t(residuals),L),ncol=M,byrow=TRUE))
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
    mf <- glmnet::glmnet(cbind(V[,2:4],covs),Y/sigmaY,weights=weights/alpha_K,standardize=FALSE,penalty.factor=c(0,0,0,loadings/alpha_L))
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