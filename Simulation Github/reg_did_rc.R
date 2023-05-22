reg_did_rc <-function(y, post, D, covariates, i.weights = NULL,
                      boot = FALSE, boot.type = "weighted", nboot = NULL,
                      inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # post as vector
  post <- as.vector(post)
  # Sample size
  n <- length(D)
  # outcome of interested
  y <- as.vector(y)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  
  
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                         subset = ((D==0) & (post==0)),
                                         weights = i.weights))
  if(anyNA(reg.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity of covariates is probably the reason for it.")
  }
  out.y.pre <-   as.vector(tcrossprod(reg.coeff.pre, int.cov))
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                          subset = ((D==0) & (post==1)),
                                          weights = i.weights))
  if(anyNA(reg.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is probably the reason for it.")
  }
  out.y.post <-   as.vector(tcrossprod(reg.coeff.post, int.cov))
  #-----------------------------------------------------------------------------
  #Compute the OR DID estimators
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont <- i.weights * D
  
  reg.att.treat.pre <- w.treat.pre * y
  reg.att.treat.post <- w.treat.post * y
  reg.att.cont <- w.cont * (out.y.post - out.y.pre)
  
  eta.treat.pre <- mean(reg.att.treat.pre) / mean(w.treat.pre)
  eta.treat.post <- mean(reg.att.treat.post) / mean(w.treat.post)
  eta.cont <- mean(reg.att.cont) / mean(w.cont)
  
  reg.att <- (eta.treat.post - eta.treat.pre) - eta.cont
  
  return(reg.att)
}