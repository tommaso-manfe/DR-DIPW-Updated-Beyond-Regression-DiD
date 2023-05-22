drdid_imp_rc_hayek <-function(y, post, D, covariates, i.weights = NULL,
                              boot = FALSE, boot.type =  "weighted", nboot = NULL,
                              inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
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
  #Compute the Pscore using the pscore.cal
  pscore.ipt <- pscore.cal(D, int.cov, i.weights = i.weights, n = n)
  ps.fit <- as.vector(pscore.ipt$pscore)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  
  #Compute the Pscore using the pscore.cal
  int.cov.tscore=cbind(D, int.cov, int.cov[,2:ncol(int.cov)]*D)
  tscore.ipt <- pscore.cal(post, int.cov.tscore, i.weights = i.weights, n = n)
  ts.fit <- as.vector(tscore.ipt$pscore)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  hy_weights=i.weights * ps.fit * (1 - D) * (1 - post)/((1 - ps.fit)*(1-ts.fit))
  reg.cont.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                              subset = ((D==0) & (post==0)),
                                              weights = hy_weights))
  if(anyNA(reg.cont.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.pre <-   as.vector(tcrossprod(reg.cont.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the post-treatment period, using ols.
  hy_weights=i.weights * ps.fit * (1 - D) * (post)/((1 - ps.fit)*(ts.fit))
  reg.cont.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==0) & (post==1)),
                                               weights = hy_weights))
  if(anyNA(reg.cont.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.post <-   as.vector(tcrossprod(reg.cont.coeff.post, int.cov))
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  hy_weights=i.weights * D * (1 - post)/((1-ts.fit))
  reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==1) & (post==0)),
                                               weights = hy_weights))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  hy_weights=i.weights * D * post/(ts.fit)
  reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                subset = ((D==1) & (post==1)),
                                                weights = hy_weights))
  out.y.treat.post <-   as.vector(tcrossprod(reg.treat.coeff.post, int.cov))
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  # w.treat.pre <- i.weights * D * (1 - post)
  # w.treat.post <- i.weights * D * post
  # w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  # w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.treat.pre <- i.weights * D * (1 - post)/((1-ts.fit))
  w.treat.post <- i.weights * D * post/(ts.fit)
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/((1 - ps.fit)*(1-ts.fit))
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/((1 - ps.fit)*(ts.fit))
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  # eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  # eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  # eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  # eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- (w.treat.post/mean(w.treat.post)) * (out.y.treat.post - out.y.cont.post)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- (w.treat.pre/mean(w.treat.pre)) * (out.y.treat.pre - out.y.cont.pre)
  
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  
  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)
  
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
}
