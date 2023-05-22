#EXPERIMENT 1C: NON-RANDOMIZED EXPERIMENT WITH X-SPECIFIC TRENDS
#              PROPENSITY SCORES CORRECTLY SPECIFIED, OUTCOME REGRESSION NOT CORRECTLY SPECIFIED



#############################################################
# CREATE FUNCTIONS FOR LATER USE
#############################################################
file.name="parallel"


pscore.cal <- function(D, int.cov, i.weights, n){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  pslogit <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial",
                                         weights = i.weights))
  
  init.gamma <- suppressWarnings(stats::coef(pslogit))
  
  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))
  
  flag <- ifelse(pscore.cal$converged, 0, 1)
  
  gamma.cal <- try(pscore.cal$argument)
  
  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {
    
    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)
    
    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }
  
  #Compute fitted pscore and weights for regression
  pscore.index <- tcrossprod(gamma.cal, int.cov)
  pscore <- as.numeric(stats::plogis(pscore.index))
  
  if(flag==1) {
    warning("trust algorithm did not converge when estimating propensity score. \n Used IPT algorithm a la Graham et al (2012)")
  }
  if(flag==2) {
    warning(" Used glm algorithm to estimate propensity score as trust and IPT method did not converge")
    
    if(pslogit$converged == FALSE){
      warning(" glm algorithm did not converge")
    }
    
  }
  
  if(anyNA(pscore)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  
  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))
  
}


loss.ps.cal <- function(gam, D, int.cov, iw){
  #gam: argument (pscore parameters)
  #D: Treatment/Group indicator
  #int.cov: covariate matrix; n x k (with intercept!)
  # iw: sampling weights (useful for weighted bootstrap, too)
  
  n <- dim(int.cov)[1]
  k <- dim(int.cov)[2]
  
  if (!any(is.na(gam))) {
    ps.ind <- as.vector(int.cov %*% gam)
    #exp.ps.ind <- exp(-ps.ind)
    exp.ps.ind <- exp(ps.ind)
    #val <- base::mean(ifelse(D, exp.ps.ind, ps.ind) * iw)
    val <- - base::mean(ifelse(D, ps.ind, -exp.ps.ind) * iw)
    
    #grad <- apply(ifelse(D, -exp.ps.ind, 1) * iw * int.cov, 2, mean)
    grad <- - apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean)
    
    #hess <- (t(int.cov) %*% (ifelse(D, exp.ps.ind, 0) * iw * int.cov))/n
    hess <- - (t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov))/n
    
  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}

wols_rc <- function(y, post, D, int.cov, pscore, i.weights, pre = NULL, treat = F){
  #-----------------------------------------------------------------------------
  # Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  or.weights <- as.vector(i.weights * pscore/(1 - pscore))
  
  if((pre == T) & (treat==T)) {
    subs <- (D==1)*(post==0)
  }
  if((pre == F) & (treat==T)) {
    subs <- (D==1)*(post==1)
  }
  if((pre == T) & (treat==F)) {
    subs <- (D==0)*(post==0)
  }
  if((pre == F) & (treat==F)) {
    subs <- (D==0)*(post==1)
  }
  #Run weighted OLS
  beta.wls <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                    subset = subs==1,
                                    weights = or.weights))
  if(anyNA(beta.wls)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }
  
  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.wls, int.cov))
  
  # return fitted values
  return(list(out.reg = out.delta))
  
}



dipw <- function(y, Time, Treatment, X, method="hayek"){
  # D as vector
  Treatment <- as.vector(Treatment)
  # Sample size
  n <- length(Treatment)
  # y as vector
  y <- as.vector(y)
  # post as vector
  Time <- as.vector(Time)
  # Add constant to covariate vector
  # poly(as.matrix(X), degree=2, raw=TRUE)
  int.cov <- as.matrix(X)
  
  mylogit <- glm(Treatment ~ int.cov, family = binomial(link = "logit"))
  treatscore=mylogit$fitted.values
  
  mylogit <- glm(Time ~ Treatment+int.cov+int.cov*Treatment, family = binomial(link = "logit"))
  timescore=mylogit$fitted.values
  
  if (method=='horvitz'){
    
    Y_11_est=rep(NA, n)
    Y_10_est=rep(NA, n)
    Y_01_est=rep(NA, n)
    Y_00_est=rep(NA, n)
    
    
    Y_11_est=y*Treatment*Time/(treatscore*timescore)
    Y_10_est=y*Treatment*(1-Time)/(treatscore*(1-timescore))
    Y_01_est=y*(1-Treatment)*Time/((1-treatscore)*timescore)
    Y_00_est=y*(1-Treatment)*(1-Time)/((1-treatscore)*(1-timescore))
    
    
    est=rep(NA, n)
    est[Treatment==1 & Time==1]=Y_11_est[Treatment==1 & Time==1]
    est[Treatment==1 & Time==0]=-Y_10_est[Treatment==1 & Time==0]
    est[Treatment==0 & Time==1]=-Y_01_est[Treatment==0 & Time==1]
    est[Treatment==0 & Time==0]=Y_00_est[Treatment==0 & Time==0]
    est=est*treatscore/mean(Treatment)
    
  } else if (method=='hayek'){
    
    w_11_est=rep(NA, n)
    w_10_est=rep(NA, n)
    w_01_est=rep(NA, n)
    w_00_est=rep(NA, n)
    
    w_11_est=Treatment*Time/(timescore)
    w_10_est=Treatment*(1-Time)/((1-timescore))
    w_01_est=(1-Treatment)*Time*treatscore/((1-treatscore)*timescore)
    w_00_est=(1-Treatment)*(1-Time)*treatscore/((1-treatscore)*(1-timescore))
    
    w_11_est=y*w_11_est/mean(w_11_est)
    w_10_est=y*w_10_est/mean(w_10_est)
    w_01_est=y*w_01_est/mean(w_01_est)
    w_00_est=y*w_00_est/mean(w_00_est)
    
    est=rep(NA, n)
    est[Treatment==1 & Time==1]=w_11_est[Treatment==1 & Time==1]
    est[Treatment==1 & Time==0]=-w_10_est[Treatment==1 & Time==0]
    est[Treatment==0 & Time==1]=-w_01_est[Treatment==0 & Time==1]
    est[Treatment==0 & Time==0]=w_00_est[Treatment==0 & Time==0]
    
  } else{
    stop('Please define a valid weighting scheme: "hayek" and "horvitz" are the 
         two valid alternative methods')
  }
  
  return(mean(est))
}


drdid_rc1_hayek <-function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  mylogit <- glm(post ~ -1 + D+ int.cov +int.cov*D, family = binomial(link = "logit"))
  
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  ts.fit=as.vector(mylogit$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                         subset = ((D==0) & (post==0)),
                                         weights = i.weights))
  if(anyNA(reg.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.pre <-   as.vector(tcrossprod(reg.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                          subset = ((D==0) & (post==1)),
                                          weights = i.weights))
  if(anyNA(reg.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.post <-   as.vector(tcrossprod(reg.coeff.post, int.cov))
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  
  
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)/((1-ts.fit))
  w.treat.post <- i.weights * D * post/(ts.fit)
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/((1 - ps.fit)*(1-ts.fit))
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/((1 - ps.fit)*(ts.fit))
  
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(dr.att)
}

drdid_rc1 <-function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                         subset = ((D==0) & (post==0)),
                                         weights = i.weights))
  if(anyNA(reg.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.pre <-   as.vector(tcrossprod(reg.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                          subset = ((D==0) & (post==1)),
                                          weights = i.weights))
  if(anyNA(reg.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.post <-   as.vector(tcrossprod(reg.coeff.post, int.cov))
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  
  # Asymptotic linear representation of OLS parameters in pre-period
  weights.ols.pre <- i.weights * (1 - D) * (1 - post)
  wols.x.pre <- weights.ols.pre * int.cov
  wols.eX.pre <- weights.ols.pre * (y - out.y.pre) * int.cov
  XpX.inv.pre <- solve(crossprod(wols.x.pre, int.cov)/n)
  asy.lin.rep.ols.pre <-  wols.eX.pre %*% XpX.inv.pre
  
  # Asymptotic linear representation of OLS parameters in post-period
  weights.ols.post <- i.weights * (1 - D) * post
  wols.x.post <- weights.ols.post * int.cov
  wols.eX.post <- weights.ols.post * (y - out.y.post) * int.cov
  XpX.inv.post <- solve(crossprod(wols.x.post, int.cov)/n)
  asy.lin.rep.ols.post <-  wols.eX.post %*% XpX.inv.post
  
  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- stats::vcov(pscore.tr) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  
  # Estimation effect from beta hat from post and pre-periods
  # Derivative matrix (k x 1 vector)
  M1.post <- - base::colMeans(w.treat.post * post * int.cov)/mean(w.treat.post)
  M1.pre <- - base::colMeans(w.treat.pre * (1 - post) * int.cov)/mean(w.treat.pre)
  
  # Now get the influence function related to the estimation effect related to beta's
  inf.treat.or.post <- asy.lin.rep.ols.post %*% M1.post
  inf.treat.or.pre <- asy.lin.rep.ols.pre %*% M1.pre
  inf.treat.or <- inf.treat.or.post + inf.treat.or.pre
  
  # Influence function for the treated component
  inf.treat <- inf.treat.post - inf.treat.pre + inf.treat.or
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect from nuisance parameters
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)
  
  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2.pre <- base::colMeans(w.cont.pre *(y - out.y - att.cont.pre) * int.cov)/mean(w.cont.pre)
  M2.post <- base::colMeans(w.cont.post *(y - out.y - att.cont.post) * int.cov)/mean(w.cont.post)
  # Now the influence function related to estimation effect of pscores
  inf.cont.ps <- asy.lin.rep.ps %*% (M2.post - M2.pre)
  
  # Estimation effect from beta hat from post and pre-periods
  # Derivative matrix (k x 1 vector)
  M3.post <- - base::colMeans(w.cont.post * post * int.cov) / mean(w.cont.post)
  M3.pre <- - base::colMeans(w.cont.pre * (1 - post) * int.cov) / mean(w.cont.pre)
  
  # Now get the influence function related to the estimation effect related to beta's
  inf.cont.or.post <- asy.lin.rep.ols.post %*% M3.post
  inf.cont.or.pre <- asy.lin.rep.ols.pre %*% M3.pre
  inf.cont.or <- inf.cont.or.post + inf.cont.or.pre
  
  # Influence function for the control component
  inf.cont <- inf.cont.post - inf.cont.pre + inf.cont.ps + inf.cont.or
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + 1.96 * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - 1.96 * se.dr.att
    #Create this null vector so we can export the bootstrap draws too.
    dr.boot <- NULL
  }
  
  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      dr.boot <- mboot.did(dr.att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.dr.att <- stats::IQR(dr.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower doundary of 95% CI
      lci <- dr.att - cv * se.dr.att
    } else {
      # do weighted bootstrap
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_rc1,
                               n = n, y = y, post = post,
                               D = D, int.cov = int.cov, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.dr.att <- stats::IQR((dr.boot - dr.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((dr.boot - dr.att)/se.dr.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower doundary of 95% CI
      lci <- dr.att - cv * se.dr.att
      
    }
  }
  
  
  if(inffunc == FALSE) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = FALSE,
    estMethod = "trad2",
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "dr"
  )
  ret <- (list(ATT = dr.att,
               se = se.dr.att,
               uci = uci,
               lci = lci,
               boots = dr.boot,
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"
  
  # return the list
  return(ret)
}

drdid_rc1_horvitz <-function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  mylogit <- glm(post ~ -1 + D+ int.cov +int.cov*D, family = binomial(link = "logit"))
  
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  ts.fit=as.vector(mylogit$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                         subset = ((D==0) & (post==0)),
                                         weights = i.weights))
  if(anyNA(reg.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.pre <-   as.vector(tcrossprod(reg.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                          subset = ((D==0) & (post==1)),
                                          weights = i.weights))
  if(anyNA(reg.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.post <-   as.vector(tcrossprod(reg.coeff.post, int.cov))
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  
  
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)/((1-ts.fit))
  w.treat.post <- i.weights * D * post/(ts.fit)
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/((1 - ps.fit)*(1-ts.fit))
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/((1 - ps.fit)*(ts.fit))
  
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(D)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(D)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(D)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(D)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(dr.att)
}

drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.cont.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                              subset = ((D==0) & (post==0)),
                                              weights = i.weights))
  if(anyNA(reg.cont.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.pre <-   as.vector(tcrossprod(reg.cont.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the post-treatment period, using ols.
  reg.cont.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==0) & (post==1)),
                                               weights = i.weights))
  if(anyNA(reg.cont.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.post <-   as.vector(tcrossprod(reg.cont.coeff.post, int.cov))
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==1) & (post==0)),
                                               weights = i.weights))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                subset = ((D==1) & (post==1)),
                                                weights = i.weights))
  out.y.treat.post <-   as.vector(tcrossprod(reg.treat.coeff.post, int.cov))
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
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
  
  return(dr.att)
}

drdid_rc_horvitz <-function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  mylogit <- glm(post ~ -1 + D+ int.cov +int.cov*D, family = binomial(link = "logit"))
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  ts.fit=as.vector(mylogit$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.cont.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                              subset = ((D==0) & (post==0)),
                                              weights = i.weights))
  if(anyNA(reg.cont.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.pre <-   as.vector(tcrossprod(reg.cont.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the post-treatment period, using ols.
  reg.cont.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==0) & (post==1)),
                                               weights = i.weights))
  if(anyNA(reg.cont.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.post <-   as.vector(tcrossprod(reg.cont.coeff.post, int.cov))
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==1) & (post==0)),
                                               weights = i.weights))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                subset = ((D==1) & (post==1)),
                                                weights = i.weights))
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
  
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(D)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(D)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(D)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(D)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
  
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

drdid_rc_hayek <-function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  mylogit <- glm(post ~ -1 + D+ int.cov +int.cov*D, family = binomial(link = "logit"))
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  ts.fit=as.vector(mylogit$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.cont.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                              subset = ((D==0) & (post==0)),
                                              weights = i.weights))
  if(anyNA(reg.cont.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.pre <-   as.vector(tcrossprod(reg.cont.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the post-treatment period, using ols.
  reg.cont.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==0) & (post==1)),
                                               weights = i.weights))
  if(anyNA(reg.cont.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.cont.post <-   as.vector(tcrossprod(reg.cont.coeff.post, int.cov))
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==1) & (post==0)),
                                               weights = i.weights))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                subset = ((D==1) & (post==1)),
                                                weights = i.weights))
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

drdid_imp_rc1 <- function(y, post, D, covariates, i.weights = NULL,
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
  #Compute the Outcome regression for the control group
  out.y.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.pre <-  as.vector(out.y.pre$out.reg)
  out.y.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.post <-  as.vector(out.y.post$out.reg)
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)
  
  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)
  
  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)

  return(dr.att)
  
}

drdid_imp_rc <- function(y, post, D, covariates, i.weights = NULL, boot = FALSE,
                         boot.type =  "weighted",  nboot = NULL, inffunc = FALSE){
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
  #Compute the Outcome regression for the control group
  out.y.cont.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.cont.pre <-  as.vector(out.y.cont.pre$out.reg)
  out.y.cont.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.cont.post <-  as.vector(out.y.cont.post$out.reg)
  # Combine the ORs
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==1) & (post==0)),
                                               weights = i.weights))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                subset = ((D==1) & (post==1)),
                                                weights = i.weights))
  out.y.treat.post <-   as.vector(tcrossprod(reg.treat.coeff.post, int.cov))
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)
  
  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)
  
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
  
  return(dr.att)
  
  
}

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


std_ipw_did_rc <-function(y, post, D, covariates, i.weights = NULL,
                          boot = FALSE, boot.type = "weighted", nboot = NULL,
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
  #Pscore estimation (logit) and also its fitted values
  PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights))
  ps.fit <- as.vector(PS$fitted.values)
  # Do not divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * y / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * y / mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * y / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * y / mean(w.cont.post)

  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)

  # ATT estimator
  ipw.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  return(ipw.att)
}

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

DMLDiD=function(Y,D,p,T){
  N=length(Y)
  B=100
  set.seed(10000+i)
  random=sample(1:1000,B)
  
  thetabar=c(0)
  for (l in 1:B){
    k=2
    samplesplit=function(k,N){
      c1=1:N
      smp_size <- floor((1/k) * length(c1))
      
      ## set the seed to make your partition reproducible
      set.seed(random[l])
      train_ind <- sample(seq_len(length(c1)), size = smp_size)
      
      k1 <- c1[train_ind]
      k2 <- c1[-train_ind]
      return(rbind(k1,k2))
    }
    K=samplesplit(k,N)
    
    thetaDML=c(0)
    
    for (q in 1:k){
      ##Trimming
      set.seed(i)
      CV=cv.glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1)
      fit=glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1,lambda=CV$lambda.1se)
      beta1hat=fit$beta
      beta1hat <- as.numeric(as.character(beta1hat))
      
      ghat=1/(1+exp(-p[K[q,],]%*%beta1hat))
      
      index1=K[q,][which(ghat<0.95 & ghat>0.05)]
      
      ##Estimation
      ghat=1/(1+exp(-p[index1,]%*%beta1hat))
      
      lambda=mean(T[-K[q,]])
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=p[-K[q,],]
      XX=XX[index,]
      
      
      
      set.seed(i)
      CV=cv.glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1)
      fit=glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1,lambda=CV$lambda.1se)
      beta2hat=fit$beta
      beta2hat <- as.numeric(as.character(beta2hat))
      
      ellhat2=p[index1,]%*%beta2hat
      
      s=((T[index1]-lambda)*Y[index1]-ellhat2)*(D[index1]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[index1])
      s=s[which(s<abs(min(s)))]
      
      thetaDML[q]=mean(s)
      
    }
    
    thetabar[l]=mean(thetaDML)
    
    
  }
  finaltheta=mean(thetabar)
  finaltheta
  return(c(finaltheta))
}

lasso_drdid_rc_hayek <-function(y, post, D, covariates, i.weights = NULL,
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
  # Add second order interactions
  int.cov <- poly(as.matrix(covariates), degree=3, raw=TRUE)
  
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by LASSO
  lasso.fit <- cv.glmnet(int.cov, D, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights)
  
  ps.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov, type='response')
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ps.fit <- pmax(ps.fit, 1e-16)
  
  int.cov.ts=cbind(int.cov, D, int.cov[,2:ncol(int.cov)]*D)
  lasso.fit <- cv.glmnet(int.cov.ts, post, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights)
  
  ts.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov.ts, type='response')
  # Avoid divide by zero
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  ts.fit <- pmax(ts.fit, 1e-16)
  
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ts.fit_C0=ts.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/((1-ps.fit_C0)*(1-ts.fit_C0))
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0*ipw_C0)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ts.fit_C1=ts.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/((1-ps.fit_C1)*ts.fit_C1)
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1*ipw_C1)
  out.y.cont.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using LASSO.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  ts.fit_T0=ts.fit[(D==1) & (post==0)]
  ipw_T0=1/(1-ts.fit_T0)
  reg.treat.coeff.pre <- cv.glmnet(int.cov_T0, y_T0, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_T0*ipw_T0)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
   ts.fit_T1=ts.fit[(D==1) & (post==1)]
  ipw_T1=1/(ts.fit_T1)
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1*ipw_T1)
  out.y.treat.post=predict(reg.treat.coeff.post, s=reg.treat.coeff.post$lambda.1se, newx=int.cov)
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)/((1-ts.fit))
  w.treat.post <- i.weights * D * post/(ts.fit)
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/((1 - ps.fit)*(1-ts.fit))
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/((1 - ps.fit)*(ts.fit))
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  
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
  
  
  return(dr.att)
}

randforest_drdid_rc_hayek <-function(y, post, D, covariates, i.weights = NULL,
                                boot = FALSE, boot.type =  "weighted", nboot = NULL,
                                inffunc = FALSE){
  #----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- covariates
  
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by RANDOM FOREST
  data=as.data.frame(cbind(post, D, int.cov, int.cov*D))
  colnames(data)[7:10]=c('x1D', 'x2D','x3D','x4D')
  model <- randomForest(post ~ ., data=data)
  ts.fit=predict(model, data, type="response")
  
  data=as.data.frame(cbind(D, int.cov))
  model <- randomForest(D ~ ., data=data)
  ps.fit=predict(model, data, type="response")
  
  
  
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ts.fit <- pmin(ts.fit, 1 - 1e-16)
  ps.fit <- pmax(ps.fit, 1e-16)
  ts.fit <- pmax(ts.fit, 1e-16)
  
  
  #Compute the Outcome regression for the control group at the pre-treatment period, using RANDOM FOREST.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  ps.fit_C0=ps.fit[(D==0) & (post==0)]
  ts.fit_C0=ts.fit[(D==0) & (post==0)]
  ipw_C0=ps.fit_C0/((1-ps.fit_C0)*(1-ts.fit_C0))
  dataor=as.data.frame(cbind(y_C0, int.cov_C0))
  model <- randomForest(y_C0 ~ ., data=dataor, weights=i.weights_C0*ipw_C0)
  out.y.cont.pre=predict(model, data, type="response")
  
  #Compute the Outcome regression for the control group at the post-treatment period, using RANDOM FOREST.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ts.fit_C1=ts.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/((1-ps.fit_C1)*ts.fit_C1)
  dataor=as.data.frame(cbind(y_C1, int.cov_C1))
  model <- randomForest(y_C1 ~ ., data=dataor, weights=i.weights_C1*ipw_C1)
  out.y.cont.post=predict(model, data, type="response")
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using RANDOM FOREST.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  ts.fit_T0=ts.fit[(D==1) & (post==0)]
  ipw_T0=1/(1-ts.fit_T0)
  dataor=as.data.frame(cbind(y_T0, int.cov_T0))
  model <- randomForest(y_T0 ~ ., data=dataor, weights=i.weights_T0*ipw_T0)
  out.y.treat.pre <- predict(model, data, type="response")
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using RANDOM FOREST.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  ts.fit_T1=ts.fit[(D==1) & (post==1)]
  ipw_T1=1/(ts.fit_T1)
  dataor=as.data.frame(cbind(y_T1, int.cov_T1))
  model <- randomForest(y_T1 ~ ., data=dataor, weights=i.weights_T1*ipw_T1)
  out.y.treat.post=predict(model, data, type="response")
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)/((1-ts.fit))
  w.treat.post <- i.weights * D * post/(ts.fit)
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/((1 - ps.fit)*(1-ts.fit))
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/((1 - ps.fit)*(ts.fit))
  
  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)
  
  
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
  
  
  return(dr.att)
}
#create rep.row function for later use

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


################################################################################
# START OF THE SIMULATION
################################################################################

## Parameters and seed

#set.seed(1)      # Seed
library(doParallel)

Nslots <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
print(sprintf("%d slots were allocated", Nslots))
cl <- parallel::makeCluster(Nslots)
doParallel::registerDoParallel(cl)
M=Nslots

library(glmnet)


# Results Storage 
#TWFE methods
twfe<- rep(0,M)
twfe_time<- rep(0,M)
twfe_corr<- rep(0,M)
twfe_corr_time<- rep(0,M)

#IPW methods
dipw_hy <- rep(0,M)
dipw_hy_time <- rep(0,M)
ipw<- rep(0,M)
ipw_time<- rep(0,M)
chang<- rep(0,M)
chang_time<- rep(0,M)

#OR Methods
ordid<- rep(0,M)
ordid_time<- rep(0,M)

#DRDID methods
drdd_rc = rep(0,M)
drdd_rc_time = rep(0,M)
drdd_rc_hy = rep(0,M)
drdd_rc_hy_time = rep(0,M)

drdd_imp_rc_hy= rep(0,M)
drdd_imp_rc_hy_time= rep(0,M)
drdd_imp_rc= rep(0,M)
drdd_imp_rc_time= rep(0,M)

#ML methods
lasso_drdd_rc_hy = rep(0,M)
lasso_drdd_rc_hy_time = rep(0,M)
rf_drdd_rc_hy = rep(0,M)
rf_drdd_rc_hy_time = rep(0,M)


ATT<- rep(0,M)


# START OF THE SIMULATION LOOP

ptime <- system.time({
res=foreach(i=1:10000, .combine = 'rbind', .packages = c("glmnet", "trust", "randomForest")) %dopar% {
  
  
  # Sample size
  n <- 1000
  # pscore index (strength of common support)
  Xsi.ps <- .75
  # Proportion in each period
  lambda <- 0.5
  # NUmber of bootstrapped draws
  nboot <- 199
  #-----------------------------------------------------------------------------
  # Mean and Std deviation of Z's without truncation
  mean.z1 <- exp(0.25/2)
  sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
  mean.z2 <- 10
  sd.z2 <- 0.54164
  mean.z3 <- 0.21887
  sd.z3 <-   0.04453
  mean.z4 <- 402
  sd.z4 <-  56.63891
  #-----------------------------------------------------------------------------
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)
  
  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2
  
  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4
  
  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  #-----------------------------------------------------------------------------
  # Gen treatment groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (- z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4))
  d  <- as.numeric(runif(n) <= pi)
  #-----------------------------------------------------------------------------
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
  
  # Create heterogenenous effects for the ATT, which is set approximately equal to zero
  index.unobs.het <- d * (index.lin)
  index.att <- 0
  
  #This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4*x1 + 13.7*(x2 + x3 + x4)
  #v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)
  
  #Gen realized outcome at time 0
  y00 <- index.lin + v + stats::rnorm(n)
  y10 <- index.lin + v + stats::rnorm(n)
  
  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y01 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend #this is for the trend based on X
  
  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend + #this is for the trend based on X
    index.att # This is the treatment effects
  

  
  # Generate "T"
  #post <- as.numeric(stats::runif(n) <= lambda)
  ti_nt<- 0.5
  ti_t <- 0.5
  ti=d*ti_t+(1-d)*ti_nt
  post  <- as.numeric(runif(n) <= ti)
  
  
  y=d*post*y11+(1-d)*post*y01+(1-d)*(1-post)*y00+(d)*(1-post)*y10
  #-----------------------------------------------------------------------------
  #Gen id
  id <- 1:n
  #-----------------------------------------------------------------------------
  # Put in a long data frame
  dta_long <- as.data.frame(cbind(id = id, y = y, post = post, d = d,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- dta_long[order(dta_long$id),]  

  ################################################################################

  #Standard TWFE
  start.time <- Sys.time()
  twfe_i <- lm(y ~ x1 + x2 + x3 + x4+post+d+post*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  twfe_time[i]=round(as.numeric(substr(tk,0,6)),5)
  twfe[i] <- twfe_i$coefficients["post:d"] 

  #Corrected TWFE
  start.time <- Sys.time()
  twfe_corr_i <- lm(y ~ post+d+post*d+
                 +x1 + x2 + x3 + x4+
                 +x1*d + x2*d + x3*d + x4*d+
                 +x1*post + x2*post + x3*post + x4*post+
                 +x1*post*d + x2*post*d + x3*post*d + x4*post*d, data=dta_long)
  end.time <-Sys.time()
  tk <- end.time - start.time
  twfe_corr_time[i]=round(as.numeric(substr(tk,0,6)),5)
  twfe_corr[i] <- twfe_corr_i$coefficients["post:d"]+
                +twfe_corr_i$coefficients["post:d:x1"]*mean(dta_long$x1[dta_long$d==1])+
                +twfe_corr_i$coefficients["post:d:x2"]*mean(dta_long$x2[dta_long$d==1])+
                +twfe_corr_i$coefficients["post:d:x3"]*mean(dta_long$x3[dta_long$d==1])+ 
                +twfe_corr_i$coefficients["post:d:x4"]*mean(dta_long$x4[dta_long$d==1])

 # DIPW Hayek
  
  y=dta_long$y
  Time=dta_long$post
  Treatment=dta_long$d
  X=cbind(dta_long$x1, dta_long$x2, dta_long$x3, dta_long$x4)
  
  start.time <- Sys.time()
  dipw_hy_i <- dipw(y=y,Time=post, Treatment=d, X=X, method='hayek')
  end.time <-Sys.time()
  tk <- end.time - start.time
  dipw_hy_time[i]=round(as.numeric(substr(tk,0,6)),5)
  dipw_hy[i] <- dipw_hy_i
  
  
  y=dta_long$y
  post=dta_long$post
  D=dta_long$d
  covariates=dta_long[,5:8]
  i.weights = rep(1, nrow(dta_long))
  boot = FALSE
  boot.type =  "weighted"
  nboot = NULL
  inffunc = FALSE

  # LOCALLY EFFICIENT SANT'ANNA 
  
  start.time <- Sys.time()
  drdd_rc_i=drdid_rc(y=y, post = post,  D = D,
                     covariates = covariates,
                     boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  drdd_rc_time[i]=round(as.numeric(substr(tk,0,6)),5)
  drdd_rc[i] <- drdd_rc_i
  
 
  # LOCALLY EFFICIENT DR-DIPW HAYEK
  start.time <- Sys.time()
  drdd_rc_hy_i=drdid_rc_hayek(y=y, post = post,  D = D,
                                covariates = covariates,
                                boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  drdd_rc_hy_time[i]=round(as.numeric(substr(tk,0,6)),5)
  drdd_rc_hy[i] <- drdd_rc_hy_i

  # IPW ABADIE
  start.time <- Sys.time()
  ipw_i=std_ipw_did_rc(y=y, post = post,  D = D,
                                covariates = covariates,
                                boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ipw_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ipw[i] <- ipw_i

  # OUTCOME REGRESSION HECKMAN
  start.time <- Sys.time()
  ordid_i=reg_did_rc(y=y, post = post,  D = D,
                                covariates = covariates,
                                boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  ordid_time[i]=round(as.numeric(substr(tk,0,6)),5)
  ordid[i] <- ordid_i

  # IMPROVED LC DR-DIPW
   start.time <- Sys.time()
  drdd_imp_rc_hy_i=drdid_imp_rc_hayek(y=y, post = post,  D = D,
                               covariates = covariates,
                               boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  drdd_imp_rc_hy_time[i]=round(as.numeric(substr(tk,0,6)),5)
  drdd_imp_rc_hy[i] <- drdd_imp_rc_hy_i
  
  
  # IMPROVED LC DRDID RC SANT'ANNA
  
  start.time <- Sys.time()
  drdd_imp_rc_i=drdid_imp_rc(y=y, post = post,  D = D,
                             covariates = covariates,
                             boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  drdd_imp_rc_time[i]=round(as.numeric(substr(tk,0,6)),5)
  drdd_imp_rc[i] <- drdd_imp_rc_i

  #LASSO DR_DIPW
  lasso_drdd_rc_hy_i=lasso_drdid_rc_hayek(y=y, post = post,  D = D,
                              covariates = covariates,
                              boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  lasso_drdd_rc_hy_time[i]=round(as.numeric(substr(tk,0,6)),5)
  lasso_drdd_rc_hy[i] <- lasso_drdd_rc_hy_i

  #RF DR_DIPW
  rf_drdd_rc_hy_i=randforest_drdid_rc_hayek(y=y, post = post,  D = D,
                                          covariates = covariates,
                                          boot = FALSE, nboot = NULL, inffunc = FALSE)
  end.time <-Sys.time()
  tk <- end.time - start.time
  rf_drdd_rc_hy_time[i]=round(as.numeric(substr(tk,0,6)),5)
  rf_drdd_rc_hy[i] <- rf_drdd_rc_hy_i

  #Chang 
  p=as.matrix(dta_long[,5:8])
  
  
  start.time <- Sys.time()
  chang_i <- DMLDiD(Y=y,D=d,p=p,T=post)
  end.time <-Sys.time()
  tk <- end.time - start.time
  chang_time[i]=round(as.numeric(substr(tk,0,6)),5)
  chang[i] <- chang_i


  data.frame(index=i, ATT=ATT[i], twfe=twfe[i], twfe_time=twfe_time[i],
  twfe_corr=twfe_corr[i], twfe_corr_time=twfe_corr_time[i], 
  ipw=ipw[i], ipw_time=ipw_time[i], ordid=ordid[i], ordid_time=ordid_time[i],
  dipw_hy=dipw_hy[i], dipw_hy_time=dipw_hy_time[i], drdd_rc=drdd_rc[i],
   drdd_rc_time=drdd_rc_time[i], drdd_rc_hy = drdd_rc_hy[i], drdd_rc_hy_time = drdd_rc_hy_time[i],
  chang=chang[i], chang_time=chang_time[i], drdd_imp_rc=drdd_imp_rc[i], 
  drdd_imp_rc_time=drdd_imp_rc_time[i], drdd_imp_rc_hy=drdd_imp_rc_hy[i], 
  drdd_imp_rc_hy_time=drdd_imp_rc_hy_time[i], lasso_drdd_rc_hy=lasso_drdd_rc_hy[i],
  lasso_drdd_rc_hy_time=lasso_drdd_rc_hy_time[i], rf_drdd_rc_hy=rf_drdd_rc_hy[i],
  rf_drdd_rc_hy_time=rf_drdd_rc_hy_time[i] )
}
})
ptime
print(res)
saveRDS(res, "EXP_1C_storage.rds")

parallel::stopCluster(cl)