load("bribes.RData")

Working_data=as.data.frame(Working_data[complete.cases(Working_data), ])
save(Working_data, file = "bribes.RData")
Y=Working_data[,1]
D=Working_data[,2]
T=Working_data[,3]
X=as.matrix(Working_data[,5:ncol(Working_data)])

Working_data$lvalue_tonnage[Working_data$D==1 & Working_data$T==1]=Working_data$lvalue_tonnage[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$lvalue_tonnage[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$lvalue_tonnage[Working_data$D==1 & Working_data$T==0])
      +mean(Working_data$lvalue_tonnage[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$lvalue_tonnage[Working_data$D==0 & Working_data$T==0]))

      


Working_data$hc_group2[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group2[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group2[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group2[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group2[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group2[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group3[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group3[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group3[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group3[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group3[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group3[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group4[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group4[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group4[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group4[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group4[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group4[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group5[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group5[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group5[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group5[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group5[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group5[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group6[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group6[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group6[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group6[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group6[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group6[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group7[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group7[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group7[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group7[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group7[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group7[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group8[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group8[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group8[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group8[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group8[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group8[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group9[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group9[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group9[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group9[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group9[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group9[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group10[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group10[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group10[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group10[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group10[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group10[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group11[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group11[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group11[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group11[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group11[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group11[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group12[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group12[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group12[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group12[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group12[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group12[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group13[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group13[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group13[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group13[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group13[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group13[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group14[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group14[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group14[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group14[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group14[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group14[Working_data$D==0 & Working_data$T==0]))

Working_data$hc_group15[Working_data$D==1 & Working_data$T==1]=Working_data$hc_group15[Working_data$D==1 & Working_data$T==1]+
  -(mean(Working_data$hc_group15[Working_data$D==1 & Working_data$T==1])+
      -mean(Working_data$hc_group15[Working_data$D==1 & Working_data$T==0])
    +mean(Working_data$hc_group15[Working_data$D==0 & Working_data$T==1])+
      -mean(Working_data$hc_group15[Working_data$D==0 & Working_data$T==0]))

Y2=Working_data[,1]
D2=Working_data[,2]
T2=Working_data[,3]
X2=as.matrix(Working_data[,5:ncol(Working_data)])

i.weights = NULL
boot = FALSE
boot.type =  "weighted"
nboot = NULL
inffunc = FALSE


wboot_lasso_imp_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
                                    boot = FALSE, boot.type =  "weighted", nboot = NULL,
                                    inffunc = FALSE){
  #-----------------------------------------------------------------------------
  
  
  
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  
  #bootstrap weights
  v <- stats::rexp(n)
  i.weights <- as.vector(i.weights * v)
  
  # Add second order interactions
  int.cov <- poly(as.matrix(covariates[,1:8]), degree=2, raw=TRUE)
  int.cov=cbind(int.cov, covariates[,9:ncol(covariates)])
  
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
  ipw_C1=ps.fit_C1/((1-ps.fit_C1)*(ts.fit_C1))
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
  dr.att.b <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
  
  
  return(dr.att.b)
}

lasso_imp_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
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
  int.cov <- poly(as.matrix(covariates[,1:12]), degree=2, raw=TRUE)
  int.cov=cbind(int.cov, covariates[,13:ncol(covariates)])
  
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Compute the Pscore by LASSO
  lasso.fit <- cv.glmnet(int.cov, D, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights, nfolds=30)
  
  ps.fit=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=int.cov, type='response')
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  ps.fit <- pmax(ps.fit, 1e-16)
  
  int.cov.ts=cbind(int.cov, D, int.cov[,2:ncol(int.cov)]*D)
  lasso.fit <- cv.glmnet(int.cov.ts, post, type.measure="deviance", 
                         alpha=1, family="binomial", weights=i.weights, nfolds=30)
  
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
                                  alpha=1, family="gaussian", weights=i.weights_C0*ipw_C0, nfolds=30)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ts.fit_C1=ts.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/((1-ps.fit_C1)*(ts.fit_C1))
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1*ipw_C1, nfolds=30)
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
                                   alpha=1, family="gaussian", weights=i.weights_T0*ipw_T0, nfolds=30)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  ts.fit_T1=ts.fit[(D==1) & (post==1)]
  ipw_T1=1/(ts.fit_T1)
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1*ipw_T1, nfolds=30)
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
  
  #-----------------------------------------------------------------------------
  # Bootstrapped standard error
  #dr.boot <- unlist(lapply(1:nboot, wboot_lasso_imp_drdid_rc,
   #                          y = y, post = post,
    #                       D = D, covariates = covariates, i.weights = i.weights))
  
  # get bootstrap std errors based on IQR
  #se.dr.att <- stats::IQR((dr.boot - dr.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
  # get symmtric critival values
  #cv <- stats::quantile(abs((dr.boot - dr.att)/se.dr.att), probs = 0.95)
  # Estimate of upper boudary of 95% CI
  #uci <- dr.att + cv * se.dr.att
  # Estimate of lower doundary of 95% CI
  #lci <- dr.att - cv * se.dr.att
  
  
  return(dr.att)
}


library(doParallel)

Nslots <- as.numeric(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
print(sprintf("%d slots were allocated", Nslots))
cl <- parallel::makeCluster(Nslots)
doParallel::registerDoParallel(cl)
M=Nslots

res=foreach(i=1:1000, .combine = 'rbind', .packages = c("glmnet", "trust", "randomForest")) %dopar% {
    n=length(Y)
    
     i.weights <- stats::rexp(n)
     
lasso_hat=lasso_imp_drdid_rc(y=Y, post=T, D=D, covariates=X, i.weights = i.weights, boot = FALSE,
                           boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)

      data.frame(index=i, lasso_hat)
}


saveRDS(res, "bribes_lasso_hat_se.rds")

parallel::stopCluster(cl)