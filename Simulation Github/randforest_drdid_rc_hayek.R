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
  ipw_C0=ps.fit_C0/(1-ps.fit_C0)
  dataor=as.data.frame(cbind(y_C0, int.cov_C0))
  model <- randomForest(y_C0 ~ ., data=dataor)
  out.y.cont.pre=predict(model, data, type="response")
  
  #Compute the Outcome regression for the control group at the post-treatment period, using RANDOM FOREST.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  ps.fit_C1=ps.fit[(D==0) & (post==1)]
  ipw_C1=ps.fit_C1/(1-ps.fit_C1)
  dataor=as.data.frame(cbind(y_C1, int.cov_C1))
  model <- randomForest(y_C1 ~ ., data=dataor)
  out.y.cont.post=predict(model, data, type="response")
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using RANDOM FOREST.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  dataor=as.data.frame(cbind(y_T0, int.cov_T0))
  model <- randomForest(y_T0 ~ ., data=dataor)
  out.y.treat.pre <- predict(model, data, type="response")
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using RANDOM FOREST.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  dataor=as.data.frame(cbind(y_T1, int.cov_T1))
  model <- randomForest(y_T1 ~ ., data=dataor)
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