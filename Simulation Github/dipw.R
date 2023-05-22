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