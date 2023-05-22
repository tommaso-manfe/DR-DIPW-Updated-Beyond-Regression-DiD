#APPLICATION: the effect of tariff reduction on corruption behaviors by using the
#bribe payment data collected by Sequeira (2016) between South Africa and Mozambique.
#The work is built on Chang (2020) who analyzed the same dataset
rm(list = ls())
#Set your WD
#setwd("C:/Users/tommy/OneDrive/Desktop/Tesi/Application/bribes/Replication_Final_original")

#Libraries
library(readstata13)
library(glmnet)
library(DRDID)
library(stargazer)
#Raw Data
data <- read.dta13 ("Bribes_Regression.dta")

#Prepare dataset regression
Y=data$lba
D=data$tariff_change_2008
T=data$post_2008

clear_agent2=(data$clear_agent==2)
clear_agent3=(data$clear_agent==3)
clear_agent4=(data$clear_agent==4)
clear_agent5=(data$clear_agent==5)
clear_agent6=(data$clear_agent==6)
clear_agent7=(data$clear_agent==7)
clear_agent8=(data$clear_agent==8)

hc_group2=(data$hc_group==2)
hc_group3=(data$hc_group==3)
hc_group4=(data$hc_group==4)
hc_group5=(data$hc_group==5)
hc_group6=(data$hc_group==6)
hc_group7=(data$hc_group==7)
hc_group8=(data$hc_group==8)
hc_group9=(data$hc_group==9)
hc_group10=(data$hc_group==10)
hc_group11=(data$hc_group==11)
hc_group12=(data$hc_group==12)
hc_group13=(data$hc_group==13)
hc_group14=(data$hc_group==14)
hc_group15=(data$hc_group==15)

tariff2007=data$tariff2007
lvalue_tonnage=data$lvalue_tonnage
differentiated=data$differentiated
agri=data$agri
perishable=data$perishable
dfs=data$dfs
day_w_arrival=data$day_w_arrival
monitor=data$monitor
psi=data$psi
hc_4digits=data$hc_4digits
term=data$term
rsa=data$rsa
X=cbind(tariff2007,lvalue_tonnage, differentiated, agri, perishable,dfs, day_w_arrival, monitor, psi, hc_4digits, term, rsa,clear_agent2,clear_agent3,clear_agent4,clear_agent5,clear_agent6,clear_agent7,clear_agent8,hc_group2,hc_group3,hc_group4,hc_group5,hc_group6,hc_group7,hc_group8,hc_group9,hc_group10,hc_group11,hc_group12,hc_group13,hc_group14,hc_group15)
#X=cbind(tariff2007,lvalue_tonnage, differentiated, agri, perishable, dfs, day_w_arrival, monitor, psi,rsa, data$hc_group, data$clear_agent, hc_4digits, term)

#Working data
Working_data=cbind(Y,D,T,D*T,X)

Working_data=as.data.frame(Working_data[complete.cases(Working_data), ])

Y=Working_data[,1]
D=Working_data[,2]
T=Working_data[,3]
X=as.matrix(Working_data[,5:ncol(Working_data)])



#Replication Equation 1 Table 9 Siqueira(2016)
twfe=miceadds::lm.cluster( data=Working_data, formula=Y ~ . ,cluster="hc_4digits")
twfe_coeff=summary(twfe)[,1]['V4']
twfe_se=summary(twfe)[,2]['V4']

#Replication Equation 2 Table 9 Siqueira(2016)
twfe_int=miceadds::lm.cluster( data=Working_data, formula=Y ~ . 
                               +differentiated*T+agri*T+lvalue_tonnage*T+perishable*T+day_w_arrival*T+dfs*T+psi*T+tariff2007*T
                              ,cluster="hc_4digits")
twfe_int_coeff=summary(twfe_int)[,1]['V4']
twfe_int_se=summary(twfe_int)[,2]['V4']

#Replication Equation 1 Table 9 Siqueira(2016)
twfe_int2=miceadds::lm.cluster( data=Working_data, formula=Y ~ .+differentiated*T+agri*T+lvalue_tonnage*T+perishable*T+day_w_arrival*T+dfs*T+psi*T+tariff2007*T
                                +differentiated*D+agri*D+lvalue_tonnage*D+perishable*D+day_w_arrival*D+dfs*D+psi*D+tariff2007*D,cluster="hc_4digits")
twfe_int2_coeff=summary(twfe_int2)[,1]['V4']
twfe_int2_se=summary(twfe_int2)[,2]['V4']

Working_data=Working_data[,-4]

#Test alternative methods


tripleIPWRA2 <- function(id,t, d, X, Y, method='logit'){
  
  # triple DiD propensity score weighting following Bludell(2004)
  # Y is the dependent variable
  # X is a group of covariates
  # X matrix of covariates
  # t is the time dummy
  # d is the treatment group dummy
  # T1 TO CO C1 are four dummies which take 1 if the individual is part of the category
  # T0 is treated group at time 0
  # T1 is treated group at time 1
  # C0 is control group at time 0
  # C1 is control group at time 1
  # method selects the estimation method for the propensity score. The available
  #estimation methods are:"logit", "probit", "lasso" and "randomforest".
  
  
  
  
  # create four groups and prepare dataset
  T1=ifelse(t==1& d==1,1,0)
  T0=ifelse(t==0& d==1,1,0)
  C1=ifelse(t==1& d==0,1,0)
  C0=ifelse(t==0& d==0,1,0)
  
  data=as.data.frame(cbind(Y,id, t, d, C0, C1,T0, T1, X))
  len=nrow(data)
  
  #First propensity score matching between T1 and C1
  data1=subset(data, T1==1 | C1==1)
  K=nrow(data1)
  X1=X[T1==1|C1==1,]
  
  #Estimate propensity score
  if (method=='logit'){
    mylogit <- glm(T1 ~ X1, family = binomial(link = "logit"), 
                   data = data1)
    data1$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X1, family = binomial(link = "probit"), 
                    data = data1)
    data1$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    
    cov1 <- poly(as.matrix(data1[,9:(ncol(data1)-21)]), degree=2, raw=TRUE)
    cov1=as.matrix(cbind(cov1,data1[,(ncol(data1)-21):(ncol(data1))]))
    lasso.fit <- cv.glmnet(cov1, data1$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data1$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov1, type='response')
    
  } else if (method=='randomforest'){
    data1rf=data1[,8:ncol(data1)] 
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X1))))
    fr.fit=cforest(T1~ .,data=data1rf, controls=mycontrols)
    data1$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data1$pscore[data1$pscore>1]=1
  data1$pscore[data1$pscore<0]=0
  
  
  #Second propensity score matching between T1 and C0
  
  data2=subset(data, T1==1 | C0==1)
  K=nrow(data2)
  X2=X[T1==1|C0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X2, family = binomial(link = "logit"), 
                   data = data2)
    data2$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X2, family = binomial(link = "probit"), 
                    data = data2)
    data2$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # LASSO linear regression
    # perform LASSO linear regression with 10-fold cross validation
    cov2 <- poly(as.matrix(data2[,9:(ncol(data2)-21)]), degree=2, raw=TRUE)
    cov2=as.matrix(cbind(cov2,data2[,(ncol(data1)-21):(ncol(data2))]))
    
    lasso.fit <- cv.glmnet(cov2, data2$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data2$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov2, type="response")
    
  } else if (method=='randomforest'){
    data2rf=data2[,8:ncol(data2)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X2))))
    fr.fit=cforest(T1~ .,data=data2rf, controls=mycontrols)
    data2$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  #Check that the p-value is in between 0 and 0
  data2$pscore[data2$pscore>1]=1
  data2$pscore[data2$pscore<0]=0
  
  
  #Third propensity score matching, between T1 and T0
  
  data3=subset(data, T1==1 | T0==1)
  K=nrow(data3)
  X3=X[T1==1|T0==1,]
  
  #Estimate propensity score
  
  if (method=='logit'){
    mylogit <- glm(T1 ~ X3, family = binomial(link = "logit"), 
                   data = data3)
    data3$pscore=mylogit$fitted.values
    
  } else if(method=='probit') {
    myprobit <- glm(T1 ~ X3, family = binomial(link = "probit"), 
                    data = data3)
    data3$pscore=myprobit$fitted.values
    
  } else if(method=='lasso') {
    
    # Estimate propensity score with LASSO linear regression 
    # perform LASSO linear regression with 10-fold cross validation
    cov3 <- poly(as.matrix(data3[,9:(ncol(data3)-21)]), degree=2, raw=TRUE)
    cov3=as.matrix(cbind(cov3,data3[,(ncol(data3)-21):(ncol(data3))]))
    
    lasso.fit <- cv.glmnet(cov3, data3$T1, type.measure="deviance", 
                           alpha=1, family="binomial")
    
    data3$pscore=predict(lasso.fit, s=lasso.fit$lambda.1se, newx=cov3, type="response")
    
  } else if (method=='randomforest'){
    
    #Estimate propensity score with random forest
    data3rf=data3[,8:ncol(data3)]
    mycontrols <- cforest_unbiased(ntree=100, mtry=as.integer(sqrt(ncol(X3))))
    fr.fit=cforest(T1~ .,data=data3rf, controls=mycontrols)
    data3$pscore=as.numeric(predict(fr.fit, type='prob'))
    
  } else{
    stop('Please define a valid estimation method for the propensity score: "logit", "probit", "lasso"
          and "randomforest" are the four valid alternatives method')
  }
  
  
  #Check that the p-value is in between 0 and 0
  data3$pscore[data3$pscore>1]=1
  data3$pscore[data3$pscore<0]=0
  
  
  #Merge propensity score
  data=rbind(data1, data2[data2$T1!=1,],data3[data3$T1!=1,])
  data$pscore[data$T1==1]=(data1$pscore[data1$T1==1]+data2$pscore[data2$T1==1]+
                             data3$pscore[data3$T1==1])/3
  
  #trimming for pscore value close to 0 and to 1
  cs_data=subset(data, pscore>=0.01 & pscore<=0.99)
  
  #Preparing data for regression
  cs_len=nrow(cs_data)
  Xmatrix=as.matrix(cs_data[9:(ncol(cs_data)-1)])
  
  #Computing weights
  cs_data$w_att=rep(1,cs_len)
  cs_data$w_att=ifelse(cs_data$T1==0, cs_data$pscore/(1-cs_data$pscore),1)
  
  
  #Regression with propensity score weights
  e_tripleIPWRA2=miceadds::lm.cluster( data=cs_data, formula=Y ~ t*d+t+d+Xmatrix,,weights=cs_data$w_att, cluster="hc_4digits")
  
  return(e_tripleIPWRA2)
  
}

id=1:nrow(Working_data)
lasso_3IPWRA=tripleIPWRA2(id=id, t=T, d=D, X=X, Y=Y, method='lasso')
lasso_3IPWRA_coeff=lasso_3IPWRA$lm_res$coefficients['t:d']
lasso_3IPWRA_se=summary(lasso_3IPWRA)[,2]['t:d']
rf_3IPWRA=tripleIPWRA2(id=id, t=T, d=D, X=X, Y=Y, method='randomforest')
rf_3IPWRA_coeff=rf_3IPWRA$lm_res$coefficients['t:d']
rf_3IPWRA_se=summary(rf_3IPWRA)[,2]['t:d']


Working_data$id=id
#drdid<- drdid(yname="Y", tname = "T", idname = "id",  dname = "D",
#                 xformla= ~ V5 + V6 + V7+V8 + V9 + V10 + V11,
#                 data = Working_data, panel = FALSE, estMethod = "imp")

# Bootrstapped "Improved" Doubly Robust Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

# Bootrstapped Improved Doubly Robust Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

wboot_drdid_imp_rc1 <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the Pscore using the pscore.cal
  ps.b <- pscore.cal(D, int.cov, i.weights = b.weights, n = n)
  ps.b <- as.vector(ps.b$pscore)
  ps.b <- pmin(ps.b, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.pre.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = T, treat = F)
  out.y.pre.b <-  as.vector(out.y.pre.b$out.reg)
  out.y.post.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = F, treat = F)
  out.y.post.b <-  as.vector(out.y.post.b$out.reg)
  
  # Combine the ORs
  out.y.b <- post * out.y.post.b + (1 - post) * out.y.pre.b
  
  # Compute AIPW estimator
  att.b <- aipw_did_rc1(y, post, D, ps.b, out.y.b, b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}

###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


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


# Loss function for estimation of the calibrated PS, using trust Based on Tan (2019).
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

# Bootrstapped Improved and locally efficient Doubly Robust Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

wboot_drdid_imp_rc <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the Pscore using the pscore.cal
  ps.b <- pscore.cal(D, int.cov, i.weights = b.weights, n = n)
  ps.b <- as.vector(ps.b$pscore)
  ps.b <- pmin(ps.b, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.cont.pre.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = T, treat = F)
  out.y.cont.pre.b <-  as.vector(out.y.cont.pre.b$out.reg)
  out.y.cont.post.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = F, treat = F)
  out.y.cont.post.b <-  as.vector(out.y.cont.post.b$out.reg)
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                 subset = ((D==1) & (post==0)),
                                                 weights = b.weights))
  out.y.treat.pre.b <-   as.vector(tcrossprod(reg.treat.coeff.pre.b, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                  subset = ((D==1) & (post==1)),
                                                  weights = b.weights))
  out.y.treat.post.b <-   as.vector(tcrossprod(reg.treat.coeff.post.b, int.cov))
  
  # Compute AIPW estimator
  att.b <- aipw_did_rc(y, post, D, ps.b,
                       out.y.treat.post.b, out.y.treat.pre.b,
                       out.y.cont.post.b, out.y.cont.pre.b,
                       b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}

lasso_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
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
  int.cov=as.matrix(cbind(int.cov,covariates[,13:(ncol(covariates))]))
  
  
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
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0)
  
  out.y.cont.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1)
  out.y.cont.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  
  
  #Compute the Outcome regression for the treated group at the pre-treatment period, using LASSO.
  y_T0=y[(D==1) & (post==0)]
  int.cov_T0=int.cov[(D==1) & (post==0),]
  i.weights_T0=i.weights[(D==1) & (post==0)]
  reg.treat.coeff.pre <- cv.glmnet(int.cov_T0, y_T0, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_T0)
  
  out.y.treat.pre=predict(reg.treat.coeff.pre, s=reg.treat.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the treated group at the post-treatment period, using LASSO.
  y_T1=y[(D==1) & (post==1)]
  int.cov_T1=int.cov[(D==1) & (post==1),]
  i.weights_T1=i.weights[(D==1) & (post==1)]
  reg.treat.coeff.post <- cv.glmnet(int.cov_T1, y_T1, type.measure="mse",
                                    alpha=1, family="gaussian", weights=i.weights_T1)
  out.y.treat.post=predict(reg.treat.coeff.post, s=reg.treat.coeff.post$lambda.1se, newx=int.cov)
  
  
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
  
  
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  # Influence function for the treated component
  inf.treat <- inf.treat.post - inf.treat.pre
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect from nuisance parameters
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)
  
  # Influence function for the control component
  inf.cont <- inf.cont.post - inf.cont.pre
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func1 <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  # Now, we only need to get the influence function of the adjustment terms
  # First, the terms as if all OR parameters were known
  inf.eff1 <- eta.d.post - w.d * att.d.post/mean(w.d)
  inf.eff2 <- eta.dt1.post - w.dt1 * att.dt1.post/mean(w.dt1)
  inf.eff3 <- eta.d.pre - w.d * att.d.pre/mean(w.d)
  inf.eff4 <- eta.dt0.pre - w.dt0 * att.dt0.pre/mean(w.dt0)
  inf.eff <- (inf.eff1 - inf.eff2) - (inf.eff3 - inf.eff4)
  #-----------------------------------------------------------------------------
  #get the influence function of the locally efficient DR estimator (put all pieces together)
  dr.att.inf.func <- dr.att.inf.func1 + inf.eff
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
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_imp_rc,
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
    estMethod = "imp",
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

lasso_drdid=lasso_drdid_rc(y=Y, post=T, D=D, covariates=X, i.weights = NULL, boot = FALSE,
                           boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)
lasso_drdid_coeff=lasso_drdid$ATT
lasso_drdid_se=lasso_drdid$se

lasso_triple_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
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
  int.cov=as.matrix(cbind(int.cov,covariates[,13:(ncol(covariates))]))
  
  
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
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0)
  
  out.y.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1)
  out.y.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)* ps.fit/(1 - ps.fit)
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




lasso_triple_drdid_rc <-function(y, post, D, covariates, i.weights = NULL,
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
  int.cov=as.matrix(cbind(int.cov,covariates[,13:(ncol(covariates))]))
  
  
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
  #Compute the Outcome regression for the control group at the pre-treatment period, using LASSO.
  y_C0=y[(D==0) & (post==0)]
  int.cov_C0=int.cov[(D==0) & (post==0),]
  i.weights_C0=i.weights[(D==0) & (post==0)]
  reg.cont.coeff.pre <- cv.glmnet(int.cov_C0, y_C0, type.measure="mse",
                                  alpha=1, family="gaussian", weights=i.weights_C0)
  
  out.y.pre=predict(reg.cont.coeff.pre, s=reg.cont.coeff.pre$lambda.1se, newx=int.cov)
  
  #Compute the Outcome regression for the control group at the post-treatment period, using LASSO.
  y_C1=y[(D==0) & (post==1)]
  int.cov_C1=int.cov[(D==0) & (post==1),]
  i.weights_C1=i.weights[(D==0) & (post==1)]
  reg.cont.coeff.post <- cv.glmnet(int.cov_C1, y_C1, type.measure="mse",
                                   alpha=1, family="gaussian", weights=i.weights_C1)
  out.y.post=predict(reg.cont.coeff.post, s=reg.cont.coeff.post$lambda.1se, newx=int.cov)
  
  # Combine the ORs for control group
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  
  
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)* ps.fit/(1 - ps.fit)
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
  
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  # Influence function for the treated component
  inf.treat <- inf.treat.post - inf.treat.pre
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect from nuisance parameters
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)
  
  # Influence function for the control component
  inf.cont <- inf.cont.post - inf.cont.pre
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
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_imp_rc1,
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
    estMethod = "imp2",
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
  
  #warning("This estimator is not locally efficient. Consider using 'drdid_imp_rc' instead.")
  # return the list
  return(ret)
  
}

lasso_3wdrdid=lasso_triple_drdid_rc(y=Y, post=T, D=D, covariates=X, i.weights = NULL, boot = FALSE,
                                    boot.type =  "weighted",  nboot = NULL, inffunc = FALSE)

lasso_3wdrdid_coeff=lasso_3wdrdid$ATT
lasso_3wdrdid_se=lasso_3wdrdid$se

#construct a summary table

coeff=cbind(twfe_coeff, twfe_int_coeff,lasso_3IPWRA_coeff,rf_3IPWRA_coeff, lasso_drdid_coeff,lasso_3wdrdid_coeff)
se=cbind(twfe_se, twfe_int_se,lasso_3IPWRA_se,rf_3IPWRA_se, lasso_drdid_se,lasso_3wdrdid_se)


tab <- matrix(1:12, ncol=6, byrow=TRUE)
tab[1,]=coeff
tab[2,]=se
tab=round(tab, digits = 3)


colnames(tab) <- c('TWFE', 'TWFE_2','lasso 3IPWRA', 'RF 3IPWRA', 'lasso DRDID','lasso 3WDRDID')
rownames(tab) <- c('Coefficient','St.Err.')
latextable=stargazer(tab)
tab
#save table
write.table(tab, file=paste('sebribes','.txt',sep = ""))

#commands to read table and transfer to latex
#stargazer(read.table(paste('EXP0B','.txt',sep = "")), summary=FALSE)