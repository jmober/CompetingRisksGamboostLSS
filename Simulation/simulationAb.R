##############################################################################
#####         Modeling postoperative mortality in older patients         ##### 
#####          by boosting discrete-time competing risks models          #####
##############################################################################
#####                    Electronic Supplement                           #####
#####                    Author: Moritz Berger                           #####
##############################################################################
#####               Content: R-Code of the Simulations                   #####
#####                        presented in Section 5                      #####
##############################################################################

# The following contains R-Code to conduct additional simulation 
# scenario (b) with correlated covariates.

### Function one_sim() ### 

## Input 
# seed:   integer passed to set.seed() 
# n:      number of observations 
# k:      number of discrete time points 
# b:      parameter affecting the degree of censoring
# noninf: number of non-informative variables 
# mmax:   maximal number of boosting iterations 
# method: index of the method evaluated 

## Description
# The function 
# - generates three data sets according to the design described in Section 5,
#   a learning sample, validation sample and test sample 
# - fits each of the alternative methods (method 0 - 7) 
# - computes (i) the predictive log-likelihood (not reported in the manuscript),
#   (ii) the C-index, (iii) the time-dependent prediction error, 
#   (iv) the integrated prediction error, (v) the C-index without accounting 
#   for ties (not reported in the manuscript)
# 
# The function requires the R add-on packages "discSurv" , "gamboostLSS" 
# (included as tar.gz), "DCRvalidation" (included as tar.gz), 
# "MRSP" (included as tar.gz), "VGAM", "MASS", the self-implemented family function 
# "CSHM_family.R", functions to fit the tree-based approach by 
# Berger et. al (2019) and the functions implemented in "functionAb.R".

## Output 
# List with five elements: 
# - predictive log-likelihood value 
# - C-index value 
# - matrix of time-dependent prediction error 
# - integrated prediction error value 
# - C-index value without accounting for ties 

one_sim <- function(seed, n, k, b, noninf, mmax, method){
  
  # (optional) install packages 
  # install.packages("discSurv")
  # install.packages("gamboostLSS")
  # install.packages("http://cran.r-project.org/src/contrib/Archive/MRSP/MRSP_0.4.3.tar.gz", repos = NULL, type="source")
  # install.packages("VGAM")
  # install.packages("MASS")
  
  # load packages
  library("discSurv")
  library("gamboostLSS")
  library("DCRvalidation")
  library("MRSP")
  library("VGAM")
  library("MASS")
  
  source("./Packages_Functions/CSHM_family.R")
  source("./Packages_Functions/CRTreeDisc_all.R")
  source("./Packages_Functions/CRTreeDisc_fit.R")
  source("./Simulation/functionAb.R")
  
  # set coefficients
  beta1 <- c(0.8, -0.8, 0.6, -0.6)
  beta2 <- c(0.8, -0.8, 0.6, -0.6)
  
  set.seed(seed)
  gamma01 <- runif(k, -3, -2)
  gamma02 <- runif(k, -3, -2)
  
  # generate data sets
  ok  <- FALSE
  while(!ok){
    dat <- oneDat(n, k, beta1, beta2, b, gamma01, gamma02, noninf)
    if(length(unique(dat$time))==k & length(unique(dat$status))==3 && table(dat$time, dat$status)[k,1]>0){
      ok <- TRUE
    }
  }
  
  ok  <- FALSE
  while(!ok){
    dat_val <- oneDat(n, k, beta1, beta2, b, gamma01, gamma02, noninf)
    if(length(unique(dat_val$time))==k & length(unique(dat_val$status))==3){
      ok <- TRUE
    }
  } 
  
  ok  <- FALSE
  while(!ok){
    dat_test <- oneDat(n, k, beta1, beta2, b, gamma01, gamma02, noninf)
    if(length(unique(dat_test$time))==k & length(unique(dat_test$status))==3){
      ok <- TRUE
    }
  } 
  
  datLong <- dataLongCompRisks(dataShort=dat, timeColumn="time", eventColumns=c("event1", "event2"), timeAsFactor=T)
  datLong$y <- apply(datLong[,c(3:5)],1,function(x) which(x==1)-1)
  datLong$y <- factor(datLong$y)
  datLong_val <- dataLongCompRisks(dataShort=dat_val, timeColumn="time", eventColumns=c("event1", "event2"), timeAsFactor=T)
  yaug_val    <- datLong_val[,c(3:5)]
  datLong_test <- dataLongCompRisks(dataShort=dat_test, timeColumn="time", eventColumns=c("event1", "event2"), timeAsFactor=T)
  yaug_test    <- datLong_test[,c(3:5)]
  datLong_pred <- dataLongCompRisks(dataShort=dat_test, timeColumn="time", eventColumns=c("event1", "event2"), timeAsFactor=T,
                                    aggTimeFormat=T, lastTheoInt=10, responseAsFactor=T)
  
  names <- paste0("x",1:(8+noninf))
  datLong[,c(names)] <- scale(datLong[,c(names)], scale=F)
  datLong_val[,c(names)] <- scale(datLong_val[,c(names)], scale=F)
  datLong_test[,c(names)] <- scale(datLong_test[,c(names)], scale=F)
  datLong_pred[,c(names)] <- scale(datLong_pred[,c(names)], scale=F)
  
  # censoring distribution 
  formC    <- formula(paste("y~","timeInt"))
  datC     <- dataCensoring(dataShort=dat, timeColumn="time", eventColumns=c("event1", "event2"))
  datLongC <- dataLong(dataShort=datC, timeColumn="timeCens", eventColumn="yCens", timeAsFactor=T)
  modC     <- glm(formC, data=datLongC, family=binomial(link="logit"))
  predC    <- predict(modC, type="response", newdata=datLong_pred)
  lambdaC  <- matrix(predC, nrow=n, ncol=k, byrow=T)
  SC       <- cbind(1,t(apply(1-lambdaC, 1, cumprod)))
  
  # null-model
  margFit <- VGAM::vgam(cbind(e0, e1, e2) ~ timeInt,
                  family=VGAM::multinomial(parallel=FALSE, refLevel=1), 
                  data = datLong)
  margHaz <- VGAM::predictvglm(margFit, datLong_pred[1:k,], type="response")[,-1]
  overallHaz <- rowSums(margHaz)
  Smarg   <- cumprod(1-overallHaz)
  margProbs <- sapply(1:2, function(j) margHaz[,j]*c(1,Smarg[1:(k-1)]))
  
  
  
  ### True
  if(method==0){
    lambda10 <- attributes(dat_test)$lambda1
    lambda20 <- attributes(dat_test)$lambda2
    lambda0  <- cbind(1-c(t(lambda10))-c(t(lambda20)),c(t(lambda10)),c(t(lambda20)))
    colnames(lambda0) <- c("e0","e1","e2") 
    C0       <- calc_CI(lambda10, lambda20, dat_test, dat_test)
    C0H      <- calc_CIH(lambda0, datLong_pred, dat_test)
    e1_true <- unlist(sapply(1:n, function(j) attributes(dat_test)$lambda1[j,1:dat_test$time[j]]))
    e2_true <- unlist(sapply(1:n, function(j) attributes(dat_test)$lambda2[j,1:dat_test$time[j]]))
    e0_true <- 1-e1_true-e2_true
    e_true  <- cbind(e0_true, e1_true, e2_true)
    ll0     <- sum(yaug_test*log(e_true))
    PE0     <- calc_PE(lambda10, lambda20, dat_test, SC)
    IPE0    <- calc_IPE(PE0, margProbs, dat_test)
    to_return <- list(ll0,C0H,PE0,IPE0,C0)
  }
  
  ### GB logL
  if(method==1){
    datLong_fit <- rbind(datLong, datLong)
    datLong_fit$y <- c(datLong$e1, datLong$e2)

    form1  <- paste0("bols(", names, ",intercept=F)", collapse="+")
    form1 <- formula(paste0("y~bols(timeInt, intercept=T, contrasts.arg='contr.dummy')+", form1))
    
    mod1 <- gamboostLSS(formula = list(eta1 = form1,
                                       eta2 = form1),
                        families = CSHM(stabilization = "none"), data = datLong_fit, method = "noncyclic", 
                        control = boost_control(mstop = mmax, trace = F, nu = c(eta1 = 0.1, eta2 = 0.1)))
    eta11_val <- sapply(1:mmax, function(j) predict(mod1[j], type="link", newdata=datLong_val, parameter="eta1"))
    eta21_val <- sapply(1:mmax, function(j) predict(mod1[j], type="link", newdata=datLong_val, parameter="eta2"))
    lambda11_val <- sapply(1:mmax, function(j) exp(eta11_val[,j])/(1+exp(eta11_val[,j])+exp(eta21_val[,j])))
    lambda21_val <- sapply(1:mmax, function(j) exp(eta21_val[,j])/(1+exp(eta11_val[,j])+exp(eta21_val[,j])))
    ll1_val      <- sapply(1:mmax, function(j) sum(yaug_val[,2]*log(lambda11_val[,j])+yaug_val[,3]*log(lambda21_val[,j])
                                                   +yaug_val[,1]*log(1-lambda11_val[,j]-lambda21_val[,j])))
    mstop1  <- which.max(ll1_val)
    modopt1 <- mod1[mstop1]
    eta11 <- predict(modopt1, type="link", newdata=datLong_test, parameter="eta1")
    eta21 <- predict(modopt1, type="link", newdata=datLong_test, parameter="eta2")
    lambda11 <- exp(eta11)/(1+exp(eta11)+exp(eta21))
    lambda21 <- exp(eta21)/(1+exp(eta11)+exp(eta21))
    ll1      <- sum(yaug_test[,2]*log(lambda11)+yaug_test[,3]*log(lambda21)+yaug_test[,1]*log(1-lambda11-lambda21))
    eta11_C <- predict(modopt1, type="link", newdata=datLong_pred, parameter="eta1")
    eta21_C <- predict(modopt1, type="link", newdata=datLong_pred, parameter="eta2")
    lambda11_C <- exp(eta11_C)/(1+exp(eta11_C)+exp(eta21_C))
    lambda21_C <- exp(eta21_C)/(1+exp(eta11_C)+exp(eta21_C))
    lambda_hat1  <- cbind(1-lambda11_C-lambda21_C,lambda11_C,lambda21_C)
    colnames(lambda_hat1) <- c("e0","e1","e2")
    lambda11_C <- matrix(lambda11_C, nrow=n, ncol=k, byrow=T)
    lambda21_C <- matrix(lambda21_C, nrow=n, ncol=k, byrow=T)
    C1         <- calc_CI(lambda11_C, lambda21_C, dat, dat_test)
    C1H        <- calc_CIH(lambda_hat1, datLong_pred, dat_test)
    PE1        <- calc_PE(lambda11_C, lambda21_C, dat_test, SC)
    IPE1       <- calc_IPE(PE1, margProbs, dat_test)
    rm(mod1)
    to_return <- list(ll1,C1H,PE1,IPE1,C1)
  }
  
  ### GBg
  if(method==2){
    datLong_fit2 <- as.list(datLong)
    datLong_fit2$e  <- factor(c("1","2"))
    datLong_val2 <- as.list(datLong_val)
    datLong_val2$e <- factor(c("1","2"))
    datLong_test2 <- as.list(datLong_test)
    datLong_test2$e <- factor(c("1","2"))
    datLong_pred2 <- as.list(datLong_pred)
    datLong_pred2$e <- factor(c("1","2"))
    
    form2  <- paste0("bols(", names, ", intercept=F) %O% bols(e, contrasts.arg='contr.dummy')", collapse="+")
    form2 <- formula(paste0("y~bols(timeInt, intercept=T, contrasts.arg='contr.dummy')  %O% bols(e, contrasts.arg='contr.dummy')+", form2))
    
    mod2 <- gamboost(form2, data=datLong_fit2, family=Multinomial(), 
                     control = boost_control(mstop = mmax, trace = F, nu = 0.1))
    lambda_hat2_val <- lapply(1:mmax, function(j) predict(mod2[j], type="response", newdata=datLong_val2))
    ll2_val         <- sapply(1:mmax, function(j) sum(yaug_val*log(lambda_hat2_val[[j]])))
    mstop2 <- which.max(ll2_val)
    modopt2 <- mod2[mstop2]
    lambda_hat2 <- predict(modopt2, type="response", newdata=datLong_test2)
    ll2 <- sum(yaug_test*log(lambda_hat2))
    lambda_hat2 <- predict(modopt2, type="response", newdata=datLong_pred2)
    colnames(lambda_hat2) <- c("e0","e1","e2")
    lambda12_C  <- matrix(lambda_hat2[,2], nrow=n, ncol=k, byrow=T)
    lambda22_C  <- matrix(lambda_hat2[,3], nrow=n, ncol=k, byrow=T)
    C2          <- calc_CI(lambda12_C, lambda22_C, dat, dat_test)
    C2H         <- calc_CIH(lambda_hat2, datLong_pred, dat_test)
    PE2         <- calc_PE(lambda12_C, lambda22_C, dat_test, SC)
    IPE2        <- calc_IPE(PE2, margProbs, dat_test)
    rm(mod2)
    to_return <- list(ll2,C2H,PE2,IPE2,C2)
  }
  
  ### Null
  if(method==3){
    mod3 <- vgam(cbind(e0, e1, e2) ~ timeInt,
                 family=multinomial(parallel=FALSE, refLevel=1), 
                 data = datLong)
    lambda_hat3  <- predict(mod3, datLong_test[,c(2,10:(17+noninf))], type="response")
    ll3          <- sum(yaug_test*log(lambda_hat3))
    lambda_hat3 <- predict(mod3, type="response", newdata=datLong_pred)
    colnames(lambda_hat3) <- c("e0","e1","e2")
    lambda13_C  <- matrix(lambda_hat3[,2], nrow=n, ncol=k, byrow=T)
    lambda23_C  <- matrix(lambda_hat3[,3], nrow=n, ncol=k, byrow=T)
    C3          <- calc_CI(lambda13_C, lambda23_C, dat, dat_test) 
    C3H         <- calc_CIH(lambda_hat3, datLong_pred, dat_test)
    PE3         <- calc_PE(lambda13_C, lambda23_C, dat_test, SC)
    IPE3        <- calc_IPE(PE3, margProbs, dat_test)
    rm(mod3)
    to_return <- list(ll3,C3H,PE3,IPE3,C3)
  }
  
  ### PLg
  if(method==4){
    form4 <- paste0(names,collapse="+")
    form4 <- formula(paste("cbind(e0, e1, e2) ~","timeInt+", form4))
    
    fit4  <- MRSP(formula=form4, data=datLong, model=multinomlogit(), constr=1, 
                  penalty=TRUE, group.classes=TRUE, group.dummies=TRUE, adaptive="ML")
    lambda_hat4_val <- lapply(1:50, function(j) predict_MRSP(fit4[[j]], datLong_val[,c(2,10:(17+noninf))], noninf))
    ll4_val         <- sapply(1:50, function(j) sum(yaug_val*log(lambda_hat4_val[[j]])))
    lambdaopt4 <- which.max(ll4_val)
    mod4       <- fit4[[lambdaopt4]]
    lambda_hat4 <- predict_MRSP(mod4, datLong_test[,c(2,10:(17+noninf))], noninf)
    ll4         <- sum(yaug_test*log(lambda_hat4))
    lambda_hat4 <- predict_MRSP(mod4, datLong_pred, noninf)
    colnames(lambda_hat4) <- c("e0","e1","e2")
    lambda14_C  <- matrix(lambda_hat4[,2], nrow=n, ncol=k, byrow=T)
    lambda24_C  <- matrix(lambda_hat4[,3], nrow=n, ncol=k, byrow=T)
    C4          <- calc_CI(lambda14_C, lambda24_C, dat, dat_test) 
    C4H         <- calc_CIH(lambda_hat4, datLong_pred, dat_test)
    PE4        <- calc_PE(lambda14_C, lambda24_C, dat_test, SC)
    IPE4        <- calc_IPE(PE4, margProbs, dat_test)
    rm(mod4)
    to_return <- list(ll4,C4H,PE4,IPE4,C4)
  }
  
  ### PL
  if(method==5){
    form5 <- paste0(names,collapse="+")
    form5 <- formula(paste("cbind(e0, e1, e2) ~","timeInt+", form5))
    
    fit5  <- MRSP(formula=form5, data=datLong, model=multinomlogit(), constr=1, 
                  penalty=TRUE, group.classes=FALSE, group.dummies=TRUE, adaptive="ML")
    lambda_hat5_val <- lapply(1:50, function(j) predict_MRSP(fit5[[j]], datLong_val[,c(2,10:(17+noninf))], noninf))
    ll5_val         <- sapply(1:50, function(j) sum(yaug_val*log(lambda_hat5_val[[j]])))
    lambdaopt5 <- which.max(ll5_val)
    mod5       <- fit5[[lambdaopt5]]
    lambda_hat5 <- predict_MRSP(mod5, datLong_test[,c(2,10:(17+noninf))], noninf)
    ll5         <- sum(yaug_test*log(lambda_hat5))
    lambda_hat5 <- predict_MRSP(mod5, datLong_pred, noninf)
    colnames(lambda_hat5) <- c("e0","e1","e2")
    lambda15_C  <- matrix(lambda_hat5[,2], nrow=n, ncol=k, byrow=T)
    lambda25_C  <- matrix(lambda_hat5[,3], nrow=n, ncol=k, byrow=T)
    C5          <- calc_CI(lambda15_C, lambda25_C, dat, dat_test) 
    C5H         <- calc_CIH(lambda_hat5, datLong_pred, dat_test)
    PE5         <- calc_PE(lambda15_C, lambda25_C, dat_test, SC)
    IPE5        <- calc_IPE(PE5, margProbs, dat_test)
    rm(mod5)
    to_return <- list(ll5,C5H,PE5,IPE5,C5)
  }
  
  ### Tree 
  if(method==6){
    form6 <- paste0(names, collapse="+")
    form6 <- formula(paste0("time~", form6))
    
    mod6 <- CRTreeDisc(formula=form6,
                       data=dat,
                       distance="Gini",
                       eventColumns=c("event1","event2"),
                       minimal_ns=seq(1,nrow(datLong), by=50), 
                       trace=F)
    datLong_val$timeInt <- as.numeric(datLong_val$timeInt)
    lambda_hat6_val  <- lapply(1:length(mod6), function(j) CRTreeDisc_predict(mod6[[j]],datLong_val[,c(2, 10:(17+noninf))])$prob)
    ll6_val          <- sapply(1:length(mod6), function(j) sum(yaug_val*log(lambda_hat6_val[[j]])))
    minimal_ns6      <- which.max(ll6_val)
    modopt6          <- mod6[[minimal_ns6]]
    datLong_test$timeInt <- as.numeric(datLong_test$timeInt)
    lambda_hat6  <- CRTreeDisc_predict(modopt6,datLong_test[,c(2, 10:(17+noninf))])$prob
    ll6 <- sum(yaug_test*log(lambda_hat6))
    lambda_hat6 <- CRTreeDisc_predict(modopt6,datLong_pred)$prob
    colnames(lambda_hat6) <- c("e0","e1","e2")
    lambda16_C  <- matrix(lambda_hat6[,2], nrow=n, ncol=k, byrow=T)
    lambda26_C  <- matrix(lambda_hat6[,3], nrow=n, ncol=k, byrow=T)
    C6          <- calc_CI(lambda16_C, lambda26_C, dat, dat_test)
    C6H         <- calc_CIH(lambda_hat6, datLong_pred, dat_test)
    PE6         <- calc_PE(lambda16_C, lambda26_C, dat_test, SC)
    IPE6        <- calc_IPE(PE6, margProbs, dat_test)
    rm(mod6)
    to_return <- list(ll6,C6H,PE6,IPE6,C6)
  }
  
  ### GB StabS with pi=0.7
  if(method==7){
    datLong_fit <- rbind(datLong, datLong)
    datLong_fit$y <- c(datLong$e1, datLong$e2)
    
    form7  <- paste0("bols(", names, ",intercept=F)", collapse="+")
    form7 <- formula(paste0("y~bols(timeInt, intercept=T, contrasts.arg='contr.dummy')+", form7))
    
    mod7 <- gamboostLSS(formula = list(eta1 = form7,
                                       eta2 = form7),
                        families = CSHM(stabilization = "none"), data = datLong_fit, method = "noncyclic", 
                        control = boost_control(mstop = mmax, trace = F, nu = c(eta1 = 0.1, eta2 = 0.1)))
    qq <- ceiling(sqrt(10*(noninf+9)))
    stab7 <- stabsel(mod7, cutoff=0.6, q=qq, sampling.type="MB", B=50)
    sel7 <- final_fit(stab7, 0.7)
    modopt7 <- vgam(sel7[[1]], data=datLong, family=multinomial(parallel=FALSE, refLevel=1), constraints=sel7[[2]])
    lambda_hat7 <- predict(modopt7, datLong_test[,c(2,10:(17+noninf))], type="response")
    ll7         <- sum(yaug_test*log(lambda_hat7))
    lambda_hat7 <- predict(modopt7, type="response", newdata=datLong_pred)
    lambda17_C  <- matrix(lambda_hat7[,2], nrow=n, ncol=k, byrow=T)
    lambda27_C  <- matrix(lambda_hat7[,3], nrow=n, ncol=k, byrow=T)
    C7          <- calc_CI(lambda17_C, lambda27_C, dat, dat_test) 
    C7H         <- calc_CIH(lambda_hat7, datLong_pred, dat_test)
    PE7         <- calc_PE(lambda17_C, lambda27_C, dat_test, SC)
    IPE7        <- calc_IPE(PE7, margProbs, dat_test)
    rm(mod7)
    to_return <- list(ll7,C7H,PE7,IPE7,C7)
  }
  
  return(to_return)
}

###########################

### Execution ### 

# The results of additional simulation scenario (b) for the seven methods 
# GB logL, GB StabS, GBg, PL, PLg, Tree, Null and True are stored in 
# "./raw_results/simulationAb.rda" (as a list with 2400 elements). 

# The function one_sim() was called 2,400 times with values n=250, k=10, b=0.9 and the following other arguments: 
args <- expand.grid(seed=c(1:100)+160920, 
                    noninf=c(10,200,1000),
                    mmax=15000, 
                    method=c(0:7))
args$mmax[args$noninf==200] <- 12000
args$mmax[args$noninf==1000] <- 10000

## Example 1: seed=160921, noninf=10, mmax=15000, method=0 (True)
one_sim(seed=160921, n=250, k=10, b=0.9, noninf=10, mmax=15000, method=0)

# result is stored in 
load("./Simulation/raw_results/simulationAb.rda")
res[[1]]

## Example 2: seed=160922, noninf=10, mmax=15000, method=3 (Null)
one_sim(seed=160922, n=250, k=10, b=0.9, noninf=10, mmax=15000, method=3)

# result is stored in 
res[[902]]

## Example 3: seed=160924, noninf=200, mmax=120, method=1 (GB logL)
one_sim(seed=160924, n=250, k=10, b=0.9, noninf=200, mmax=120, method=1)

# corresponding result with mmax=12000 is stored in 
res[[404]]


###########################

# The whole simulation was executed on a Linux Cluster by use of R-version 
# 3.5.1 and R add-on package "batchtools". 
# The original call was the following: 

# library("batchtools")
# 
# r <- makeRegistry(file.dir="./Simulation/simulationAb")
# # r <- loadRegistry(file.dir="./Simulation/simulationAb")
# 
# args <- expand.grid(seed=c(1:100)+160920, 
#                     n=250,
#                     k=10,
#                     b=0.9,
#                     noninf=c(10,200,1000),
#                     mmax=15000, 
#                     method=c(0:6))
# args$mmax[args$noninf==200] <- 12000
# args$mmax[args$noninf==1000] <- 10000
# 
# batchMap(one_sim, args=args, reg=r)
# submitJobs(1:2400, reg = r)



