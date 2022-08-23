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

# The following contains R-Code to conduct the first part of the
# simulation study ('Comparison of variable selection strategies'). 

### Function one_sim() ### 

## Input 
# seed:   integer passed to set.seed() 
# n:      number of observations 
# k:      number of discrete time points 
# b:      parameter affecting the degree of censoring
# noninf: number of non-informative variables 
# mmax:   maximal number of boosting iterations 
# B:      number of bootstrap replication for stability selection 

## Description
# The function 
# - generates three data sets according to the design described in Section 5,
#   a learning sample, validation sample and test sample 
# - fits each of the four GB StabS models with pi=0.6,0.7,0.8,0.9
# - computes (i) the predictive log-likelihood (not reported in the manuscript),
#   (ii) the C-index, (iii) the time-dependent prediction error, 
#   (iv) the integrated prediction error, (v) the C-index without accounting 
#   for ties (not reported in the manuscript)
# 
# The function requires the R add-on packages "discSurv" , "gamboostLSS" 
# (included as tar.gz), "DCRvalidation" (included as tar.gz), 
#  VGAM", the self-implemented family function 
# "CSHM_family.R" and the functions implemented in "functionGBStabs.R".

## Output 
# List with twenty elements (five for each model): 
# - predictive log-likelihood value 
# - C-index value 
# - matrix of time-dependent prediction error 
# - integrated prediction error value 
# - C-index value without accounting for ties 


one_sim <- function(seed, n, k, b, noninf, mmax, B){
  
  # load packages
  library("discSurv")
  library("gamboostLSS")
  library("DCRvalidation")
  library("VGAM")
  
  source("./Packages_Functions/CSHM_family.R")
  source("./Simulation/functionGBStabs.R")
  
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
  
  ### initialize with gamboostLSS
  datLong_fit <- rbind(datLong, datLong)
  datLong_fit$y <- c(datLong$e1, datLong$e2)

  form1  <- paste0("bols(", names, ",intercept=F)", collapse="+")
  form1 <- formula(paste0("y~bols(timeInt, intercept=T, contrasts.arg='contr.dummy')+", form1))
    
  mod1 <- gamboostLSS(formula = list(eta1 = form1,
                                       eta2 = form1),
                      families = CSHM(stabilization = "none"), data = datLong_fit, method = "noncyclic", 
                      control = boost_control(mstop = mmax, trace = F, nu = c(eta1 = 0.1, eta2 = 0.1)))
  
  ### stability selection
  qq <- ceiling(sqrt(10*(noninf+9)))
  stab1 <- stabsel(mod1, cutoff=0.6, q=qq, sampling.type="MB", B=B)
    
  # pi=0.6
  sel11 <- final_fit(stab1, 0.6)
  mod11 <- vgam(sel11[[1]], data=datLong, family=multinomial(parallel=FALSE, refLevel=1), constraints=sel11[[2]])
  lambda_hat11 <- predict(mod11, datLong_test[,c(2,10:(17+noninf))], type="response")
  ll11         <- sum(yaug_test*log(lambda_hat11))
  lambda_hat11 <- predict(mod11, type="response", newdata=datLong_pred)
  lambda111_C  <- matrix(lambda_hat11[,2], nrow=n, ncol=k, byrow=T)
  lambda211_C  <- matrix(lambda_hat11[,3], nrow=n, ncol=k, byrow=T)
  C11          <- calc_CI(lambda111_C, lambda211_C, dat, dat_test) 
  C11H         <- calc_CIH(lambda_hat11, datLong_pred, dat_test)
  PE11         <- calc_PE(lambda111_C, lambda211_C, dat_test, SC)
  IPE11        <- calc_IPE(PE11, margProbs, dat_test)
    
  # pi=0.7
  sel12 <- final_fit(stab1, 0.7)
  mod12 <- vgam(sel12[[1]], data=datLong, family=multinomial(parallel=FALSE, refLevel=1), constraints=sel12[[2]])
  lambda_hat12 <- predict(mod12, datLong_test[,c(2,10:(17+noninf))], type="response")
  ll12         <- sum(yaug_test*log(lambda_hat12))
  lambda_hat12 <- predict(mod12, type="response", newdata=datLong_pred)
  lambda112_C  <- matrix(lambda_hat12[,2], nrow=n, ncol=k, byrow=T)
  lambda212_C  <- matrix(lambda_hat12[,3], nrow=n, ncol=k, byrow=T)
  C12          <- calc_CI(lambda112_C, lambda212_C, dat, dat_test) 
  C12H         <- calc_CIH(lambda_hat12, datLong_pred, dat_test)
  PE12         <- calc_PE(lambda112_C, lambda212_C, dat_test, SC)
  IPE12        <- calc_IPE(PE12, margProbs, dat_test)
    
  # pi=0.8
  sel13 <- final_fit(stab1, 0.8)
  mod13 <- vgam(sel13[[1]], data=datLong, family=multinomial(parallel=FALSE, refLevel=1), constraints=sel13[[2]])
  lambda_hat13 <- predict(mod13, datLong_test[,c(2,10:(17+noninf))], type="response")
  ll13         <- sum(yaug_test*log(lambda_hat13))
  lambda_hat13 <- predict(mod13, type="response", newdata=datLong_pred)
  lambda113_C  <- matrix(lambda_hat13[,2], nrow=n, ncol=k, byrow=T)
  lambda213_C  <- matrix(lambda_hat13[,3], nrow=n, ncol=k, byrow=T)
  C13          <- calc_CI(lambda113_C, lambda213_C, dat, dat_test)
  C13H         <- calc_CIH(lambda_hat13, datLong_pred, dat_test)
  PE13         <- calc_PE(lambda113_C, lambda213_C, dat_test, SC)
  IPE13        <- calc_IPE(PE13, margProbs, dat_test)
    
  # pi=0.9
  sel14 <- final_fit(stab1, 0.9)
  mod14 <- vgam(sel14[[1]], data=datLong, family=multinomial(parallel=FALSE, refLevel=1), constraints=sel14[[2]])
  lambda_hat14 <- predict(mod14, datLong_test[,c(2,10:(17+noninf))], type="response")
  ll14         <- sum(yaug_test*log(lambda_hat14))
  lambda_hat14 <- predict(mod14, type="response", newdata=datLong_pred)
  lambda114_C  <- matrix(lambda_hat14[,2], nrow=n, ncol=k, byrow=T)
  lambda214_C  <- matrix(lambda_hat14[,3], nrow=n, ncol=k, byrow=T)
  C14          <- calc_CI(lambda114_C, lambda214_C, dat, dat_test) 
  C14H         <- calc_CIH(lambda_hat14, datLong_pred, dat_test)
  PE14         <- calc_PE(lambda114_C, lambda214_C, dat_test, SC)
  IPE14        <- calc_IPE(PE14, margProbs, dat_test)
    
  rm(mod1)
  to_return <- list(ll11,C11H,PE11,IPE11,C11,
                    ll12,C12H,PE12,IPE12,C12,
                    ll13,C13H,PE13,IPE13,C13,
                    ll14,C14H,PE14,IPE14,C14)

  return(to_return)
}

###########################

### Execution ### 

# The results of the first part of the simulation study ('Comparison of variable selection strategies') for the four 
# GB StabS models with pi=0.6,0.7,0.8,0.9 are stored in "./raw_results/simulationGBStabs.rda" (as a list with 300 elements). 
# The results of the method GB logL are stored in "./raw_results/simulation.rda".

# The function one_sim() was called 300 times with values n=250, k=10, b=0.9 and the following other arguments: 
args <- expand.grid(seed=c(1:100)+160920, 
                    noninf=c(10,200,1000),
                    mmax=10000, 
                    B=50)

## Modified example: seed=160921, noninf=10, mmax=1000, B=10
one_sim(seed=160921, n=250, k=10, b=0.9, noninf=10, mmax=1000, B=10)

# corresponding result with mmax=10000 and B=50 is stored in 
load("./Simulation/raw_results/simulationGBStabs.rda")
res[[1]]


###########################

# The whole simulation was executed on a Linux Cluster by use of R-version 
# 3.5.1 and R add-on package "batchtools". 
# The original call was the following: 

# library("batchtools")
# 
# r <- makeRegistry(file.dir="./Simulation/simulationGBStabs")
# # r <- loadRegistry(file.dir="./Simulation/simulationGBStabs")
# 
# args <- expand.grid(seed=c(1:100)+160920, 
#                     n=250,
#                     k=10,
#                     b=0.9,
#                     noninf=c(10,200,1000),
#                     mmax=10000,
#                     B=50)
# 
# batchMap(one_sim, args=args, reg=r)
# submitJobs(1:300, reg = r)




### results logL
load("./SimulationClean/simulation.rda")
reslogL <- res[301:600]

### results StabS
load("./SimulationClean/simulationGBStabs.rda")

### Table S1 (upper) ### 
noninf <- rep(c(10,200,1000), each=100)

# logL
C_logL <- sapply(1:300, function(j) reslogL[[j]][[2]])
aggregate(C_logL~noninf, mean, data=data.frame(C_logL,noninf))
aggregate(C_logL~noninf, sd, data=data.frame(C_logL,noninf))

# StabS 0.6
C_06 <- sapply(1:300, function(j) res[[j]][[2]])
aggregate(C_06~noninf, mean, data=data.frame(C_06,noninf))
aggregate(C_06~noninf, sd, data=data.frame(C_06,noninf))

# StabS 0.7
C_07 <- sapply(1:300, function(j) res[[j]][[7]])
aggregate(C_07~noninf, mean, data=data.frame(C_07,noninf))
aggregate(C_07~noninf, sd, data=data.frame(C_07,noninf))

# StabS 0.8
C_08 <- sapply(1:300, function(j) res[[j]][[12]])
aggregate(C_08~noninf, mean, data=data.frame(C_08,noninf))
aggregate(C_08~noninf, sd, data=data.frame(C_08,noninf))

# StabS 0.9
C_09 <- sapply(1:300, function(j) res[[j]][[17]])
aggregate(C_09~noninf, mean, data=data.frame(C_09,noninf))
aggregate(C_09~noninf, sd, data=data.frame(C_09,noninf))


### Table S1 (lower) ###

# logL
IPE_logL <- sapply(1:300, function(j) reslogL[[j]][[4]])
aggregate(IPE_logL~noninf, mean, data=data.frame(IPE_logL,noninf))*100
aggregate(IPE_logL~noninf, sd, data=data.frame(IPE_logL,noninf))*100

# StabS 0.6
IPE_06 <- sapply(1:300, function(j) res[[j]][[4]])
aggregate(IPE_06~noninf, mean, data=data.frame(IPE_06,noninf))*100
aggregate(IPE_06~noninf, sd, data=data.frame(IPE_06,noninf))*100

# StabS 0.7
IPE_07 <- sapply(1:300, function(j) res[[j]][[9]])
aggregate(IPE_07~noninf, mean, data=data.frame(IPE_07,noninf))*100
aggregate(IPE_07~noninf, sd, data=data.frame(IPE_07,noninf))*100

# StabS 0.8
IPE_08 <- sapply(1:300, function(j) res[[j]][[14]])
aggregate(IPE_08~noninf, mean, data=data.frame(IPE_08,noninf))*100
aggregate(IPE_08~noninf, sd, data=data.frame(IPE_08,noninf))*100

# StabS 0.9
IPE_09 <- sapply(1:300, function(j) res[[j]][[19]])
aggregate(IPE_09~noninf, mean, data=data.frame(IPE_09,noninf))*100
aggregate(IPE_09~noninf, sd, data=data.frame(IPE_09,noninf))*100













