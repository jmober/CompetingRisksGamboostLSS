##############################################################################
#####         Modeling postoperative mortality in older patients         ##### 
#####          by boosting discrete-time competing risks models          #####
##############################################################################
#####                    Electronic Supplement                           #####
#####                    Author: Moritz Berger                           #####
##############################################################################
#####               Content: R-Code of the Application                   #####
#####                        presented in Section 6                      #####
##############################################################################

# The following contains R-Code to conduct the second part of the 
# analysis of the POSE data ('Comparison to alternative methods').

### Function one_analysis() ### 

## Input 
# seed:   integer passed to set.seed() 
# method: index of the method evaluated 
# data:   analysis data set 

## Description
# The function 
# - randomly splits the analysis data into a learning and test sample
# - fits each of the alternative methods (method 1 - 6) 
# - computes (i) the C-index, (ii) the time-dependent prediction error, 
#   (iii) the integrated prediction error, (iv) coefficient estimates,
#   confidence intervals and the P-value of the recalibration assessment. 
# 
# The function requires the R add-on packages "discSurv" , "gamboostLSS" 
# (included as tar.gz), "DCRvalidation" (included as tar.gz), 
# "MRSP" (included as tar.gz), "VGAM", the self-implemented family function 
# "CSHM_family.R", functions to fit the tree-based approach by 
# Berger et. al (2019) and the functions implemented in "functionResampling.R".

## Output 
# List with eight elements (method 2,3,4,6) 
# - C-index value 
# - matrix of time-dependent prediction error 
# - integrated prediction error value 
# - vector of coefficient estimates
# - matrix with confidence intervals 
# - P-value of recalibration test 
# - C-index value (death)
# - C-index value (discharge)
# List with nine elements (method 1 and 5)
# - C-index value 
# - matrix of time-dependent prediction error 
# - integrated prediction error value 
# - vector of coefficient estimates
# - matrix with confidence intervals 
# - P-value of recalibration test 
# - C-index value (death)
# - C-index value (discharge)
# - number of fitted coefficients 


one_analysis <- function(seed, method, data){
  
  
  # load packages 
  library("discSurv")
  library("gamboostLSS")
  library("DCRvalidation")
  library("MRSP")
  library("VGAM")
  
  source("./Packages_Functions/CSHM_family.R")
  source("./Packages_Functions/CRTreeDisc.R")
  source("./Packages_Functions/CRTreeDisc_fit.R")
  source("./Application/functionResampling.R")
  
  # randomly split data
  pose_analysis <- na.omit(data)
  
  set.seed(seed)
  train_ids0 <- sample(which(pose_analysis$status30==0), 176, replace=FALSE)
  train_ids1 <- sample(which(pose_analysis$status30==1), 141, replace=FALSE)
  train_ids2 <- sample(which(pose_analysis$status30==2), 4423, replace=FALSE)
  train_ids  <- c(train_ids0,train_ids1,train_ids2)
  test_ids   <- seq(1,7110)[-train_ids]
  dat12 <- pose_analysis[train_ids,]
  dat3  <- pose_analysis[test_ids,]
  
  # prepare training sample
  dat12$death <- (dat12$status30==1)*1
  dat12$discharge <- (dat12$status30==2)*1
  dat12$int    <- 1
  dat12L <- dataLongCompRisks(dataShort=dat12, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=F)
  dat12L[,c("gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")] <- scale(dat12L[,c("gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")], scale=F)
  
  # prepare test sample
  dat3$death <- (dat3$status30==1)*1
  dat3$discharge <- (dat3$status30==2)*1
  dat3$int    <- 1
  dat3L <- dataLongCompRisks(dataShort=dat3, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=F)
  dat3L[,c("gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")] <- scale(dat3L[,c("gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")], scale=F)
  
  # for prediction
  n <- 2370
  k <- 31
  dat3L_pred <- dataLongCompRisks(dataShort=dat3, timeColumn="time30", eventColumns=c("death","discharge"), 
                                  timeAsFactor=F, aggTimeFormat=TRUE, lastTheoInt=31, responseAsFactor=TRUE)
  dat3L_pred[,c("gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")] <- scale(dat3L_pred[,c("gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")], scale=F)
  
  # censoring distribution
  dat3L_pred0 <- dat3L_pred
  dat3L_pred0$timeInt <- as.factor(dat3L_pred0$timeInt)
  formC    <- formula(paste("y~","timeInt"))
  datC     <- dataCensoring(dataShort=dat12, timeColumn="time30", eventColumns=c("death","discharge"))
  datLC <- dataLong(dataShort=datC, timeColumn="timeCens", eventColumn="yCens", timeAsFactor=T)
  modC     <- glm(formC, data=datLC, family=binomial(link="logit"))
  predC    <- predict(modC, type="response", newdata=dat3L_pred0)
  lambdaC  <- matrix(predC, nrow=n, ncol=k, byrow=T)
  Ghat     <- cbind(1,t(apply(1-lambdaC, 1, cumprod)))
  
  # null-model
  margFit <- VGAM::vgam(cbind(e0, e1, e2) ~ timeInt,
                        family=VGAM::multinomial(parallel=FALSE, refLevel=1), 
                        data = dat12L)
  margHaz <- VGAM::predictvglm(margFit, dat3L_pred[1:k,], type="response")[,-1]
  overallHaz <- rowSums(margHaz)
  Smarg   <- cumprod(1-overallHaz)
  margProbs <- sapply(1:2, function(j) margHaz[,j]*c(1,Smarg[1:(k-1)]))

  ### GB StabS with pi=0.7
  if(method==1){
    
    dat12L[,c("timeInt")] <- scale(dat12L[,c("timeInt")], scale=F)
    dat3L[,c("timeInt")] <- scale(dat3L[,c("timeInt")], scale=F)
    dat3L_pred1 <- dat3L_pred
    dat3L_pred[,c("timeInt")] <- scale(dat3L_pred[,c("timeInt")], scale=F)
    
    # variables for smooth-baselearners
    dat12L$ages <- dat12L$age
    dat12L$timeInts <- dat12L$timeInt
    dat3L$ages <- dat3L$age
    dat3L$timeInts <- dat3L$timeInt
    
    names <- names(dat12[,c(7:8,10:12,14,16:23)])
    dat12L_fit <- rbind(dat12L, dat12L)
    dat12L_fit$y   <- c(dat12L$death, dat12L$discharge)
    
    form  <- paste0("bols(", names, ",intercept=F)", collapse="+")
    form  <- formula(paste0("y~bols(int, intercept=F)+bols(timeInt, intercept=F)+bbs(timeInts, center=T, df=1)+", form, "+bbs(ages, center=T, df=1)"))
    
    # Fit
    mod <- gamboostLSS(formula = list(eta1 = form,
                                      eta2 = form),
                       families = CR2(stabilization = "none"), data = dat12L_fit, method = "noncyclic", 
                       control = boost_control(mstop = 10000, trace = F, nu = c(eta1 = 0.1, eta2 = 0.1)))
    coef(mod)
    set.seed(seed)
    stab <- stabsel(mod, cutoff=0.7, q=15, sampling.type="MB", B=300)
    
    ff <- final_fit(stab, 0.7)
    formf <- ff[[1]]
    clist <- ff[[2]]
    mod1  <- vgam(formf, data=dat12L, family=multinomial(parallel=FALSE, refLevel=1), constraints=clist)
    
    covars1 <- coef(mod1, matrix=T)
    ncoef1  <- length(which(covars1!=0))
    
    # C-Index
    dat3L_pred$ages <- dat3L_pred$age
    dat3L_pred$timeInts <- dat3L_pred$timeInt
    lambda_hat_pred <- predict(mod1, type="response", newdata=dat3L_pred)
    lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n, ncol=k, byrow=T)
    lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n, ncol=k, byrow=T)
    C1 <- calc_CIH(lambda_hat_pred, dat3L_pred1, dat3)
    
    # Prediction error 
    PE1  <- calc_PE(lambda1_pred, lambda2_pred, dat3, Ghat)
    IPE1 <- calc_IPE(PE1, margProbs, dat3)
    
    # Recalibration model 
    eta1_pred   <- predict(mod1, newdata=dat3L)[,1]
    eta2_pred   <- predict(mod1, newdata=dat3L)[,2]
    RM1 <- calc_RM(eta1_pred, eta2_pred, dat3L)
    
    rm(mod1)
    to_return <- list(C1[1],PE1,IPE1,RM1$coef,RM1$CI,RM1$pval,C1[2],C1[3],ncoef1)
  }  
  
  ### Null
  if(method==2){
    form2 <- formula(paste0("cbind(e0, e1, e2)~factor(timeInt)"))
    mod2 <- vglm(form2, data=dat12L, family=multinomial(parallel=FALSE, refLevel=1))
    
    lambda_hat_pred <- predict(mod2, type="response", newdata=dat3L_pred0)
    lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n, ncol=k, byrow=T)
    lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n, ncol=k, byrow=T)
    
    # C-Index
    C2 <- calc_CIH(lambda_hat_pred, dat3L_pred, dat3)
    
    # Prediction error 
    PE2  <- calc_PE(lambda1_pred, lambda2_pred, dat3, Ghat)
    IPE2 <- calc_IPE(PE2, margProbs, dat3)
    
    # Recalibration model 
    eta1_pred   <- predict(mod2, newdata=dat3L)[,1]
    eta2_pred   <- predict(mod2, newdata=dat3L)[,2]
    RM2 <- calc_RM(eta1_pred, eta2_pred, dat3L)
    
    rm(mod2)
    to_return <- list(C2[1],PE2,IPE2,RM2$coef,RM2$CI,RM2$pval,C2[2],C2[3])
  }
  
  ### Full
  if(method==3){
    names <- names(dat12[,c(7:8,10:12,14,16:23)])
    form3 <- paste0(names, collapse="+")
    form3 <- formula(paste0("cbind(e0, e1, e2)~factor(timeInt)+", form3))
    mod3 <- vglm(form3, data=dat12L, family=multinomial(parallel=FALSE, refLevel=1))
    
    lambda_hat_pred <- predict(mod3, type="response", newdata=dat3L_pred0)
    lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n, ncol=k, byrow=T)
    lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n, ncol=k, byrow=T)
    
    # C-Index
    C3 <- calc_CIH(lambda_hat_pred, dat3L_pred, dat3)
    
    # Prediction error 
    PE3  <- calc_PE(lambda1_pred, lambda2_pred, dat3, Ghat)
    IPE3 <- calc_IPE(PE3, margProbs, dat3)
    
    # Recalibration model 
    eta1_pred   <- predict(mod3, newdata=dat3L)[,1]
    eta2_pred   <- predict(mod3, newdata=dat3L)[,2]
    RM3 <- calc_RM(eta1_pred, eta2_pred, dat3L)
    
    rm(mod3)
    to_return <- list(C3[1],PE3,IPE3,RM3$coef,RM3$CI,RM3$pval,C3[2],C3[3])
  }
  
  ### PLg
  if(method==4){
    dat12L$timeInt <- as.factor(dat12L$timeInt)
    dat3L$timeInt <- as.factor(dat3L$timeInt)
    
    names <- names(dat12[,c(7:8,10:12,14,16:23)])
    form40 <- paste0(names,collapse="+")
    form4 <- formula(paste("cbind(e0, e1, e2) ~","timeInt+", form40))
    fold <- 10
    ngrid <- 100
    n12 <- nrow(dat12)
    t_index <- split(sample(1:n12), rep(1:fold, length = n12))  
    l_index <- lapply(t_index, function(i) setdiff(1:n12, i))
    
    ll_hat <- matrix(NA, nrow=fold, ncol=ngrid)
    for(kk in 1:fold){
      pose12_l <- dat12[l_index[[kk]],]
      poseLong_l <- dataLongCompRisks(dataShort=pose12_l, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=T)
      pose12_t <- dat12[t_index[[kk]],]
      poseLong_t <- dataLongCompRisks(dataShort=pose12_t, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=T)
      fit      <- MRSP(formula=form4, data=poseLong_l, model=multinomlogit(), constr=1, lambdamin=20, lambdamax=4000, nrlambda=ngrid,
                       penalty=TRUE, group.classes=TRUE, group.dummies=TRUE, adaptive="ML")
      lambda_hat <- lapply(1:ngrid, function(j) predict_MRSP(fit[[j]], poseLong_t, form40))
      ll_hat[kk,] <- sapply(1:ngrid, function(j) sum(poseLong_t[,c(3:5)]*log(lambda_hat[[j]])))
    }
    ll_hat_mean <- apply(ll_hat, 2, mean, na.rm=T)
    lambda_opt  <- fit[[which.max(ll_hat_mean)]]@lambda
    
    mod4      <- MRSP(formula=form4, data=dat12L, model=multinomlogit(), constr=1, lambda=lambda_opt,
                      penalty=TRUE, group.classes=TRUE, group.dummies=TRUE, adaptive="ML")
    
    lambda_hat_pred <- predict_MRSP(mod4, dat3L_pred0, form40)
    colnames(lambda_hat_pred) <- c("e0","e1","e2") 
    lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n, ncol=k, byrow=T)
    lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n, ncol=k, byrow=T)
    
    # C-Index
    C4 <- calc_CIH(lambda_hat_pred, dat3L_pred, dat3)
    
    # Prediction error 
    PE4  <- calc_PE(lambda1_pred, lambda2_pred, dat3, Ghat)
    IPE4 <- calc_IPE(PE4, margProbs, dat3)
    
    # Recalibration model 
    form <- formula(paste("~","timeInt+", form40))
    X <- model.matrix(form, data=dat3L)
    eta1_pred <- X%*%mod4@coef[[1]][2,]
    eta2_pred <- X%*%mod4@coef[[1]][3,]
    RM4 <- calc_RM(eta1_pred, eta2_pred, dat3L)
    
    rm(mod4)
    to_return <- list(C4[1],PE4,IPE4,RM4$coef,RM4$CI,RM4$pval,C4[2],C4[3])
  }
  
  ### PL
  if(method==5){
    dat12L$timeInt <- as.factor(dat12L$timeInt)
    dat3L$timeInt <- as.factor(dat3L$timeInt)
    
    names <- names(dat12[,c(7:8,10:12,14,16:23)])
    form50 <- paste0(names,collapse="+")
    form5 <- formula(paste("cbind(e0, e1, e2) ~","timeInt+", form50))
    fold <- 10
    ngrid <- 100 
    n12 <- nrow(dat12)
    t_index <- split(sample(1:n12), rep(1:fold, length = n12))  
    l_index <- lapply(t_index, function(i) setdiff(1:n12, i))
    
    ll_hat <- matrix(NA, nrow=fold, ncol=ngrid)
    for(kk in 1:fold){
      pose12_l <- dat12[l_index[[kk]],]
      poseLong_l <- dataLongCompRisks(dataShort=pose12_l, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=T)
      pose12_t <- dat12[t_index[[kk]],]
      poseLong_t <- dataLongCompRisks(dataShort=pose12_t, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=T)
      fit      <- MRSP(formula=form5, data=poseLong_l, model=multinomlogit(), constr=1, lambdamin=20, lambdamax=4000, nrlambda=ngrid,
                       penalty=TRUE, group.classes=FALSE, group.dummies=TRUE, adaptive="ML")
      lambda_hat <- lapply(1:ngrid, function(j) predict_MRSP(fit[[j]], poseLong_t, form50))
      ll_hat[kk,] <- sapply(1:ngrid, function(j) sum(poseLong_t[,c(3:5)]*log(lambda_hat[[j]])))
    }
    ll_hat_mean <- apply(ll_hat, 2, mean, na.rm=T)
    lambda_opt  <- fit[[which.max(ll_hat_mean)]]@lambda
    
    mod5      <- MRSP(formula=form5, data=dat12L, model=multinomlogit(), constr=1, lambda=lambda_opt,
                      penalty=TRUE, group.classes=FALSE, group.dummies=TRUE, adaptive="ML")
    
    covars5 <- mod5@coef[[1]]
    ncoef5  <- length(which(covars5!=0))
    
    lambda_hat_pred <- predict_MRSP(mod5, dat3L_pred0, form50)
    colnames(lambda_hat_pred) <- c("e0","e1","e2")
    lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n, ncol=k, byrow=T)
    lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n, ncol=k, byrow=T)
    
    # C-Index
    C5 <- calc_CIH(lambda_hat_pred, dat3L_pred, dat3)
    
    # Prediction error 
    PE5  <- calc_PE(lambda1_pred, lambda2_pred, dat3, Ghat)
    IPE5 <- calc_IPE(PE5, margProbs, dat3)
    
    # Recalibration model 
    form <- formula(paste("~","timeInt+", form50))
    X <- model.matrix(form, data=dat3L)
    eta1_pred <- X%*%mod5@coef[[1]][2,]
    eta2_pred <- X%*%mod5@coef[[1]][3,]
    RM5 <- calc_RM(eta1_pred, eta2_pred, dat3L)
    
    rm(mod5)
    to_return <- list(C5[1],PE5,IPE5,RM5$coef,RM5$CI,RM5$pval,C5[2],C5[3],ncoef5)
  }
  
  ### Tree 
  if(method==6){
    names <- names(dat12[,c(7:8,10:12,14,16:23)])
    form6 <- paste0(names,collapse="+")
    form6 <- formula(paste("time30~", form6))
    
    mod6 <- CRTreeDisc(formula=form6,
                       data=dat12,
                       tuning="ll",
                       distance="Gini",
                       eventColumns=c("death","discharge"),
                       minimal_ns=seq(1,nrow(dat12L), by=50), 
                       nfolds=10,
                       trace=F)
    
    lambda_hat_pred <- CRTreeDisc_predict(mod6[[1]], dat3L_pred)$prob
    colnames(lambda_hat_pred) <- c("e0","e1","e2")
    lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n, ncol=k, byrow=T)
    lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n, ncol=k, byrow=T)
    
    # C-Index
    C6 <- calc_CIH(lambda_hat_pred, dat3L_pred, dat3)
    
    # Prediction error 
    PE6  <- calc_PE(lambda1_pred, lambda2_pred, dat3, Ghat)
    IPE6 <- calc_IPE(PE6, margProbs, dat3)
    
    # Recalibration model 
    lambda_hat <- CRTreeDisc_predict(mod6[[1]], dat3L)$prob
    eta_hat <- log(lambda_hat/(1-lambda_hat))
    eta1_pred <- eta_hat[,2]
    eta2_pred <- eta_hat[,3]
    RM6 <- calc_RM(eta1_pred, eta2_pred, dat3L)
    
    rm(mod6)
    to_return <- list(C6[1],PE6,IPE6,RM6$coef,RM6$CI,RM6$pval,C6[2],C6[3])
  }
  return(to_return)
}

###########################

### Execution ### 

# The results of the second part of the analysis of the POSE data ('Comparison to alternative methods') for the six methods 
# GB StabS, PL, PLg, Tree, Null and Full are stored in "./raw_results/Resampling.rda" (as a list with 600 elements). 

# Note: For confidentiality reasons, the data set presented here is not 
# identical with the data used in the paper. Instead, the following R-code 
# employs an anonymized data set that was generated from the original data 
# in order to produce results that are similar (but not identical!) to those 
# in the paper. In particular, data lines in the anonymized data set do not 
# refer to real individuals. The results in "./raw_results/Resampling.rda" 
# correspond to the results obtained from the original data used in the paper. 

# The function one_analysis() was called 600 times with the arguments: 
args <- expand.grid(seed=c(1:100)+160920, 
                    method=c(1:6))

# load learning sample (n=4740; 2/3 of the anonymized data set)
load("./Application/data12_altered.rda")
# load test sample (n=2370; 1/3 of the anonymized data set)
load("./Application/data3_altered.rda")

# data set used for resampling
pose_analysis <- rbind(pose12_altered, pose3_altered)


## Example 1: seed=160921, method=2 (Null)
one_analysis(seed=160921, method=2, data=pose_analysis)

# the corresponding results obtained from the original data is stored in 
load("./Application/Resampling.rda")
res[[101]]

## Example 2: seed=160922, method=3 (Full)
one_analysis(seed=160922, method=3, data=pose_analysis)

# the corresponding results obtained from the original data is stored in 
res[[202]]


###########################

# The whole simulation was executed on a Linux Cluster by use of R-version 
# 3.5.1 and R add-on package "batchtools". 
# The original call was the following: 
#
# library("batchtools")
#
# r <- makeRegistry(file.dir="./Application/Resampling")
# # r <- loadRegistry(file.dir="./Application/Resampling")
# 
# args <- expand.grid(seed=c(1:100)+160920, 
#                     method=c(1:6), 
#                     data=pose_analysis)
# 
# batchMap(one_analysis, args=args, reg=r)
# submitJobs(1:600, reg = r)









