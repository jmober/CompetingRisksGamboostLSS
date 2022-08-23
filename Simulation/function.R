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

# The following contains auxiliary functions to conduct the second part of the
# simulation study ('Comparison to alternative methods'). 


# one discrete data set 
oneDat <- function(n, k, beta1, beta2, b, gamma01, gamma02, noninf=NULL){

  x1 <- rnorm(n,0,1)
  x2 <- rnorm(n,0,1)
  x3 <- rnorm(n,0,1)
  x4 <- rnorm(n,0,1)
  x5 <- rnorm(n,0,1)
  x6 <- rnorm(n,0,1)
  x7 <- rnorm(n,0,1)
  x8 <- rnorm(n,0,1)
  X  <- cbind(x1,x2,x3,x4,x5,x6,x7,x8)
  if(!is.null(noninf)){
    for(i in 1:noninf){
      x <- rnorm(n,0,1)
      X <- cbind(X,x)
    }
    colnames(X)[9:(8+noninf)] <- paste0("x",9:(8+noninf))
  }

  eta1 <- sapply(1:k, function(j) gamma01[j]+beta1[1]*x1+beta1[2]*x2+beta1[3]*x5+beta1[4]*x6)
  eta2 <- sapply(1:k, function(j) gamma02[j]+beta2[1]*x3+beta2[2]*x4+beta2[3]*x7+beta2[4]*x8) 

  lambda1 <- exp(eta1)/(1+exp(eta1)+exp(eta2))
  lambda2 <- exp(eta2)/(1+exp(eta1)+exp(eta2))
  lambda  <- cbind(lambda1+lambda2,1)
  Ssim    <- t(apply(lambda,1,function(x) cumprod(1-x)))
  probs   <- cbind(1,Ssim[,1:k])*lambda

  t_disc  <- sapply(1:n, function(j) sample(k+1,1,prob=probs[j,]))

  cens_dist <- function(b){
    b_pot <- b^((k+1):1)
    ps    <- b_pot/sum(b_pot)
  }

  C_disc <- sample(1:(k+1), n, prob=cens_dist(b), replace=TRUE)

  time   <- sapply(1:n, function(j) min(C_disc[j], t_disc[j]))
  status <- ifelse(C_disc<t_disc, 0, 1)
  
  event <- sapply(1:n, function(j) {
    if(time[j]<=k & status[j]==1){
      p1 <- lambda1[j,t_disc[j]]/(lambda1[j,t_disc[j]]+lambda2[j,t_disc[j]])
      p2 <- lambda2[j,t_disc[j]]/(lambda1[j,t_disc[j]]+lambda2[j,t_disc[j]])
      ev <- rbinom(1,1,p2)+1
    }else{
      ev <- 0 
    }
    return(ev)
  })
  time[time==(k+1)] <- k
  
  t_data <- data.frame(time=time, status=event)
  e_data <- model.matrix(~as.factor(status), data=t_data)[,-1]
  
  disc_data <- cbind(e_data, t_data, X)
  disc_data <- as.data.frame(disc_data)
  names(disc_data)[c(1:2)] <- c("event1","event2")
  
  attr(disc_data, "lambda1") <- lambda1
  attr(disc_data, "lambda2") <- lambda2
  return(disc_data)
}

# calculate C-index based on R-package discSurv (not accounting for ties)
calc_CI <- function(lambda1, lambda2, data, data_test){

  cl_hat1 <- apply(lambda1, 1, cumsum)
  cl_hat2 <- apply(lambda2, 1, cumsum)
  marker1 <- apply(cl_hat1, 2, sum)
  marker2 <- apply(cl_hat2, 2, sum)
  C1      <- cIndex(marker1, data_test$time, data_test$event1,
                              data$time, data$event1)
  C2      <- cIndex(marker2, data_test$time, data_test$event2,
                              data$time, data$event2)
  nevents <- sum(data_test$event1)+sum(data_test$event2)
  Cind    <- C1*sum(data_test$event1)/nevents+C2*sum(data_test$event2)/nevents
  return(Cind)
}

# calculate C-index based on R-package DCRvalidation (accounting for ties)
calc_CIH <- function(lambda, data_pred, data_test){
  
  data_pred$timeInt <- as.numeric(data_pred$timeInt)
  
  data_pred$responses2 <- data_pred$responses
  for(i in 1:250){
    timei  <- data_pred$time[data_pred$obj==i][1]
    indexi <- (c(1:10)+10*(i-1))[(timei+1):10]
    data_pred$responses2[indexi] <- data_pred$responses[indexi[1]-1]
  }
  
  C1  <- discrimination(predictedProbs = lambda, 
                        time_points = 1:10, 
                        event_of_interest = "e1",
                        data = data_pred, 
                        name_outcome = 'responses2',
                        name_timeVar = 'timeInt')
  C2  <- discrimination(predictedProbs = lambda, 
                        time_points = 1:10, 
                        event_of_interest = "e2",
                        data = data_pred, 
                        name_outcome = 'responses2',
                        name_timeVar = 'timeInt')
  
  nevents <- sum(data_test$event1)+sum(data_test$event2)
  Cind    <- C1$Cindex*sum(data_test$event1)/nevents+C2$Cindex*sum(data_test$event2)/nevents
  return(Cind)
}

# calculate time-dependent prediction error  
calc_PE <- function(lambda1, lambda2, data_test, Ghat){

  n      <- nrow(data_test)
  k      <- ncol(lambda1)
  l_hats <- lambda1+lambda2
  S_hat  <- t(apply(1-l_hats, 1, cumprod))
  F_hat1 <- t(sapply(1:n, function(j) cumsum(c(1,S_hat[j,1:(k-1)])*lambda1[j,])))
  F_hat2 <- t(sapply(1:n, function(j) cumsum(c(1,S_hat[j,1:(k-1)])*lambda2[j,])))
  
  weights <- t(sapply(1:n, function(j) I(data_test$status[j]>0)*((data_test[j,"time"]<=c(1:k))*1)/Ghat[j,data_test[j,"time"]]+
                        (data_test[j,"time"]>(1:k))*1/Ghat[j,-1]))
  
  PE_1   <- weights*t(sapply(1:n, function(j) (c(data_test[j,"time"]<=c(1:k) & data_test[j,"status"]==1)*1- F_hat1[j,])^2))
  PE_1   <- apply(PE_1,2,mean)
  PE_2   <- weights*t(sapply(1:n, function(j) (c(data_test[j,"time"]<=c(1:k) & data_test[j,"status"]==2)*1- F_hat2[j,])^2))
  PE_2   <- apply(PE_2,2,mean)
  return(cbind(PE_1,PE_2))
}

# calculate integrated prediction error 
calc_IPE <- function(PE, marg_probs, data_test){
  IPEs <- colSums(PE*marg_probs)
  nevents <- sum(data_test$event1)+sum(data_test$event2)
  IPE  <- IPEs[1]*sum(data_test$event1)/nevents+IPEs[2]*sum(data_test$event2)/nevents
  return(IPE)
}

# calculate predictions from object of class MRSP (penalized ML approach)
predict_MRSP <- function(model, Xnew, noninf){
  form <- formula(paste0("~timeInt+",paste0("x",1:(8+noninf),collapse="+")))
  X <- model.matrix(form, data=Xnew)
  eta1 <- X%*%model@coef[[1]][2,]
  eta2 <- X%*%model@coef[[1]][3,]
  lambda1 <- exp(eta1)/(1+exp(eta1)+exp(eta2))
  lambda2 <- exp(eta2)/(1+exp(eta1)+exp(eta2))
  lambda  <- 1-lambda1-lambda2
  return(cbind(lambda,lambda1,lambda2))
}

# function to calculate C-index (also part of the R-package DCRvalidation)
discrimination <- function (predictedProbs, time_points, event_of_interest, data, 
                            name_outcome, name_timeVar, conditional_to = NULL) 
{
  AUC_tp <- sapply(time_points, function(t) {
    t_cases <- which(data[, name_timeVar] == t & data[, name_outcome] == 
                       event_of_interest)
    t_controls <- which(data[, name_timeVar] == t & data[, 
                                                         name_outcome] != event_of_interest)
    cases_probs <- predictedProbs[t_cases, event_of_interest]
    controls_probs <- predictedProbs[t_controls, event_of_interest]
    if (length(cases_probs >= 1)) {
      ci_auc <- biostatUZH::confIntAUC(cases = cases_probs, 
                                       controls = controls_probs)
      return(ci_auc$AUC[1])
    }
    else return(NA)
  })
  AUC_tp_forC <- AUC_tp
  AUC_tp_forC[which(is.na(AUC_tp))] <- 0.5
  df <- data[, c(name_outcome, name_timeVar)]
  model <- VGAM::vglm(formula = paste0(name_outcome, "~ns(", 
                                       name_timeVar, ",df=4)"), data = df, family = multinomial(refLevel = 1))
  new.df <- data.frame((unique(df[, name_timeVar])))
  colnames(new.df) <- name_timeVar
  Pr <- predictvglm(model, newdata = new.df, type = "response")[, 
                                                                -1]
  Pr_r <- Pr[, event_of_interest]
  Pr_r_t <- (Pr_r * c(1, cumprod(1 - rowSums(Pr))[-nrow(Pr)]))[new.df[, 
                                                                      name_timeVar] %in% time_points]
  Pr_largert <- cumprod(1 - rowSums(Pr))[new.df[, name_timeVar] %in% 
                                           time_points]
  if (!is.null(conditional_to)) {
    Pr_r_t <- Pr_r * c(1, cumprod(1 - rowSums(Pr))[-nrow(Pr)])
    Pr_r_t <- Pr_r_t/c(rep(1, conditional_to), Pr_r_t[-((length(Pr_r_t) - 
                                                           conditional_to + 1):length(Pr_r_t))])
    Pr_largert <- cumprod(1 - rowSums(Pr))
    Pr_largert <- Pr_largert/c(rep(1, conditional_to), Pr_largert[-((length(Pr_largert) - 
                                                                       conditional_to + 1):length(Pr_largert))])
  }
  C <- sum(AUC_tp_forC * Pr_r_t * Pr_largert)/sum(Pr_r_t * 
                                                    Pr_largert)
  return(list(time_dependent_auc = AUC_tp, Cindex = C))
}

