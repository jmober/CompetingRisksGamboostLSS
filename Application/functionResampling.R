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

# The following contains auxiliary functions to conduct the second part of the 
# analysis of the POSE data ('Comparison to alternative methods').

# evaluate selection of variables with GB StabS
final_fit <- function(stab, q){
  
  sel     <- names(which(stab$max>q))
  sel     <- strsplit(sel,"[.]")
  selvar  <- sapply(sel, `[[`, 1)
  selvarI <- selvar[!selvar%in%"int"]
  seleta  <- as.numeric(substr(sapply(sel, `[[`, 2)[!selvar%in%"int"],4,4))
  smooth  <- selvarI%in%c("timeInts","ages")
  selvarI[smooth] <- paste0("sm.ps(", selvarI[smooth],", outer.ok = T)")
  forml <- paste0(unique(selvarI[!smooth]), collapse="+")
  forms <- paste0(unique(selvarI[smooth]), collapse="+")
  if(forms==""){
    formf  <- formula(paste0("cbind(e0,e1,e2)~", forml))
  } else{
    formf  <- formula(paste0("cbind(e0,e1,e2)~", forml, "+", forms))
  }
  clist  <- vector(mode="list", length=length(unique(selvarI))+1)
  names(clist) <- c("(Intercept)", unique(selvarI[!smooth]), unique(selvarI[smooth]))
  clist[[1]] <- diag(2)
  for(i in names(clist)[-1]){
    if(sum(i==selvarI)==2){
      clist[[i]] <- diag(2)
    } else{ 
      if(seleta[i==selvarI]==1){
        clist[[i]] <- rbind(1,0)
      } else{
        clist[[i]] <- rbind(0,1)
      }
    } 
  }
  return(list(formf,clist))
}

# calculate C-index based on R-package discSurv (not accounting for ties)
calc_CI <- function(lambda1, lambda2, data, data_test){
  
  cl_hat1 <- apply(lambda1, 1, cumsum)
  cl_hat2 <- apply(lambda2, 1, cumsum)
  marker1 <- apply(cl_hat1, 2, sum)
  marker2 <- apply(cl_hat2, 2, sum)
  C1      <- cIndex(marker1, data_test$time30, data_test$death,
                    data$time30, data$death)
  C2      <- cIndex(marker2, data_test$time30, data_test$discharge,
                    data$time30, data$discharge)
  nevents <- sum(data_test$death)+sum(data_test$discharge)
  Cind    <- C1*sum(data_test$death)/nevents+C2*sum(data_test$discharge)/nevents
  return(c(Cind,C1,C2))
}

# calculate C-index based on R-package DCRvalidation (accounting for ties)
calc_CIH <- function(lambda, data_pred, data_test){
  
  data_pred$responses2 <- data_pred$responses
  for(i in 1:2370){
    timei  <- data_pred$time30[data_pred$obj==i][1]
    indexi <- (c(1:31)+31*(i-1))[(timei+1):31]
    data_pred$responses2[indexi] <- data_pred$responses[indexi[1]-1]
  }
  
  C1  <- discrimination(predictedProbs = lambda, 
                        time_points = 1:31, 
                        event_of_interest = "e1",
                        data = data_pred, 
                        name_outcome = 'responses2',
                        name_timeVar = 'timeInt')
  C2  <- discrimination(predictedProbs = lambda, 
                        time_points = 1:31, 
                        event_of_interest = "e2",
                        data = data_pred, 
                        name_outcome = 'responses2',
                        name_timeVar = 'timeInt')
  
  nevents <- sum(data_test$death)+sum(data_test$discharge)
  Cind    <- C1$Cindex*sum(data_test$death)/nevents+C2$Cindex*sum(data_test$discharge)/nevents
  return(c(Cind,C1$Cindex,C2$Cindex))
}

# calculate time-dependent prediction error 
calc_PE <- function(lambda1, lambda2, data_test, Ghat){
  
  n      <- nrow(data_test)
  k      <- ncol(lambda1)
  l_hats <- lambda1+lambda2
  S_hat  <- t(apply(1-l_hats, 1, cumprod))
  F_hat1 <- t(sapply(1:n, function(j) cumsum(c(1,S_hat[j,1:(k-1)])*lambda1[j,])))
  F_hat2 <- t(sapply(1:n, function(j) cumsum(c(1,S_hat[j,1:(k-1)])*lambda2[j,])))
  
  weights <- t(sapply(1:n, function(j) I(data_test$status30[j]>0)*((data_test[j,"time30"]<=c(1:k))*1)/Ghat[j,data_test[j,"time30"]]+
                        (data_test[j,"time30"]>(1:k))*1/Ghat[j,-1]))
  
  PE_1   <- weights*t(sapply(1:n, function(j) (c(data_test[j,"time30"]<=c(1:k) & data_test[j,"status30"]==1)*1- F_hat1[j,])^2))
  PE_1   <- apply(PE_1,2,mean)
  PE_2   <- weights*t(sapply(1:n, function(j) (c(data_test[j,"time30"]<=c(1:k) & data_test[j,"status30"]==2)*1- F_hat2[j,])^2))
  PE_2   <- apply(PE_2,2,mean)
  return(cbind(PE_1,PE_2))
}

# calculate integrated prediction error
calc_IPE <- function(PE, marg_probs, data_test){
  IPEs <- colSums(PE*marg_probs)
  nevents <- sum(data_test$death)+sum(data_test$discharge)
  IPE  <- IPEs[1]*sum(data_test$death)/nevents+IPEs[2]*sum(data_test$discharge)/nevents
  return(IPE)
}

# fit recalibration models 
calc_RM <- function(eta1, eta2, data_test_long){
  
  clistr <- list("(Intercept)" = diag(2),
                 "eta1" = rbind(1, 0),
                 "eta2" = rbind(0, 1))
  recal_mod_c  <- vglm(cbind(e0, e1, e2) ~ eta1+eta2, constraints=clistr,
                       data=data_test_long, family=multinomial(parallel=FALSE, refLevel=1))

  # likelihood under H_0 
  li <- sapply(1:nrow(data_test_long), function(x) sum(data_test_long[x,4:5]*c(eta1[x], eta2[x]))-log(1+exp(eta1[x])+exp(eta2[x])))
  ll0 <- sum(li)

  LRtest <- -2*(ll0-logLik(recal_mod_c))
  pval <- pchisq(LRtest, 4, lower.tail=F)

  return(list("coef"=coef(recal_mod_c), "CI"=confint(recal_mod_c), "pval"=pval))
}
  
# calculate predictions from object of class MRSP (penalized ML approach)
predict_MRSP <- function(model, Xnew, form){
  form <- formula(paste("~","timeInt+", form))
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
