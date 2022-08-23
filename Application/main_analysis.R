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

# The following contains R-Code to conduct the main part of the
# analysis of the POSE data.

# The R-Code provided can be used to reproduce Figure 3, Figure 4, Figure 5 and 
# Table 1 of the manuscript and Figure S7 of the supplementary material. 

# Note: For confidentiality reasons, the data set presented here is not 
# identical with the data used in the paper. Instead, the following R-code 
# employs an anonymized data set that was generated from the original data 
# in order to produce results that are similar (but not identical!) to those 
# in the paper. In particular, data lines in the anonymized data set do not 
# refer to real individuals.



### load packages
library("gamboostLSS")
library("discSurv")
library("VGAM")
library("e1071") 
library("DCRvalidation")

source("./Packages_Functions/CSHM_family.R")

# load learning sample (n=4740; 2/3 of the anonymized data set)
load("./Application/data12_altered.rda")
# load test sample (n=2370; 1/3 of the anonymized data set)
load("./Application/data3_altered.rda")

### prepare learning sample 
pose12_alteredL <- dataLongCompRisks(dataShort=pose12_altered, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=F)
names <- names(pose12_altered[,c(7:8,10:12,14,16:23)])

### prepare test sample 
pose3_altered$int <- 1
pose3_alteredL <- dataLongCompRisks(dataShort=pose3_altered, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=F)
### for prediction 
n_test <- 2370
k <- 31
pose3_altered_pred <- pose3_altered[rep(1:n_test, each=k), ]
pose3_altered_pred <- data.frame("timeInt"=rep(1:k, n_test), pose3_altered_pred)

### Fit GB with StabS 

# prepare data and formula
pose12_alteredL_fit <- rbind(pose12_alteredL, pose12_alteredL)
pose12_alteredL_fit$y   <- c(pose12_alteredL$death, pose12_alteredL$discharge)
pose12_alteredL_fit$int <- 1

pose12_alteredL_fit[,c("timeInt","gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")] <- scale(pose12_alteredL_fit[,c("timeInt","gender","age","frailty","plasma","platelets","red_blood","multimorbidity","timed_up_go")], scale=F)
form  <- paste0("bols(", names, ",intercept=F)", collapse="+")
form  <- formula(paste0("y~bols(int, intercept=F)+bols(timeInt, intercept=F)+bbs(timeInt, center=T, df=1)+", form, "+bbs(age, center=T, df=1)"))

# gamboostLSS with mstop=10000
mod <- gamboostLSS(formula = list(eta1 = form,
                                  eta2 = form),
                   families = CR2(stabilization = "none"), data = pose12_alteredL_fit, method = "noncyclic", 
                   control = boost_control(mstop = 10000, trace = T, nu = c(eta1 = 0.1, eta2 = 0.1)))

# apply stability selection with pi=0.7, q=15 and B=300.  
set.seed(11012021)
# stab <- stabsel(mod, cutoff=0.7, q=15, sampling.type="MB", B=300)
# Note: Depending on the system the call of stabsel might take 12-24 hours. The result is stored in "./Application/stab.rda".

# load results 
load("./Application/stab.rda")

#####

# Figure S7. The plot presents the selection frequencies for the 25 most frequently selected base-learners, as obtained from fitting
# the proposed GB model. The vertical line indicates the threshold value pi=0.7.   

par(mar=c(4,18,1,1))
plot(stab, main="", yaxt="n", xlab=bquote(hat(pi)), cex.axis=1.3, cex.lab=1.5, np=25)
axis(2, c(expression(paste("Transfusion of plasma (", discharge, ")")),
          expression(paste("Type of intervention (", discharge, ")")),
          expression(paste("Frailty (", discharge, ")")),
          expression(paste("Urgency (", discharge, ")")),
          expression(paste("Time, linear (", discharge, ")")),
          expression(paste("Time, smooth (", discharge, ")")), 
          expression(paste("Intercept (", discharge, ")")), 
          expression(paste("Transfusion of plasma (", death, ")")),
          expression(paste("Type of intervention (", death, ")")), 
          expression(paste("Frailty (", death, ")")),
          expression(paste("Intercept (", death, ")")), 
          expression(paste("Anaesthesia technique (", death, ")")),
          expression(paste("Severity (", death, ")")), 
          expression(paste("Facility (", discharge, ")")),
          "...",
          "...",
          "..."), at=seq(25,9), las=2, line=0, cex.axis=1.3)

#####

### Fit final model with VGAM 
formf <- formula("cbind(e0, e1, e2)~sm.ps(timeInt)+plasma+category9+frailty+technique+urgency+severity+facility")
clist <- list("(Intercept)" = diag(2),
              "sm.ps(timeInt)" = rbind(0,1),
              "plasma" = diag(2),
              "category9" = diag(2),
              "frailty" = diag(2),
              "technique" = rbind(1,0),
              "urgency" = rbind(0,1),
              "severity"=rbind(1,0),
              "facility" = rbind(0,1))
modf  <- vgam(formf, data=pose12_alteredL, family=multinomial(parallel=FALSE, refLevel=1), constraints=clist)

### Compute confidence intervals 
# Note: Depending on the system the following loop might take some time. The result is stored in "./Application/ci.rda"
res_ci <- array(NA, dim=c(30,2,300))
f_discharge <- matrix(NA, nrow=300, ncol=31)
for(seed in 1:300){
  cat(seed)
  set.seed(seed+18052021)
  okay <- FALSE
  while(!okay){
    ind_boot0 <- sample(which(pose12_altered$status30==0), table(pose12_altered$status30)[1], replace=TRUE)
    ind_boot1 <- sample(which(pose12_altered$status30==1), table(pose12_altered$status30)[2], replace=TRUE)
    ind_boot2 <- sample(which(pose12_altered$status30==2), table(pose12_altered$status30)[3], replace=TRUE)
    ind_boot  <- c(ind_boot0,ind_boot1,ind_boot2)
    data_boot <- pose12_altered[ind_boot,]
    if(all(table(data_boot$category9, data_boot$status30)>0)){
      okay <- TRUE
    }
  }
  data_boot_long <- dataLongCompRisks(dataShort=data_boot, timeColumn="time30", eventColumns=c("death","discharge"), timeAsFactor=F)
  mod_boot       <- vgam(formf, data=data_boot_long, family=multinomial(parallel=FALSE, refLevel=1), constraints=clist)
  res_ci[,,seed] <- coef(mod_boot, matrix=T)
  fid <- data_boot$Study.Subject.ID[data_boot$time30==31][1]
  f_discharge[seed,] <- mod_boot@x[which(data_boot_long$Study.Subject.ID==fid)[1:31],2:10]%*%res_ci[2:10,2,seed]
}

#####

# Table 1. The table presents the coefficient estimates and 95% confidence intervals obtained from fitting a 
# discrete cause-specific hazard model to the learning sample. 

coeff <- coef(modf, matrix=T)
lb <- apply(res_ci, c(1,2), quantile, probs=0.025)
ub <- apply(res_ci, c(1,2), quantile, probs=0.975)
res1 <- cbind(coeff[,1], lb[,1], ub[,1], exp(coeff[,1]))
res2 <- cbind(coeff[,2], lb[,2], ub[,2], exp(coeff[,2]))

round(res1[-c(1:10),], 4) # Death 
round(res2[-c(1:10),], 4) # Discharge

#####

# Figure 3. Time-dependent baseline coefficients for discharge alive obtained from fitting the discrete cause-specific hazard model 
# to the learning sample of the POSE data.

time_discharge <- modf@x[which(pose12_alteredL$Study.Subject.ID=="090-007-062")[1:31],2:10]%*%res2[2:10,1]
time_discharge_l <- apply(f_discharge, 2, quantile, probs=0.025)
time_discharge_u <- apply(f_discharge, 2, quantile, probs=0.975)
par(mar=c(5,5,1,1))
plot(1:30, time_discharge[-31], type="s", ylim=c(-1.2,0.4), lwd=3, xlab="Days", ylab=bquote(f["0,discharge"]),
     cex.axis=1.8, cex.lab=2, col=0)
lines(1:30, time_discharge_l[-31], lty="solid", lwd=2.5, col="darkgrey", type="s")
lines(1:30, time_discharge_u[-31], lty="solid", lwd=2.5, col="darkgrey", type="s")
lines(1:30, time_discharge[-31], type="s", lwd=3)

#####

# Figure 4. The figure shows the estimated cumulative incidence functions for death (panels (a) and (c)) and 
# discharge alive (panels (b) and (d)) referring to the covariates of a randomly selected patient 
# (`elective' and `gynaecologic and urological' intervention, no transfusion of plasma, not frail, 
# general anaesthesia technique, no premedication). Panels (a) and (b) refer to the frailty of the patient, 
# panels (c) and (d) refer to different degrees of urgency of the intervention.

lambda_hat_pred <- predict(modf, type="response", newdata=pose3_altered_pred)
lambda1_pred    <- matrix(lambda_hat_pred[,2], nrow=n_test, ncol=k, byrow=T)
lambda2_pred    <- matrix(lambda_hat_pred[,3], nrow=n_test, ncol=k, byrow=T)

pose_selected <- pose3_altered[,c(11,12,14,17,20,22)]
ids41 <- which(pose_selected$urgency=="elective" & pose_selected$category9=="gyn_uro" & pose_selected$frailty==0 &  
              pose_selected$plasma==0 & pose_selected$technique=="General" & pose_selected$premedication=="None")
ids42 <- which(pose_selected$urgency=="elective" & pose_selected$category9=="gyn_uro" & pose_selected$frailty==1 &  
                 pose_selected$plasma==0 & pose_selected$technique=="General" & pose_selected$premedication=="None")
ids43 <- which(pose_selected$urgency=="urgent" & pose_selected$category9=="gyn_uro" & pose_selected$frailty==0 &  
                 pose_selected$plasma==0 & pose_selected$technique=="General" & pose_selected$premedication=="None")
ids44 <- which(pose_selected$urgency=="emergency" & pose_selected$category9=="gyn_uro" & pose_selected$frailty==0 &  
                 pose_selected$plasma==0 & pose_selected$technique=="General" & pose_selected$premedication=="None")

l_hats <- lambda1_pred+lambda2_pred
S_hat  <- t(apply(1-l_hats, 1, cumprod))
F_hat1 <- t(sapply(1:n_test, function(j) cumsum(c(1,S_hat[j,1:(k-1)])*lambda1_pred[j,])))
F_hat2 <- t(sapply(1:n_test, function(j) cumsum(c(1,S_hat[j,1:(k-1)])*lambda2_pred[j,])))

### Figure 4 a) and b) 
par(mfrow=c(1,2), mar=c(5,5.5,1,1))
plot(0:31, c(0,F_hat1[ids41[1],]), type="s", ylim=c(0,0.035), col=grey(0.8), lwd=2, xlab="Days", ylab=bquote(hat(F)[death]),
     cex.axis=1.5, cex.lab=1.5) # frailty = 0
lines(0:31, c(0,F_hat1[ids42[1],]), type="s", col=grey(0.2), lwd=2) # frailty = 1
legend("topleft", legend=c("frail", "not frail"), lwd=2, col=c(grey(0.2),grey(0.8)), bty="n", cex=1.2)

plot(0:31, c(0,F_hat2[ids42[1],]), type="s", ylim=c(0,1), col=grey(0.2), lwd=2, xlab="Days", ylab=bquote(hat(F)[discharge]),
     cex.axis=1.5, cex.lab=1.5) # frailty = 1
lines(0:31, c(0,F_hat2[ids41[1],]), type="s", col=grey(0.8), lwd=2) # frailty = 0 
legend("bottomright", legend=c("frail", "not frail"), lwd=2, col=c(grey(0.2),grey(0.8)), bty="n", cex=1.2)

### Figure 4 c) and d) 
par(mfrow=c(1,2), mar=c(5,5.5,1,1))
plot(0:31, c(0,F_hat1[ids41[1],]), type="s", ylim=c(0,0.035), col=grey(0.8), lwd=2, xlab="Days", ylab=bquote(hat(F)[death]),
     cex.axis=1.5, cex.lab=1.5) # urgency = elective
lines(0:31, c(0,F_hat1[ids43[1],]), type="s", col=grey(0.5), lwd=2) # urgency = urgent
lines(0:31, c(0,F_hat1[ids44[1],]), type="s", col=grey(0.2), lwd=2) # urgency = emergency
legend("topleft", legend=c("elective", "urgent", "emergency"), lwd=2, col=c(grey(0.8),grey(0.5),grey(0.2)), bty="n", cex=1.2)

plot(0:31, c(0,F_hat2[ids41[1],]), type="s", ylim=c(0,1), col=grey(0.8), lwd=2, xlab="Days", ylab=bquote(hat(F)[discharge]),
     cex.axis=1.5, cex.lab=1.5) # urgency = elective
lines(0:31, c(0,F_hat2[ids43[1],]), type="s", col=grey(0.5), lwd=2) # urgency = urgent
lines(0:31, c(0,F_hat2[ids44[1],]), type="s", col=grey(0.2), lwd=2) # urgency = emergency
legend("bottomright", legend=c("elective", "urgent", "emergency"), lwd=2, col=c(grey(0.8),grey(0.5),grey(0.2)), bty="n", cex=1.2)

#####

# Figure 5. Calibration plots for death and discharge alive using G=22 subsets according the rule 
# suggested by Berger and Schmid (2022). 

lambda_hatf <- predict(modf, type="response", newdata=pose3_alteredL)

par(mfrow=c(1,2), mar=c(5,5,1,1))
### death 
lambda_hat1 <- lambda_hatf[,2]
N     <- nrow(pose3_alteredL)
const <- sqrt(6*(N-2)/(N+1)/(N+3))
nbin  <- 1+floor(log2(N)+log2(1+abs(skewness(lambda_hat1))/const))

quantile_hat <- quantile(lambda_hat1, probs=seq(0, 1, length.out=22))
subset_ID    <- cut(lambda_hat1, breaks=quantile_hat)
subsets      <- split(pose3_alteredL$e1, subset_ID)
obs_hazards1  <- sapply(subsets, mean)
subsets_l    <- split(lambda_hat1, subset_ID)
fit_hazards1  <- sapply(subsets_l, mean)

cal   <- list()
cal$x <- fit_hazards1
cal$y <- obs_hazards1
plot(cal$x, cal$y, xlab="", ylab="",  pch=19, cex=1.2, cex.axis=1.5, ylim=c(0,0.022), xlim=c(0,0.022))
mtext(quote(bar(hat(lambda))["g,death"]), side=1, line=4, cex=1.5)
mtext(bquote(bar(y)["g,death"]), side=2, line=3, las=0, cex=1.5)
lines(c(0,0.03), c(0,0.03), lty="dashed", lwd=3, col="grey")

### discharge alive 
lambda_hat2 <- lambda_hatf[,3]
N     <- nrow(pose3_alteredL)
const <- sqrt(6*(N-2)/(N+1)/(N+3))
nbin  <- 1+floor(log2(N)+log2(1+abs(skewness(lambda_hat2))/const))

quantile_hat <- quantile(lambda_hat2, probs=seq(0, 1, length.out=22))
subset_ID    <- cut(lambda_hat2, breaks=quantile_hat)
subsets      <- split(pose3_alteredL$e2, subset_ID)
obs_hazards2  <- sapply(subsets, mean)
subsets_l    <- split(lambda_hat2, subset_ID)
fit_hazards2  <- sapply(subsets_l, mean)

cal   <- list()
cal$x <- fit_hazards2
cal$y <- obs_hazards2
plot(cal$x, cal$y, xlab="", ylab="",  pch=19, cex=1.2, cex.axis=1.5, ylim=c(0,0.5), xlim=c(0,0.5))
mtext(quote(bar(hat(lambda))["g,discharge"]), side=1, line=4, cex=1.5)
mtext(bquote(bar(y)["g,discharge"]), side=2, line=3, las=0, cex=1.5)
lines(c(0,0.5), c(0,0.5), lty="dashed", lwd=3, col="grey")

#####
