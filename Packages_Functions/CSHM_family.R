##############################################################################
#####         Modeling postoperative mortality in older patients         ##### 
#####          by boosting discrete-time competing risks models          #####
##############################################################################
#####                    Electronic Supplement                           #####
#####                    Author: Moritz Berger                           #####
##############################################################################
#####       Content: Family function CSHM() for fitting the              #####
#####                proposed model using the R package gamboostLSS      #####
##############################################################################

CSHM_lambda1 <- function(eta1, eta2, stabilization){
  
  loss <- function(eta2, y, f, w = 1){
    f  <- pmin(abs(f), 36) * sign(f)
    y1  <- y
    y1[seq(length(y1)/2+1,length(y1), by=1)] <- 0 
    y2  <- y
    y2[seq(1,length(y2)/2, by=1)] <- 0 
    fy <- y1*f+y2*eta2-log(1+exp(f)+exp(eta2))
    return(-fy)
  }
  
  risk <- function(y, f, w = 1){
    sum(w*loss(y = y, f = f, eta2 = eta2))
  }
  
  ngradient <- function(y, f, w = 1){
    f   <- pmin(abs(f), 36) * sign(f)
    l1  <- exp(f)/(1+exp(f)+exp(eta2))
    y1  <- y
    y1[seq(length(y1)/2+1,length(y1), by=1)] <- 0 
    ngr <- y1-l1
    ngr <- gamboostLSS:::stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr) 
  }
  
  offset <- function(y, w = 1) {
    0
  }
  
  response <- function(f){
    f
  }
  
  Family(ngradient = ngradient, risk = risk, offset = offset, response = response, name = "CR_lambda1")
  
}

CSHM_lambda2 <- function(eta1, eta2, stabilization){
  
  loss <- function(eta1, y, f, w = 1){
    f  <- pmin(abs(f), 36) * sign(f)
    y1  <- y
    y1[seq(length(y1)/2+1,length(y1), by=1)] <- 0 
    y2  <- y
    y2[seq(1,length(y2)/2, by=1)] <- 0 
    fy <- y1*eta1+y2*f-log(1+exp(eta1)+exp(f))
    return(-fy)
  }
  
  risk <- function(y, f, w = 1){
    sum(w*loss(y = y, f = f, eta1 = eta1))
  }
  
  ngradient <- function(y, f, w = 1){
    f   <- pmin(abs(f), 36) * sign(f)
    l2  <- exp(f)/(1+exp(eta1)+exp(f))
    y2  <- y
    y2[seq(1,length(y2)/2, by=1)] <- 0 
    ngr <- y2-l2
    ngr <- gamboostLSS:::stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr) 
  }
  
  offset <- function(y, w = 1) {
    0
  }
  
  response <- function(f){
    f
  }
  
  Family(ngradient = ngradient, risk = risk, offset = offset, response = response, name = "CR_lambda2")
  
}

CSHM <- function (eta1 = NULL, eta2 = NULL,
                  stabilization = c("none", "MAD", "L2")){
                  stabilization <- gamboostLSS:::check_stabilization(stabilization)
                  Families(eta1 = CSHM_lambda1(eta1, eta2, stabilization = stabilization),
                           eta2 = CSHM_lambda2(eta1, eta2, stabilization = stabilization))
                  }
