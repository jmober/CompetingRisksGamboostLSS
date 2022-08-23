##############################################################################
##### A CLASSIFICATION TREE APPROACH FOR THE MODELING OF COMPETING RISKS ##### 
#####                          IN DISCRETE TIME                          #####
#####                        Electronic Supplement                       #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Function to fit a competing risks discrete survival tree for a sequence of 
# minimal node sizes. 

### Input: 
# formula: description of the model, like for lm() or glm()
# data: data in short format 
# distance: splitting criterion Hellinger distance or Gini impurity 
# eventColumns: for dataLong() in discSurv 
# minimal_ns: sequence of minimal node sizes (tuning parameters)

### Description: 
# The function requires the R add-on package discSurv and sampling, and the self-implemented 
# function CRTreeDisc_fit.R 

### Output: 
# est: list of fitted tree as returned by CRTreeDisc_fit()


CRTreeDisc <- function(formula, data, 
                       distance=c("Hellinger","Gini"),
                       eventColumns, minimal_ns, 
                       trace = FALSE){
  


  require(discSurv)
  
  # adjust input 
  form       <- paste(formula)
  timeColumn <- form[2]
  form       <- formula(paste("y~timeInt+", form[3]))
  
  distance <- match.arg(distance)

  nevents       <- length(eventColumns)
  dataShortAll  <- data
  dataLongAll   <- dataLongCompRisks(data, timeColumn, eventColumns,
                                    timeAsFactor = FALSE)
  dataLongAll$y <- apply(dataLongAll[,3:(3+nevents)],1,function(x) which(x==1)-1)
  dataLongAll$y <- factor(dataLongAll$y)
  
  if(missing(minimal_ns)){
    minimal_ns <- 1:(floor(nrow(dataLongAll) / 2))
  }
  
  # One Tree
  one_tree <- function(mb, data){
    
    y     <- data$y
    X     <- data[,all.vars(formula[[3]])]
    model <- CRTreeDisc_fit(X,y,mb,distance)
    
    p       <- CRTreeDisc_predict(model,X)$prob
    yaug    <- data[,3:(3+nevents)]
    ll      <- sum(yaug*log(p))
    nsplits <- sum(grepl("terminal",names(unlist(model))))-1

    AIC <- - 2 * ll + 2 * nsplits
    BIC <- - 2 * ll + log(nrow(data)) * nsplits
    
    return(list("model" = model, "AIC" = AIC, "BIC" = BIC, "p" = p))
  }
  
  est <- list()
  
  k <- 0 
  for(MB in minimal_ns){
    k   <- k + 1
    est[[k]] <- one_tree(MB, dataLongAll)$model
    if(trace & k == 1){
      cat("minimal node size:\n")
    }
    if(trace){
      cat(MB, "\n")
    }
  }
  
  return(est)
}