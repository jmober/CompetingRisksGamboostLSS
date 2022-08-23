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

# The following contains R-Code to evaluate the second part of the 
# analysis of the POSE data ('Comparison to alternative methods').

# The R-Code provided can be used to reproduce Table S8 of the 
# supplementary material. 

# Note: The results stored in "./raw_results/Resampling.rda" correspond to the 
# results obtained from the original data used in the paper.

### results
load("./Application/Resampling.rda")
library("ggplot2")

args <- expand.grid(seed=c(1:100)+160920, 
                    method=c(1:6))

### Table S8 (C-index)
Cs <- sapply(1:600, function(j) res[[j]][[1]])
for_table <- data.frame(values=c(Cs), Approach=factor(args$method))
Css <- aggregate(values~Approach, mean, data=for_table)
round(Css[,-1], 3)
CssSD <- aggregate(values~Approach, sd, data=for_table)
round(CssSD[,-1],3)

### Table S8 (Recalibration model)
## Coefficients
Coefs <- t(sapply(1:600, function(j) res[[j]][[4]]))
for_table <- data.frame(values=c(Coefs), Approach=rep(factor(args$method),4),Coef=rep(1:4, each=600))
Coefss <- aggregate(values~Approach+Coef, mean, data=for_table)
round(Coefss[,-1], 3)

## lower bound 
lCI <- t(sapply(1:600, function(j) res[[j]][[5]][,1]))
for_table <- data.frame(values=c(lCI), Approach=rep(factor(args$method),4),Coef=rep(1:4, each=600))
lCIs <- aggregate(values~Approach+Coef, mean, data=for_table)
round(lCIs[,-1], 3)

## upper bound 
uCI <- t(sapply(1:600, function(j) res[[j]][[5]][,2]))
for_table <- data.frame(values=c(uCI), Approach=rep(factor(args$method),4),Coef=rep(1:4, each=600))
uCIs <- aggregate(values~Approach+Coef, mean, data=for_table)
round(uCIs[,-1], 3)

## P-values 
PValues <- sapply(1:600, function(j) res[[j]][[6]])
for_table <- data.frame(values=PValues, tvalues=-log10(PValues), Approach=factor(args$method))
Pss <- aggregate(values~Approach, mean, data=for_table)
round(Pss[,-1], 3)
