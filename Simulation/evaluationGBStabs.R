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

# The following contains R-Code to evaluate the results of the first part of the
# simulation study ('Comparison of variable selection strategies'). 

# The R-Code provided can be used to reproduce Table S1 
# of the supplementary material. 


### results GB logL
load("./Simulation/raw_results/simulation.rda")
reslogL <- res[301:600]

### results GB StabS
load("./Simulation/raw_results/simulationGBStabs.rda")

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
aggregate(IPE_logL~noninf, mean, data=data.frame(IPE_logL*100,noninf))
aggregate(IPE_logL~noninf, sd, data=data.frame(IPE_logL*100,noninf))

# StabS 0.6
IPE_06 <- sapply(1:300, function(j) res[[j]][[4]])
aggregate(IPE_06~noninf, mean, data=data.frame(IPE_06*100,noninf))
aggregate(IPE_06~noninf, sd, data=data.frame(IPE_06*100,noninf))

# StabS 0.7
IPE_07 <- sapply(1:300, function(j) res[[j]][[9]])
aggregate(IPE_07~noninf, mean, data=data.frame(IPE_07*100,noninf))
aggregate(IPE_07~noninf, sd, data=data.frame(IPE_07*100,noninf))

# StabS 0.8
IPE_08 <- sapply(1:300, function(j) res[[j]][[14]])
aggregate(IPE_08~noninf, mean, data=data.frame(IPE_08*100,noninf))
aggregate(IPE_08~noninf, sd, data=data.frame(IPE_08*100,noninf))

# StabS 0.9
IPE_09 <- sapply(1:300, function(j) res[[j]][[19]])
aggregate(IPE_09~noninf, mean, data=data.frame(IPE_09*100,noninf))
aggregate(IPE_09~noninf, sd, data=data.frame(IPE_09*100,noninf))
