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

# The following contains R-Code to evaluate the results of the additional simulation 
# scenario (c) with nonlinear covariate effects. 

# The R-Code provided can be used to reproduce Figure S5, Figure S6 and 
# Table S5, Table S6 of the supplementary material. 

### (optional) install libraries
# install.packages("ggplot2", "viridis")

### libraries
library(ggplot2)
library(viridis)

### results 
load("./Simulation/raw_results/simulationAc.rda")

### Figure S5 ###
Cs      <- sapply(c(301:600,2101:2400,601:900,1501:1800,1201:1500,1801:2100,901:1200,1:300), function(j) res[[j]][[2]])

args_plot <- expand.grid(seed=c(1:100)+160920, 
                         n=250,
                         k=10,
                         b=0.9,
                         noninf=c(10,200,1000),
                         mmax=15000, 
                         method=c(0:7))
noninf <- c("10"="low-dimensional","200"="moderate-dimensional","1000"="high-dimensional")

for_plot <- data.frame(values=Cs, Approach=factor(args_plot$method), noninf=args_plot$noninf)
medCs <- aggregate(values~noninf+Approach, median, data=for_plot[!for_plot$Approach==7,])
for_plot$hlines <- rep(rep(aggregate(values~noninf, max, data=medCs)$values, each=100),8)

ggplot(for_plot, aes(x=Approach,y=values))+
  geom_jitter(mapping = aes(color=Approach), position=position_jitter(0.15))+
  geom_boxplot(mapping = aes(fill=Approach), notch=TRUE, outlier.shape=NA)+
  facet_grid(~noninf,labeller=labeller(noninf=noninf))+
  scale_fill_viridis(discrete=TRUE, labels=c("GB logL","GB StabS", "GBg", "PL","PLg","Tree","Null","True"), option="D")+
  scale_color_viridis(discrete=TRUE, labels=c("GB logL","GB StabS","GBg", "PL","PLg","Tree","Null","True"), option="D")+
  xlab("Approach")+
  ylab("C-index")+
  ylim(0.4,0.9)+
  scale_x_discrete(labels=c("GB \n logL ","GB \n StabS","GBg", "PL","PLg","Tree","Null","True"))+
  geom_hline(aes(yintercept = hlines), size=1, lty="dashed")+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size= 18),
        legend.title = element_text(size=24, face="bold"),
        legend.text  = element_text(size=18),
        legend.position = "top",
        legend.justification = "left",
        strip.text.x = element_text(size= 18),
        strip.text.y = element_text(size= 18))


### Table S5 ### 
IPEs     <- sapply(c(301:600,2101:2400, 601:900,1501:1800,1201:1500,1801:2100,901:1200,1:300), function(j) res[[j]][[4]])

for_table <- data.frame(values=IPEs, method=args_plot$method, noninf=args_plot$noninf)
mean_S5 <- aggregate(values~noninf+method, mean, data=for_table)
mean_S5[,3]*100
sd_S5 <- aggregate(values~noninf+method, sd, data=for_table)
sd_S5[,3]*100

### Figure S6 ### 
PEs1 <- c(sapply(c(1:300,901:1200,1801:2100,1201:1500,1501:1800,601:900,2101:2400,301:600), function(j) res[[j]][[3]][,1]))
PEs2 <- c(sapply(c(1:300,901:1200,1801:2100,1201:1500,1501:1800,601:900,2101:2400,301:600), function(j) res[[j]][[3]][,2]))

for_plot1 <- data.frame(values=PEs1, time=factor(rep(1:10, 2400)), method=factor(rep(args_plot$method, each=10)), 
                        noninf=rep(args_plot$noninf, each=10))
for_plot1 <- aggregate(values~time+method+noninf, mean, data=for_plot1)
for_plot1$values[for_plot1$time==10] <- for_plot1$values[for_plot1$time==9]

for_plot2 <- data.frame(values=PEs2, time=factor(rep(1:10, 2400)), method=factor(rep(args_plot$method, each=10)), 
                        noninf=rep(args_plot$noninf, each=10))
for_plot2 <- aggregate(values~time+method+noninf, mean, data=for_plot2)
for_plot2$values[for_plot2$time==10] <- for_plot2$values[for_plot2$time==9]

for_plot <- rbind(for_plot1,for_plot2)
for_plot$event <- rep(c(1,2), each=240)
event  <- c("1" = "event 1", "2" = "event 2")

ggplot(for_plot, aes(x=time,y=values, color=method, group=method))+
  geom_point()+geom_step(size=0.8)+
  scale_color_viridis(discrete=TRUE, labels=c("True", "Null","Tree","PLg","PL","GBg","GB StabS","GB logL"), option="D", 
                      direction=-1, alpha=0.75)+
  facet_grid(event~noninf,labeller=labeller(event=event, noninf=noninf))+
  xlab("Time")+
  ylab(bquote(PE[j](t)))+
  ylim(0.075,0.25)+
  labs(color="Approach")+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size= 18),
        legend.title = element_text(size=18, face="bold"),
        legend.text  = element_text(size=15),
        legend.position = "top",
        legend.justification = "left",
        strip.text.x = element_text(size= 18),
        strip.text.y = element_text(size= 18))


### Table S6 ### 
coefs7 <- lapply(2101:2400, function(j) res[[j]][[5]])
vars7  <- lapply(1:300, function(j) gsub("(.*)(x[0-9]{1,}s?)(.*)", "\\2", rownames(coefs7[[j]])))

get_effects7 <- function(index, var, eta){
  li <- sapply(index, function(j) (var%in%vars7[[j]]) && (coefs7[[j]][var,eta] !=0))*1
  sm <- sapply(index, function(j) (paste0(var,"s") %in% vars7[[j]]) && any(coefs7[[j]][which(vars7[[j]]==paste0(var,"s")),eta]!=0))*1
  return(c(sum(!li & !sm),sum(li & !sm), sum(sm)))
}

get_effects7(1:100, "x1", 1)
get_effects7(101:200, "x1", 1)
get_effects7(201:300, "x1", 1)

get_effects7(1:100, "x1", 2)
get_effects7(101:200, "x1", 2)
get_effects7(201:300, "x1", 2)

get_effects7(1:100, "x2", 1)
get_effects7(101:200, "x2", 1)
get_effects7(201:300, "x2", 1)

get_effects7(1:100, "x2", 2)
get_effects7(101:200, "x2", 2)
get_effects7(201:300, "x2", 2)

get_effects7(1:100, "x3", 1)
get_effects7(101:200, "x3", 1)
get_effects7(201:300, "x3", 1)

get_effects7(1:100, "x3", 2)
get_effects7(101:200, "x3", 2)
get_effects7(201:300, "x3", 2)

get_effects7(1:100, "x4", 1)
get_effects7(101:200, "x4", 1)
get_effects7(201:300, "x4", 1)

get_effects7(1:100, "x4", 2)
get_effects7(101:200, "x4", 2)
get_effects7(201:300, "x4", 2)

get_effects7(1:100, "x5", 1)
get_effects7(101:200, "x5", 1)
get_effects7(201:300, "x5", 1)

get_effects7(1:100, "x5", 2)
get_effects7(101:200, "x5", 2)
get_effects7(201:300, "x5", 2)

get_effects7(1:100, "x6", 1)
get_effects7(101:200, "x6", 1)
get_effects7(201:300, "x6", 1)

get_effects7(1:100, "x6", 2)
get_effects7(101:200, "x6", 2)
get_effects7(201:300, "x6", 2)

get_effects7(1:100, "x7", 1)
get_effects7(101:200, "x7", 1)
get_effects7(201:300, "x7", 1)

get_effects7(1:100, "x7", 2)
get_effects7(101:200, "x7", 2)
get_effects7(201:300, "x7", 2)

get_effects7(1:100, "x8", 1)
get_effects7(101:200, "x8", 1)
get_effects7(201:300, "x8", 1)

get_effects7(1:100, "x8", 2)
get_effects7(101:200, "x8", 2)
get_effects7(201:300, "x8", 2)
