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
# scenario (a) with 60% censoring. 

# The R-Code provided can be used to reproduce Figure S1, Figure S2 and 
# Table S3 of the supplementary material. 


### libraries
library(ggplot2)
library(viridis)

### results 
load("./Simulation/raw_results/simulationAa.rda")

### Figure S1 ###
Cs      <- sapply(c(301:600,2101:2400,601:900,1501:1800,1201:1500,1801:2100,901:1200,1:300), function(j) res[[j]][[2]])

args_plot <- expand.grid(seed=c(1:100)+160920, 
                         n=250,
                         k=10,
                         b=1.8,
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

### Table S3 ### 
IPEs     <- sapply(c(301:600,2101:2400, 601:900,1501:1800,1201:1500,1801:2100,901:1200,1:300), function(j) res[[j]][[4]])

for_table <- data.frame(values=IPEs, method=args_plot$method, noninf=args_plot$noninf)
mean_S3 <- aggregate(values~noninf+method, mean, data=for_table)
mean_S3[,3]*100
sd_S3 <- aggregate(values~noninf+method, sd, data=for_table)
sd_S3[,3]*100

### Figure S2 ### 
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
  ylim(0.08,0.25)+
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





