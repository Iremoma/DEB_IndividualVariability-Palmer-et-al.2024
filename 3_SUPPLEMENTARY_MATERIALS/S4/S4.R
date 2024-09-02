#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 4
#---------------------------------------
# Last update: 30 January 2024

#-----
# 1: Loading libraries, data and results
#-----
remove(list=ls())
load("input.RData")  

library(cmdstanr) 
library(bayesplot) 
library(posterior) 
library(ggplot2)  
library(deSolve)  
library(gridExtra)
library(ggplot2)
library(cowplot)


# Loading results with crowding:

fit = readRDS("out.RDS")  
draws=fit$draws(format="matrix")

#-----
# Figure S5: Expected effects of crowding on assimilation rate. The correction 
# factor C_B ([0,1]) is plotted against biomass (left panel) and fish age (right panel).
#----- 

t=seq(t0,age[n])
W0=27652.82
W1=median(data.frame(as_draws_df(fit$draws(paste("W1"))))[,1])
biomass = 3.173425e+03-3.720773e+02*t+2.362952e+00*t^2-1.431395e-03*t^3;
W = 1*(biomass<=W0)+
  (1+(-W0/(W0-W1))+(biomass/(W0-W1)))*(biomass>W0)*(biomass<W1)+
  0*(biomass>=W1);

data = data.frame(t = t, W = W)

# Biomass vs. Correction factor:

A= ggplot(data, aes(x = biomass / 1000 / 20, y = W)) +
  geom_line(color = "#333399", size = 1) +
  labs(x = bquote("Biomass (Kg "*m^-3*")"), y = "Correction factor") +
  theme(
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))


# Fish age vs. Correction factor:

B= ggplot(data, aes(x = t, y = W)) +
  geom_line(color = "#333399", size = 1) +
  labs(x = "Fish age (days)", y = "Correction factor") +
  theme(
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))



ggdraw (
  plot_grid (A,
             B,
             ncol=2,align='h'))

# Saving: 

# ggsave(filename=paste("S5.png"),path=paste0(getwd(),"/Out",sep=""),
#        width=12,
#        height=4,
#        units="cm")



#-----
# 2: Preparing data for Figure S4 and S6
#-----

# WITH CROWDING:

temp=array(NA,dim=c(N,2,(n+1)))
residuals=temp
for (i in 1:N){
  print(i)
  temp[i,1,1]=median(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",1]",sep=""))))[,1])
  temp[i,2,1]=median(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",2]",sep=""))))[,1])
  residuals[i,1,1]=(obs0[i,1]-temp[i,1,1])/temp[i,1,1]
  residuals[i,2,1]=(obs0[i,2]-temp[i,2,1])/temp[i,2,1]
  for (j in 1:n){
    temp[i,1,(j+1)]=median(data.frame(as_draws_df(fit$draws(paste("obs_hat[",i,",1,",j,"]",sep=""))))[,1])
    temp[i,2,(j+1)]=median(data.frame(as_draws_df(fit$draws(paste("obs_hat[",i,",2,",j,"]",sep=""))))[,1])
    residuals[i,1,(j+1)]=(obs[i,1,j]-temp[i,1,(j+1)])/temp[i,1,(j+1)]
    residuals[i,2,(j+1)]=(obs[i,2,j]-temp[i,2,(j+1)])/temp[i,2,(j+1)]
  }
}
# Observed weight:
temp3=array(NA,dim=c(3,(n+1)))
temp3[,2:(n+1)]=as.array(apply(obs[,2,],2,quantile,c(0.025,0.5,0.975)))
temp3[,1]=quantile(obs0[,2],c(0.025,0.5,0.975))
temp2=c(t0,age)
# Expected weight:
temp4=apply(temp[,2,],2,quantile,c(0.025,0.5,0.975))


# WITHOUT CROWDING:  

fit = readRDS("no_crowding.RDS")  
draws=fit$draws(format="matrix")

# Expected weight:
temp=array(NA,dim=c(N,2,(n+1)))
residuals_no=temp
for (i in 1:N){ 
  print(i)
  temp[i,1,1]=median(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",1]",sep=""))))[,1])
  temp[i,2,1]=median(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",2]",sep=""))))[,1])
  residuals_no[i,1,1]=(obs0[i,1]-temp[i,1,1])/temp[i,1,1]
  residuals_no[i,2,1]=(obs0[i,2]-temp[i,2,1])/temp[i,2,1]
    for (j in 1:n){
    temp[i,1,(j+1)]=median(data.frame(as_draws_df(fit$draws(paste("obs_hat[",i,",1,",j,"]",sep=""))))[,1])
    temp[i,2,(j+1)]=median(data.frame(as_draws_df(fit$draws(paste("obs_hat[",i,",2,",j,"]",sep=""))))[,1])
    residuals_no[i,1,(j+1)]=(obs[i,1,j]-temp[i,1,(j+1)])/temp[i,1,(j+1)]
    residuals_no[i,2,(j+1)]=(obs[i,2,j]-temp[i,2,(j+1)])/temp[i,2,(j+1)]
  }
}

temp5=apply(temp[,2,],2,quantile,c(0.025,0.5,0.975))



#--------------
# Figure S4: Expected weight-at-age predicted by a DEB model that (1) that 
# ignores the crowding effects (red line) and (2) includes the effects of 
# crowding and predicts growth deceleration toward the end of the experiment 
# (i.e., the DEB version referred in the main text; blue line). The observed 
# weight-at-age at each sampling event are denoted by dots (median) and black 
# lines (95% quantile).
#--------------


data = data.frame(Age = temp2, WetWeight = temp3[2,])

p=ggplot(data, aes(x = Age, y = WetWeight)) +
  geom_point(aes(color = "Observed ",shape="Observed"), size = 1.2) +
  geom_line(aes(x = temp2, y = temp4[2,], color = "With crowding",  alpha = 1)) +
  geom_line(aes(x = temp2, y = temp5[2,], color = "Without crowding", alpha = 1)) +
  labs(x = "Age (days)", y = "Wet weight (gr)") +
  ylim(range(temp3)) +
  scale_color_manual(
    labels = c("Observed ","With crowding", "Without crowding"),
    values = c("darkgrey","#333399", "red"))+
  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.ticks.length.y = unit (.25, "cm"), 
      axis.ticks.length.x = unit (.25, "cm"),
      panel.border = element_rect(linetype = 1, fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent")
    )
  
p=p+guides(alpha = FALSE, shape = FALSE,size=FALSE)
p=p + theme(legend.position = c(0.9, 0.15))
p=p + theme(legend.title = element_blank())
p

# CI for observed data
# for (i in 1:(n + 1)) {
#   p = p + geom_line(
#     data = data.frame(Age = rep(temp2[i], 2), WetWeight = temp3[c(1, 3), i]),
#     aes(x = Age, y = WetWeight), color = "darkgrey"
#   )
# }
# p

#Saving: 

# ggsave(filename=paste("S4.png"),path=paste0(getwd(),"/Out",sep=""),width = 8, height = 5)



#-----------
# Figure S6: Relative residuals at the level of sampling event. Left panels: 
# DEB model versions with crowding; Right panels: without crowding. Upper panels: 
# total length; Lower panels: wet weight. The boxplots denote the distribution 
# of the residuals at the fish level. The colors indicate the different length 
# measurement devices (see main text). 
#-----------

z=1 # (including initial sampling #0)

# Length with crowding:

data_with_crowding = data.frame(
  SamplingEvent = factor(rep(0:10, each = nrow(residuals))),
  Length = c(residuals[, 1, z:11]),
  length_method = factor(rep(1 + c(1,length_method), each = nrow(residuals))),
  Crowding = "With Crowding"
)
colors = c("#9999FF", "#333399")

lengthwc=ggplot(data_with_crowding, aes(x = SamplingEvent, y = Length, fill = length_method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0)+
  labs(title = "Length With Crowding", x = "Sampling Event", y = "") +
  scale_fill_manual(values = colors,guide = FALSE) + 
  ylim(range(c(residuals[, 1, z:11], residuals_no[, 1, z:11]))) +
  guides(fill = "none")+
  theme(
    axis.title.x=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))

# Length with no crowding:

data_with_no_crowding = data.frame(
  SamplingEvent = factor(rep(0:10, each = nrow(residuals_no))),
  Length = c(residuals_no[, 1, z:11]),
  length_method = factor(rep(1 + c(1,length_method), each = nrow(residuals))),
  Crowding = "Without Crowding")

colors = c("#9999FF", "#333399")

lengthno=ggplot(data_with_no_crowding, aes(x = SamplingEvent, y = Length, fill = length_method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0)+
  labs(title = "Length With No Crowding", x = "Sampling Event", y = "") +
  scale_fill_manual(values = colors,guide = FALSE) +  
  ylim(range(c(residuals[, 1, z:11], residuals_no[, 1, z:11]))) +
  guides(fill = "none")+
  theme(
    axis.title.x=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))

# Weight with crowding:

data_with_crowding = data.frame(
  SamplingEvent = factor(rep(0:10, each = nrow(residuals))),
  Weight = c(residuals[, 2, z:11]),
  Crowding = "With Crowding")


weightwc=ggplot(data_with_crowding, aes(x = SamplingEvent, y = Weight)) +
  geom_boxplot(outlier.shape = NA,fill="grey") +
  geom_hline(yintercept = 0)+
  labs(title = "Weight With Crowding", x = "Sampling Event", y = "") +
  ylim(range(c(residuals[, 2, z:11], residuals_no[, 2, z:11]))) +
  guides(fill = "none")+
  theme(
    axis.title.x=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))

# Weight with no crowding:

data_with_no_crowding = data.frame(
  SamplingEvent = factor(rep(0:10, each = nrow(residuals_no))),
  Weight = c(residuals_no[, 2, z:11]),
  Crowding = "With Crowding")


weightno=ggplot(data_with_no_crowding, aes(x = SamplingEvent, y = Weight)) +
  geom_boxplot(outlier.shape = NA,fill="grey") +
  geom_hline(yintercept = 0)+
  labs(title = "Weight With No Crowding", x = "Sampling Event", y = "") +
  ylim(range(c(residuals[, 2, z:11], residuals_no[, 2, z:11]))) +
  guides(fill = "none")+
  theme(
    axis.title.x=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))


ggdraw (
  plot_grid (lengthwc,
             lengthno,
             weightwc,
             weightno,
             ncol=2,nrow=2,align='h'))


# Saving:

# ggsave(filename=paste("S6.png"),path=paste0(getwd(),"/Out",sep=""),
#        width=25,
#        height=15,
#        units="cm")

