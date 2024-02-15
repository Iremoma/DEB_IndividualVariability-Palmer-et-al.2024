#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# Producing table for main test: Table 3.
#
# Producing table for supplementary materials: Table S4.
#
# Producing figures for main test: Figure 1, Figure 2, Figure 3 & Figure 4.
#
# Last update: 25 January 2024
#

#-----
# Loading data, results and libraries
#-----
remove(list=ls())
library(bayesplot) 
library(posterior) 
library(ggplot2)   
library(deSolve)   

load("input.RData")  

fit = readRDS("out.RDS")  # Results from Analysis.R script
draws=fit$draws(format="matrix")

#-----
# Producing Table 3: Summary statistics for model variables derived from 
# Bayesian analysis. For each variable, the table displays the median, 
# 5th and 95th percentiles, Gelman-Rubin convergence diagnostic (Rhat), 
# and effective sample size (ESS_bulk).
#-----

table=fit$summary(c("pAm_mu","pAm_sd","v_mu","v_sd","sd_length","sd_weight","E0_mu","E0_sd","V0_mu","V0_sd","W1"))
table=data.frame(table)[,c(1,6,3,7,8,9)]
table[,c(2,3,4,5)]=round(table[,c(2,3,4,5)],5)
table[,6]=round(table[,6],0)
table[,2:5]=round(table[,2:5],4)
table
#write.table(table,"Table 3.txt")

#-----
# Producing Supplementary materials S2_Table S4: Fish-specific values for all the 
#parameters estimated at the individual level. Individuals are ordered by 
#increasing value of the estimated "f*pAm". Note that the estimates here refer 
#to the q50 of the posterior interval. Colors gradient indicates from lower
#to higher values for "f*pAm" in red and v in green.
#-----
table=array(NA,dim=c(N,4))
table=data.frame(table)
names(table)=c("pAm","v","E0","V0")
temp = apply(data.frame(as_draws_df(fit$draws("pAm")))[,1:N],2,median)
table[,1]=round(temp,2)
temp = apply(data.frame(as_draws_df(fit$draws("v")))[,1:N],2,median)
table[,2]=round(temp,5)
temp = apply(data.frame(as_draws_df(fit$draws("y0")))[,1:N],2,median)
table[,3]=round(temp,0)
temp = apply(data.frame(as_draws_df(fit$draws("y0")))[,(N+1):(2*N)],2,median)
table[,4]=round(temp,3)
table$name=ID.label
table
#write.table(table,"table.txt")

#------------------------------------------------------------------------------
# Producing Figure 1: : Observed length- and weight-at-age for 69 males or Sparus aurata,
# coming from the same breeding family and breed in a common garden. Each line 
# corresponds to the observations for each individual fish (11 sampling events). 
# The different length measurement methods and their rounding are denoted by 
# different colors.
#------------------------------------------------------------------------------

par(mfrow=c(1,2))
plot(c(t0,age),c(obs0[1,1],obs[1,1,]),
     ylim=range(c(obs0[,1],obs[,1,])),
     type="n",xlab="Age (days)",ylab="Total length (cm)")
for (i in 1:N){
  lines(c(t0,age),c(obs0[i,1],obs[i,1,]))
  points(c(t0,age),c(obs0[i,1],obs[i,1,]),pch=19,cex=0.6,
         col=c("red","red","blue","blue","blue","blue","green","green","green","green","green"))
}
legend("bottomright",legend=c("ictiometer 0.5 cm","ictiometer 0.1 cm", "digital device"),
       col=c("red","blue","green"),pch=19,
       cex=0.6)

plot(c(t0,age),c(obs0[1,2],obs[1,2,]),
     ylim=range(c(obs0[,2],obs[,2,])),
     #xlim=c(0,1100),
     type="n",xlab="Age (days)",ylab="Wet weight (gr)")
for (i in 1:N){
  lines(c(t0,age),c(obs0[i,2],obs[i,2,]))
  points(c(t0,age),c(obs0[i,2],obs[i,2,]),pch=19,cex=0.6)
}

#------------------------------------------------------------------------------
# Producing Figure 2: Expected growth pattern for a fish with average values of  f*pAm and
# v (blue line), i.e., the estimated population-level mean values 
# (pAm_mu =  380.14 J/day/cm2, v_mu = 0.335 cm/day); compared with the 
# observed length- and weight-at-age (where the medians of the sampling events 
# are connected by a red line). Vertical lines represent the 95% interval of 
# each of the 11 sampling events for the estimated (red) and the observed (blue) 
# values. 
#------------------------------------------------------------------------------

par(mfrow=c(1,2))
temp=array(NA,dim=c(N,2,(n+1)))
for (i in 1:N){ 
  print(i)
  temp[i,1,1]=median(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",1]",sep=""))))[,1])
  temp[i,2,1]=median(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",2]",sep=""))))[,1])
  for (j in 1:n){
    temp[i,1,(j+1)]=median(data.frame(as_draws_df(fit$draws(paste("obs_hat[",i,",1,",j,"]",sep=""))))[,1])
    temp[i,2,(j+1)]=median(data.frame(as_draws_df(fit$draws(paste("obs_hat[",i,",2,",j,"]",sep=""))))[,1])
  }
}

#--- Observed length

temp3=array(NA,dim=c(3,(n+1)))
temp3[,2:(n+1)]=as.array(apply(obs[,1,],2,quantile,c(0.025,0.5,0.975)))
temp3[,1]=quantile(obs0[,1],c(0.025,0.5,0.975))
temp2=c(t0,age)
plot(temp2,temp3[2,],
     xlab="Age (days)",ylab="Total Length (cm)",
     ylim=range(temp3),
     pch=19,cex=0.8)
lines(temp2,temp3[2,])
for (i in 1:(n+1)){
  lines(rep(temp2[i],2),temp3[c(1,3),i])
}
#--- Expected length

dt=2
temp4=apply(temp[,1,],2,quantile,c(0.025,0.5,0.975))
points(temp2+dt,temp4[2,],
       pch=19,cex=0.8,col="red")
lines(temp2+dt,temp4[2,],col="red")
for (i in 1:(n+1)){
  lines(rep(temp2[i],2)+dt,temp4[c(1,3),i],col="red")
} 
legend("bottomright",c("observed","estimated"),pch=19,cex=0.8,col=c("black","red"))

#--- Observed weight

temp3=array(NA,dim=c(3,(n+1)))
temp3[,2:(n+1)]=as.array(apply(obs[,2,],2,quantile,c(0.025,0.5,0.975)))
temp3[,1]=quantile(obs0[,2],c(0.025,0.5,0.975))
temp2=c(t0,age)
plot(temp2,temp3[2,],
     xlab="Age (days)",ylab="Wet weight (gr)",
     ylim=range(temp3),
     pch=19,cex=0.8)
lines(temp2,temp3[2,])
for (i in 1:(n+1)){
  lines(rep(temp2[i],2),temp3[c(1,3),i])
}
#--- Expected weight

#dt=5
temp4=apply(temp[,2,],2,quantile,c(0.025,0.5,0.975))
points(temp2+dt,temp4[2,],
       pch=19,cex=0.8,col="red")
lines(temp2+dt,temp4[2,],col="red")
for (i in 1:(n+1)){
  lines(rep(temp2[i],2)+dt,temp4[c(1,3),i],col="red")
} 
legend("bottomright",c("observed","estimated"),pch=19,cex=0.8,col=c("black","red"))

#------------------------------------------------------------------------------
# Producing Figure 3: Estimated values for f*pAm and v at the individual level.
# Each boxplot represents a fish (median and 95% Bayesian credibility interval).
# Fish are sorted in decreasing order for f*pAm in the first two panels. Predicted 
# between-individual correlation for the mean estimates for the DEB parameters 
# f*pAm and v is depicted at the third panel.
#------------------------------------------------------------------------------


par(mfrow = c(3,1)) 

pAm.hat = data.frame(as_draws_df(fit$draws("pAm")))[,1:N]
order=order(apply(pAm.hat,2,median))
par(oma = c(1, 1, 1, 1))
par(mar = c(5, 4, 0, 0))
boxplot(pAm.hat[,order], outline=F,xaxt="n",ylab="f*pAm (j/d/cm^2)")
axis(1,c(1:69),ID.label,cex.axis=0.6,las=2)
pAm.hat=apply(pAm.hat,2,median)

v.hat = data.frame(as_draws_df(fit$draws("v")))[,1:N]
boxplot(v.hat[,order], outline=F,xaxt="n",ylab="v (cm/d)")
axis(1,c(1:69),ID.label,cex.axis=0.6,las=2)
v.hat=apply(v.hat,2,median)
v.mu = data.frame(as_draws_df(fit$draws("v_mu")))[,1]

plot(pAm.hat,v.hat,pch=19,xlab="f*pAm (j/d/cm^2)",ylab="v (cm/d)")

# correaltion f*pAm vs v:
cor(pAm.hat,v.hat) 
#[1] -0.5719766 

#------------------------------------------------------------------------------
# Producing Figure 4: Observed (points) and expected (blue lines) values for 
# length- and weight-at-age predicted by DEB for 4 selected individual fish (following 
# criteria based on figure S16). The solid line is the median expectation. 
# The dashed lines denote 95% Bayesian credibility interval (this CI does not 
# include the measurement error). Note that the observed patterns result from 
# fish-specific combinations of f*pAm and v. 
#------------------------------------------------------------------------------

pAm.hat = data.frame(as_draws_df(fit$draws("pAm")))[,1:N]
pAm.hat=apply(pAm.hat,2,median)
v.hat = data.frame(as_draws_df(fit$draws("v")))[,1:N]
v.hat=apply(v.hat,2,median)
residuals=array(NA,dim=c(N,2,n+1))

layout(matrix(c(1,2,0,3,4,0,0,0,0,0,5,6,0,7,8),3,5,byrow=T),
       heights=c(4,1,4), widths=c(4,4,2,4,4))
par(oma = c(1, 1, 0, 0))
par(mar = c(1, 1, 2, 1))
#for (i in 1:N){ # all 69 fish
#    print(i)
for (i in c(24,47,69,32)){  # 4 selected fish only
  length_hat=array(NA,c(n,3))
  weight_hat=array(NA,c(n,3))
  for (j in 1:n){
    temp=as_draws_df(fit$draws(paste("obs_hat[",i,",1,",j,"]",sep="")))
    temp=data.frame(temp)[,1]
    length_hat[j,]=quantile(temp,c(0.025,0.5,0.975))
    temp=as_draws_df(fit$draws(paste("obs_hat[",i,",2,",j,"]",sep="")))
    temp=data.frame(temp)[,1]
    weight_hat[j,]=quantile(temp,c(0.025,0.5,0.975))
  }
  length0_hat=quantile(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",1]",sep=""))))[,1],c(0.025,0.5,0.975))
  weight0_hat=quantile(data.frame(as_draws_df(fit$draws(paste("obs0_hat[",i,",2]",sep=""))))[,1],c(0.025,0.5,0.975))
  
  #residuals
  residuals[i,1,1]=(length0_hat[2]-obs0[i,1])/length0_hat[2]
  residuals[i,2,1]=(weight0_hat[2]-obs0[i,2])/weight0_hat[2]
  for (j in 1:n){
    residuals[i,1,1+j]=(length_hat[j,2]-obs[i,1,j])/length_hat[j,2]
    residuals[i,2,1+j]=(weight_hat[j,2]-obs[i,2,j])/weight_hat[j,2]
  }
  #length
  plot(c(t0,age),c(length0_hat[2],length_hat[,2]),type="l",lwd=2,col="blue",
       #ylim=range(c(length0_hat[2],length_hat,obs[i,1,],obs0[i,1])),
       ylim=c(5,30),
       xlab="",
       ylab="",
       main="Length (cm)",
       yaxt="n"
  )
  axis(2,cex.axis=0.8)
  lines(c(t0,age),c(length0_hat[1],length_hat[,1]),lty=2,col="blue")
  lines(c(t0,age),c(length0_hat[3],length_hat[,3]),lty=2,col="blue")
  points(age,obs[i,1,],pch=19)
  points(t0,obs0[i,1],pch=19)
  text(500,7,
       paste(
         #i,
         #"\n",
         ID.label[i],
         "\n", 
         "pAm=",round(pAm.hat[i],0),"\n v=",round(v.hat[i],3),sep=""),cex=1)
  abline(h=median(obs[,1,n]),lty=3)
  
  #weight
  plot(c(t0,age),c(weight0_hat[2],weight_hat[,2]),type="l",lwd=2,col="blue",
       #ylim=range(c(weight0_hat[2],weight_hat,obs[i,2,],obs0[i,2])),
       ylim=c(20,710),#
       xlab="",
       ylab="",
       main="Wet weight (gr)",
       yaxt="n"
  )
  axis(2,cex.axis=0.8)
  lines(c(t0,age),c(weight0_hat[1],weight_hat[,1]),lty=2,col="blue")
  lines(c(t0,age),c(weight0_hat[3],weight_hat[,3]),lty=2,col="blue")
  points(age,obs[i,2,],pch=19)
  points(t0,obs0[i,2],pch=19)
  abline(h=median(obs[,2,n]),lty=3)
}

