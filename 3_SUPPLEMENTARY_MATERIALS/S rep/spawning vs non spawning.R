#---------------------------------------
# Estimating DEB parameters for gilthead seabream (Sparus aurata)
#---------------------------------------
# Last update: 11-08-2023
# Supplementary maeria S5: comparing 
#
#-----
# 1: Loading libraries and data
#-----
remove(list=ls())
library(cmdstanr)  # stan
library(bayesplot) # plotting library
library(posterior) # plotting library
library(ggplot2)   # plotting library
library(deSolve)   # numerical integration

load("input.Rdata")   # inout
fit = readRDS("out_v3.RDS") # STAN resukts (model in the main text)
draws=fit$draws(format="matrix")

#-----
# 2: Updating DEB parameters with the results from analysis_v2.R
#-----
parms2 = parms # cloning DEB parameters
parms2$Hp = 145400 # add-my-pet
parms2$birth_correction = birth_correction # days from 1st january
parms2$pAm = median(data.frame(as_draws_df(fit$draws("pAm_mu")))[,1])
parms2$v = median(data.frame(as_draws_df(fit$draws("v_mu")))[,1])
parms2$del_M = 0.2805034
parms2$W0 = 27652.82
parms2$W1 = median(data.frame(as_draws_df(fit$draws("W1")))[,1])#1118730
state = NULL
state["E"] = median(data.frame(as_draws_df(fit$draws("E0_mu")))[,1])
state["V"] = median(data.frame(as_draws_df(fit$draws("V0_mu")))[,1])
state["H"] = 18137.6         # set for getting Age at maturity = 365
state["R"] = 1e-30           # Reproduction energy

#-----
# 3: preliminary solving DEB with the model in the main text
#-----
source("DEB_spawning0.R")  # deb model without spawning
times=seq(t0,max(age),1)
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms2,
          method="ode45")
plot(res)
res=data.frame(res)

#-----
# 4: estimating H0 and Hp for getting puberty at age_p
#-----
times=seq(t0,max(age),.1)
parms2$age_p = 300
foo=function(x){
  state["H"]=x[1]
  parms2$Hp=x[2]
  res = ode(y=state,
            times=times,
            func=DEBmodel,
            parms=parms2,
            method="ode45")
  res=data.frame(res)
  temp=res$time[which(res$H>=parms2$Hp)][1]
  (temp-parms2$age_p)^2
}

opt=optim(c(10000,50000),foo)
#[1] 27382.78 35889.07

# testing
state["H"]=opt$par[1]
parms2$Hp=opt$par[2]
times=seq(t0,max(age),1)
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms2,
          method="ode45")
plot(res)
res=data.frame(res)
res$time[which(res$H>=parms2$Hp)][1]
#[1] 300
# saving results
res.no_spawning=res

#-----
# 5: estimating max spawn fraction for getting 
#  gonad weight = 1.6% of total wet weight
#-----
#checking equation for setting the spawning temporal window
temp=seq(t0,age[n],1)
t2=temp-(floor(temp/365.0)*365.0)
max.spawn=1
spawn=max.spawn*(exp(-(t2-330.0)^2.0/(2.0*30.0^2)) + exp(-(t2+365.0-330.0)^2/(2.0*30.0^2.0)))
plot(temp,spawn,type="l")
abline(v=c(0,1,2,3)*365+330)

source("DEB_spawning1.R")  # deb model with spawning
parms3=parms2
parms3$max.spawn=0.5
parms3$peak.spawn=330
parms3$spread.spawn=30
# checking DEB with spawning
times=seq(t0,age[n],1) # caution! crowding function does not allow extrapolation
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms3,
          method="ode45")
#plot(res)
res=data.frame(res)
plot(res$time,res$gam_rel_weight,type="l")
par(new = TRUE) 
plot(temp,spawn,type="l",axes = FALSE,xlab="",ylab="",col="red",lty=3)
axis(side = 4, at = pretty(range(spawn))) 

plot(res$time,res$wet_weight,type="l")
abline(v=c(0,1,2,3)*365+parms3$peak.spawn)


times=seq(t0,max(age),1)
#x=c(parms3$max.spawn,parms3$spread.spawn)
foo=function(x){
  parms3$max.spawn=x[1]
  parms3$peak.spawn=x[2]
  parms3$spread.spawn=x[3]
  res = ode(y=state,
            times=times,
            func=DEBmodel,
            parms=parms3,
            method="ode45")
  #plot(res)
  res=data.frame(res)
  temp=100*res$gonad_weight[length(times)] / # gonad weight at tf
  res$wet_weight[length(times)]  # wet weight at tf
  #(temp-1.6)^2
  (temp-1.0)^2 # assuming that some eggs are hydrated
}
#opt=optim(par=c(0.5,330, 30),foo,lower=c(0,300,20),upper=c(1,360,100),method="L-BFGS-B")
# it does not converge

library(dfoptim)
opt=nmkb(p=c(0.5,330,30), f=foo, lower=c(0,200,20), upper=c(1,365,50)) # Converges
#[1]   0.03034681 338.17816951  36.55759845
# $message
# [1] "Successful convergence"

# estimated peak location in julian days
as.Date("2000-01-01")+opt$par[2]+birth_correction
#"2000-12-24"

#cheking
parms3$max.spawn=opt$par[1]
parms3$peak.spawn=opt$par[2]
parms3$spread.spawn=opt$par[3]
#checking the spawning temporal window
temp=seq(t0,age[n],1) 
t2=temp-(floor(temp/365.0)*365.0)
max.spawn=1
spawn=parms3$max.spawn*(exp(-(t2-parms3$peak.spawn)^2.0/(2.0*parms3$spread.spawn^2)) + exp(-(t2+365.0-parms3$peak.spawn)^2/(2.0*parms3$spread.spawn^2.0)))
plot(temp,spawn,type="l")
abline(v=parms3$peak.spawn+365*c(1,2,3))

res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms3,
          method="ode45")
#plot(res)
res=data.frame(res)
temp=100*res$gonad_weight[length(times)] / # gonad weight at tf
  res$wet_weight[length(times)]
temp
#[1] 1.000267 (%) (CHANGED TO 1%)
plot(res$time,res$gonad_weight,type="l")

plot(res$time,res$gam_rel_weight,type="l",ylab="weigth of gametes released per day (gr)")
abline(v=c(0,1,2,3)*365+parms3$peak.spawn,)
temp=seq(t0,age[n],1)
spawn=max.spawn*(exp(-(t2-parms3$peak.spawn)^2.0/(2.0*parms3$spread.spawn^2)) + exp(-(t2+365.0-parms3$peak.spawn)^2/(2.0*parms3$spread.spawn^2.0)))
par(new = TRUE) 
plot(temp,spawn,type="l",axes = FALSE,xlab="",ylab="",col="red",lty=3)
axis(side = 4, at = pretty(range(spawn))) 

plot(res$time,res$wet_weight,type="l")
abline(v=c(0,1,2,3,4)*365+parms3$peak.spawn)

res.spawning=res

#-------------
# Figure S7
#-------------
plot(res.no_spawning$time,res.no_spawning$wet_weight,type="l",col="blue",lwd=2,xlab="Fish age (days)",ylab="Wet weight (gr)")
lines(res.no_spawning$time,res.spawning$wet_weight,col="red",lwd=2)
legend("bottomright", c("DEB without reproduction","DEB With reproduction"),
       col=c("blue","red"),lwd=2)
#abline(v=c(t0,age),lty=3)

#abline(v=2*365+(parms3$peak.spawn),lty=2)
#abline(v=parms3$peak.spawn*c(1,2)-birth_correction,lty=2)

# points(age,apply(obs[,2,],2,median),pch=19,cex=0.7)
# points(t0,median(obs0[,2]),pch=19,cex=0.7)

max(res.spawning$wet_weight-res.no_spawning$wet_weight)
#[1] 19.68681
