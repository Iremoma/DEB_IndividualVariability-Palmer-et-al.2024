#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 10
#---------------------------------------
# Last update: 31 January 2024
#
#
# 1: Loading data and libraries
#-----
remove(list=ls())
library(cmdstanr)
library(bayesplot) 
library(posterior)
library(ggplot2) 
library(hexbin)
load("input.RData") # observed data
fit = readRDS("out.RDS") # results of the analysis in the main text

#-----
# Figure S11: Trace plots for the population-level mean and the between-fish 
# standard deviations for f*pAm  and v (represented as pAm and v respectively in
# the STAN code) suggesting proper mixing and no abnormal patterns.
#-----

mcmc_trace(fit$draws(c("pAm_mu","pAm_sd","v_mu","v_sd")))

#-----
# Figure S12: Trace plots for the measurement errors (from left to right, 
# ichthyometer, digital caliper and  digital balance) suggesting proper mixing 
# and no abnormal patterns.
#-----

mcmc_trace(fit$draws(c("sd_length","sd_weight")))

#-----
# Figure S13: Pairs plot of f*pAm and v population's mean (represented as pAm_mu
# and v_mu respectively). Note that posteriors are uncorrelated, as denoted by 
# the absence of correlation or any abnormal pattern.
#-----
mcmc_pairs(fit$draws(c("pAm_mu","v_mu")), off_diag_fun = c("hex"))

#-----
# Figure S14: Pairs plots of f*pAmi and vi of a given fish (i=10). Note that posteriors
# are uncorrelated, as denoted by the absence of correlation or any abnormal pattern.
#-----

mcmc_pairs(fit$draws(c("pAm[10]","v[10]")), off_diag_fun = c("hex"))


#-----
# Figure S15: Priors (blue line) versus posteriors (red line). Note that 
# priors' distributions are always wider than posteriors' distributions, 
# suggesting that prior selection is not constraining parameters' estimates 
# (i.e., priors are virtually non-informative).
#-----

# Priors (copied from Analysis.R):

prior_E0_mu=c(50000,500000)
prior_E0_sd=c(20000,200000)
prior_V0_mu=c(30,100)
prior_V0_sd=c(10,100)
prior_pAm = c(parms$pAm_after , parms$pAm_after) 
prior_pAm_sd = c(50,100) 
prior_v = c(0, 2) 
prior_v_sd = c(0, 1) 
prior_sd_length=c(0,10) 
prior_sd_weight=c(0,100) 

t=seq(t0,age[n])
biomass = 3.173425e+03-3.720773e+02*t+2.362952e+00*t^2-1.431395e-03*t^3; 
prior_W1=c( mean(biomass), 10*mean(biomass)) 


iter=length(data.frame(as_draws_df(fit$draws("pAm")))[,1])

# Producing Figure S15:

par(mfrow=c(4,3))

#1: E0_mu

temp=rnorm(iter,prior_E0_mu[1],prior_E0_mu[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("E0_mu")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(E[0~mu]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#2: E0_sd

temp=rnorm(iter,prior_E0_sd[1],prior_E0_sd[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("E0_sd")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(E[0~sd]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#3: V0_mu

temp=rnorm(iter,prior_V0_mu[1],prior_V0_mu[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("V0_mu")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(V[0~mu]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#4: V0_sd

temp=rnorm(iter,prior_V0_sd[1],prior_V0_sd[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("V0_sd")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(V[0~sd]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#5: pAm_mu

temp=rnorm(iter,prior_pAm[1],prior_pAm[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("pAm_mu")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(bold(symbol("\042"))[~mu]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#6: pAm_sd

temp=rcauchy(iter,prior_pAm_sd[1],prior_pAm_sd[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("V0_sd")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(bold(symbol("\042"))[~sd]),xlab="",
     #xlim=c(0,max(density_prior$x)),
     xlim=c(0,1000),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#7: v_mu

temp=rcauchy(iter,prior_v[1],prior_v[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("v_mu")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(dot(nu)[~mu]),xlab="",
     #xlim=c(0,max(density_prior$x)),
     xlim=c(0,5),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#8: v_sd

temp=rcauchy(iter,prior_v_sd[1],prior_v_sd[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("v_sd")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(dot(nu)[~sd]),xlab="",
     #xlim=c(0,max(density_prior$x)),
     xlim=c(0,5),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#9: sd_length

temp=rcauchy(iter,prior_sd_length[1],prior_sd_length[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("sd_length")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(Length[~sd]),xlab="",
     #xlim=c(0,max(density_prior$x)),
     xlim=c(0,10),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
temp=data.frame(as_draws_df(fit$draws("sd_length")))[,2]
density_post=density(temp)
lines(density_post$x,density_post$y,col="red")

#10: v_sd

temp=rcauchy(iter,prior_sd_weight[1],prior_sd_weight[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("sd_weight")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(Weight[~sd]),xlab="",
     #xlim=c(0,max(density_prior$x)),
     xlim=c(0,50),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

#11: W1

temp=rnorm(iter,prior_W1[1],prior_W1[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("W1")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(W[~1]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")

