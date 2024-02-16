#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 1
#---------------------------------------
# Last update: 30 January 2024

# Fitting daily temperature to a sinusoidal function using STAN

#-----
# 1: Loading data and libraries
#-----
remove(list=ls())
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(bayesplot)
library(posterior)

# Daily data: Temperature and other environmental variables

daily = read.csv("daily data.csv",dec=",")
daily$fecha=as.Date(as.character(daily$fecha),"%d/%m/%Y")

data.sin=list(pi=pi)
temp=which(is.na(daily$temperatura))
data.sin$time=as.numeric(daily$fecha[-temp]-as.Date("2017-01-01"))
data.sin$temperature=daily$temperatura[-temp]
data.sin$n=length(data.sin$time)

#-----       
# 2: STAN model
#-----
sink("DEB.stan")
cat(" // first line
functions {
}

data {
  int n;
  array [n] real temperature;
  array [n] real time;
  
  array [2] real prior_temp_mean;
  array [2] real prior_temp_amp;
  array [2] real prior_temp_phi_days;
  array [2] real prior_temp_sd;
}

transformed data {
}

parameters {  // DEB parameters
   real <lower=10.0, upper=30.0> temp_mean;
   real <lower=0.001, upper=5.0> temp_amp;
   real <lower=0.0, upper=365.0> temp_phi_days;
   real <lower=0.0> temp_sd;
}

transformed parameters {
  real temp_phi;
  temp_phi = 2.0*pi()*temp_phi_days/365.0;
}   

model {
   array [n] real temperature_hat;
   
   temp_mean ~ normal (prior_temp_mean[1],prior_temp_mean[2]);
   temp_amp ~ normal (prior_temp_amp[1],prior_temp_amp[2]);
   temp_phi_days ~ cauchy (prior_temp_phi_days[1],prior_temp_phi_days[2]);
   temp_sd ~ cauchy (prior_temp_sd[1],prior_temp_sd[2]);
   
   for (i in 1:n){
     temperature_hat[i] = temp_mean + temp_amp*sin(2.0*pi()*time[i]/365.0 + temp_phi);
     temperature[i] ~ normal (temperature_hat[i],temp_sd);
   }
}

generated quantities {
}

" 
,fill = TRUE)
sink()

#-----
# 3: Compiling model
#-----
mod = cmdstan_model("DEB.stan")

#-----
# 4: Initializating the model
#-----
 
prior_temp_mean=c(mean(data.sin$temperature),100) 
prior_temp_amp=c(diff(range(data.sin$temperature))/2,10)
prior_temp_phi_days=c(365/2,300) 
prior_temp_sd=c(1,10)

n.chains = 2
initializer = function() list(
  "temp_mean"= mean(data.sin$temperature),
  "temp_amp"= diff(range(data.sin$temperature))/2,
  "temp_phi_days"= 365/2,
  "temp_sd"=1
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()

#-----
# 5: Running
#-----

fit = mod$sample(
  data =list (
    n=data.sin$n,
    temperature=data.sin$temperature,
    time=data.sin$time,
    
    prior_temp_mean=prior_temp_mean,
    prior_temp_amp=prior_temp_amp,
    prior_temp_phi_days=prior_temp_phi_days,
    prior_temp_sd=prior_temp_sd
  ),
  chains = n.chains,
  parallel_chains = n.chains,
  iter_warmup = 1000,
  iter_sampling = 1000,
  init = inits,
  max_treedepth = 10,
  adapt_delta = 0.8
)

#-----
# 6: Saving
#-----
# fit$save_object(file = "out_temp.RDS")
# sink(file = "out_dignosi_temp.txt")
# fit$cmdstan_diagnose()
# sink()
# file.rename("DEB.stan","out_temp_model.R")

draws=fit$draws(format="matrix")

#-----
# Figure S1: Observed daily temperature (grey line) and estimated daily 
# temperature from the sinusoidal function adjusted to the data (blue line). 
# This function was used to facilitate numerical integration when estimating 
# DEB parameters.
#-----

#Preparing data for plotting:

time_pre=seq(min(data.sin$time),max(data.sin$time),1)
temp_mean_hat=median(data.frame(as_draws_df(fit$draws("temp_mean")))[,1])
temp_amp_hat=median(data.frame(as_draws_df(fit$draws("temp_amp")))[,1])
temp_phi_hat=median(data.frame(as_draws_df(fit$draws("temp_phi")))[,1])
Temp.hat=temp_mean_hat + temp_amp_hat*sin(2*pi*time_pre/365 + temp_phi_hat)

# Create a data frame for the predicted temperature values:

predicted_data = data.frame(
  fecha = as.Date("2017-01-01") + time_pre,
  Temperature = Temp.hat
)

# Producing Figure S1

ggplot() +
  geom_line(data = daily, aes(x = fecha, y = temperatura), color = "darkgrey") +
  geom_line(data = predicted_data, aes(x = fecha, y = Temperature), color = "#333399", size =1 ) +
  labs(x = "Date", y = "Temperature (Celsius degrees)") +
  theme(
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))

# Saving:

#ggsave(filename=paste("S1.png"),path=paste0(getwd(),"/Out",sep=""), width = 8, height = 4)


