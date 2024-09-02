#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 2
#---------------------------------------
# Last update: 30 January 2024

#----
# 1: Loading libraries and data
#----
remove(list=ls())
library(cmdstanr) 
library(deSolve)
library(posterior)
library(gridExtra)
library(ggplot2)
library(cowplot)
load("input.RData") 
fit = readRDS("out.RDS") 

# cloning parms (generic DEB parameters) to parms2 and setting the estimated
#  values for the actual fish analyzed
parms2=parms
parms2$birth_correction=as.numeric(as.Date("2017-07-26")-as.Date("2017-01-01"))
parms2$pAm=median(data.frame(as_draws_df(fit$draws("pAm_mu")))[,1])
parms2$v=median(data.frame(as_draws_df(fit$draws("v_mu")))[,1])
parms2$del_M=0.2805034
state=NULL
state["E"] = median(data.frame(as_draws_df(fit$draws("E0_mu")))[,1])
state["V"] = median(data.frame(as_draws_df(fit$draws("V0_mu")))[,1])

source("DEB.R")  
# Caution! this is a version (1) without reproduction, and (2) without crowding
#  because the aim is to compare growth using the actual daily temperature
#  versus the temperature as estimated using the sinusoidal function in S1

#----
# 2: Solving using the sinusoidal function
#----
times=seq(t0,age[n],1)
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms2,
          method="ode45")
#plot(res)
res=data.frame(res)

#----
# 3: Solving using daily temperature
#----
load("Temperature.RData")
# data.sin$time # days from 2017-01-01
# data.sin$temperature # daily temperature (Celsius degrees)
# data.sin$n # number of obervations
# data.sin$pi # pi
times.forcing=data.sin$time-data.sin$time[1]
forcing.var = approxfun(times.forcing, data.sin$temperature, rule = 2)
parms_daily=parms2
parms_daily$forcing.var=forcing.var

source("DEB_daily.R")
# This is the same DEB version as above but it uses observed daily temperature as input
res2 = ode(y=state,
           times=times,
           func=DEBmodel,
           parms=parms_daily,
           method="ode45")
#plot(res2)
res2=data.frame(res2)

#-----
# Figure S2: Left panel: Observed daily temperature (blue dots) and estimated 
# daily temperature from the R interpolation function approxfun (blue line). 
# Right panel: Expected length-at-age for a single fish when temperature is 
# inputted with the sinusoidal function introduced in S1 (red line) or with the 
# interpolation function from the left panel (blue line). In both cases, the 
# vertical dashed lines denote the period during which fish were monitored for 
# length and weight.
#----- 

# Preparing data for plotting:

df_res = data.frame(
  time = res$time,
  total_length = res$total_length,
  Model = "Sinusoidal Function"
)

df_res2 = data.frame(
  time = res2$time,
  total_length = res2$total_length,
  Model = "Daily Temperature"
)

df = rbind(df_res, df_res2)

# Producing Figure S2:
data.sin$fecha=data.sin$time+as.Date("2017-01-01")
data.sin=data.frame(data.sin)

sf=ggplot(data.sin, aes(x = fecha, y = temperature)) +
  geom_point(shape = 19, size = 1,color="#333399") +
  geom_line(aes(group = 1), color = "darkgrey") +
  labs(x = "Date", y = "Temperature (Celsius degrees)") +
  theme(
    plot.margin = unit(c(4, 4, 0, 1), "lines"),
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))+
  
  geom_vline(xintercept = as.numeric(236)+as.Date("2017-01-21"), linetype = "dashed", color = "#9999FF") +
  geom_vline(xintercept = as.numeric(738)+as.Date("2017-01-21"), linetype = "dashed", color = "#9999FF") 



fit=ggplot(df, aes(x = time, y = total_length, color = Model)) +
  geom_line(size = 2) +
  labs(x = "Fish age (days)", y = "Length (cm)" ) +
  scale_color_manual(
    values = c("Sinusoidal Function" = "#333399", "Daily Temperature" = "red"),
    labels = c("Sinusoidal Function", "Daily Temperature")
  ) +
  geom_vline(xintercept = as.numeric(236), linetype = "dashed", color = "#9999FF") +
  geom_vline(xintercept = as.numeric(738), linetype = "dashed", color = "#9999FF") +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.position = "none"  
  ) +
  labs(
    color = "Modelo"
  ) +
  annotate("text", x = 600, y = 12, label = "Sinusoidal Function", color = "#333399", size = 3) +
  annotate("text", x = 600, y = 14, label = "Daily Temperature", color = "red", size = 3)



ggdraw (
  plot_grid (sf,
             fit,
             ncol=2,align='h',rel_widths = c(2, 1)))

# Saving:

# ggsave(filename=paste("S2.png"),path=paste0(getwd(),"/Out",sep=""),
#        width=25,
#        height=10,
#        units="cm")
