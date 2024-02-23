#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 5
#---------------------------------------
# Last update: 30 January 2024


#-----
# 1: Loading libraries and data
#-----
remove(list=ls())
library(cmdstanr)  
library(bayesplot) 
library(posterior) 
library(ggplot2)   
library(deSolve) 
library(dfoptim)

load("input.Rdata") 
fit = readRDS("out.RDS") 
draws=fit$draws(format="matrix")

#-----
# 2: Updating DEB parameters with the results from Analysis.R
#-----

parms2 = parms 
parms2$Hp = 145400 # add-my-pet
parms2$birth_correction = birth_correction 
parms2$pAm = median(data.frame(as_draws_df(fit$draws("pAm_mu")))[,1])
parms2$v = median(data.frame(as_draws_df(fit$draws("v_mu")))[,1])
parms2$del_M = 0.2805034 # supp mat 6
parms2$W0 = 27652.82 # supp mat 3
parms2$W1 = median(data.frame(as_draws_df(fit$draws("W1")))[,1])
state = NULL
state["E"] = median(data.frame(as_draws_df(fit$draws("E0_mu")))[,1])
state["V"] = median(data.frame(as_draws_df(fit$draws("V0_mu")))[,1])
state["H"] = 18137.6         
state["R"] = 1e-30           

#-----
# 3: Preliminary solving DEB with the model in the main text
#-----

source("DEB_spawning0.R")  # without spawning dynamics
times=seq(t0,max(age),1)
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms2,
          method="ode45")
res=data.frame(res)

#-----
# 4: Estimating H0 and Hp for getting puberty at age_p
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

state["H"]=opt$par[1]
parms2$Hp=opt$par[2]
times=seq(t0,max(age),1)
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms2,
          method="ode45")

res=data.frame(res)
res$time[which(res$H>=parms2$Hp)][1]
res.no_spawning=res

#-----
# 5: Estimating max spawn fraction for getting gonad weight = 1% of total wet weight
#-----

temp=seq(t0,age[n],1)
t2=temp-(floor(temp/365.0)*365.0)
max.spawn=1
spawn=max.spawn*(exp(-(t2-330.0)^2.0/(2.0*30.0^2)) + exp(-(t2+365.0-330.0)^2/(2.0*30.0^2.0)))


source("DEB_spawning1.R")  
parms3=parms2
parms3$max.spawn=0.5
parms3$peak.spawn=330
parms3$spread.spawn=30
times=seq(t0,age[n],1)
res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms3,
          method="ode45")

res=data.frame(res)


times=seq(t0,max(age),1)

foo=function(x){
  parms3$max.spawn=x[1]
  parms3$peak.spawn=x[2]
  parms3$spread.spawn=x[3]
  res = ode(y=state,
            times=times,
            func=DEBmodel,
            parms=parms3,
            method="ode45")

  res=data.frame(res)
  temp=100*res$gonad_weight[length(times)] /
    res$wet_weight[length(times)] 
  (temp-1.0)^2 
}


opt=nmkb(p=c(0.5,330,30), f=foo, lower=c(0,200,20), upper=c(1,365,50)) 

# estimated peak date
as.Date("2000-01-01")+opt$par[2]+birth_correction

# estimated parameters for the spawning dynamics
parms3$max.spawn=opt$par[1]
parms3$peak.spawn=opt$par[2]
parms3$spread.spawn=opt$par[3]

res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms3,
          method="ode45")

res=data.frame(res)

# checking optimization success
100*res$gonad_weight[length(times)] / 
  res$wet_weight[length(times)]

res.spawning=res

# maximum deference in wet weight 
max(res.spawning$wet_weight-res.no_spawning$wet_weight)


# plot(res.no_spawning$time,res.no_spawning$wet_weight,type="l",col="blue",lwd=2,xlab="Fish age (days)",ylab="Wet weight (gr)")
# lines(res.no_spawning$time,res.spawning$wet_weight,col="red",lwd=2)
# legend("bottomright", c("DEB without reproduction","DEB With reproduction"),
#        col=c("blue","red"),lwd=2,cex=0.8)


#-----
# Figure S7: Expected weight-at-age predicted with (1) a model that includes 
# reproductive dynamics (red line) and (2) a model ignoring reproductive dynamics 
# (blue line).
#----- 

p=ggplot() +
  geom_line(data = res.no_spawning, aes(x = time, y = wet_weight, color = "DEB without reproduction"), size = 1, alpha = 1) +
  geom_line(data = res.spawning, aes(x = time, y = wet_weight, color = "DEB with reproduction"), size = 1, alpha = 1) +
  labs(x = "Fish age (days)", y = "Wet weight (g)") +
  scale_color_manual(
    labels = c("DEB without reproduction", "DEB with reproduction"),
    values = c("DEB without reproduction" = "#333399", "DEB with reproduction" = "red")
  ) +
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

p=p+guides(alpha = FALSE, shape = FALSE)
p=p + theme(legend.position = c(0.8, 0.15))
p=p + theme(legend.title = element_blank())
p=p+guides(color = guide_legend(override.aes = list(fill = NA))) 
p=p+theme(legend.background = element_rect(fill='transparent'),
          legend.box.background = element_blank())
p


# Saving:

# ggsave(filename=paste("S7.png"),path=paste0(getwd(),"/Out",sep=""),width = 8, height = 5)







