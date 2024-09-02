#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 3
#---------------------------------------
# Last update: 30 January 2024


#-----
# 1: Loading libraries and data
#-----
remove(list=ls())
library(gridExtra)
library(ggplot2)
library(cowplot)
library(MASS)

remove(list=ls())
load("biomass.Rdata")
# age: fish age (days) at each sampling event
# biomass: fish biomass (kg) at each sampling event

# Selecting the optimal degree for a polynomial fit:

fit_1 = lm(biomass ~ age)
fit_2 = lm(biomass~poly(age,2))
fit_3 = lm(biomass~poly(age,3))
fit_4 = lm(biomass~poly(age,4))
fit_5 = lm(biomass~poly(age,5))
# print(anova(fit_1,fit_2,fit_3,fit_4,fit_5))


# fitting
lm = lm(biomass ~ I(age)+I(age^2)+I(age^3))
empty.model = lm(biomass ~ 1)
temp=stepAIC(empty.model, trace=TRUE, direction="forward",scope=formula(lm))

pre=as.numeric(predict(lm))
temp=seq(min(age),max(age),1)
pre=3.173425e+03-3.720773e+02*temp+2.362952e+00*temp^2-1.431395e-03*temp^3 

#-----
# Figure S3: Observed biomass (blue dots) and a polynomial function adjusted to
# the data (grey line). This function was used to facilitate numerical integration 
# when estimating DEB parameters in STAN.
#----- 

# Producing Figure S3:

ggplot() +
  geom_point(aes(x = age, y = biomass/1000/20), shape = 19, size = 2,color="#333399") +
  geom_line(aes(x=temp,y=pre/1000/20),color="darkgrey",size=1)+   
  ylab(bquote("Biomass (Kg? "*m^-3*")"))+
  xlab("Fish age (days)") +
  theme(
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"))


# Saving: 

# ggsave(filename=paste("S3.png"),path=paste0(getwd(),"/Out",sep=""), width = 6, height = 3)
