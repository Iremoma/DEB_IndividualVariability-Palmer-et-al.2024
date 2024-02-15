#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 7
#---------------------------------------
# Last update: 25 January 2024

#-----
# 1: Loading data, results and libraries
#-----
remove(list=ls())
load("input_sim.RData")  

library(cmdstanr)  
library(bayesplot) 
library(posterior) 
library(ggplot2) 
library(cowplot)

fit = readRDS("out_sim.RDS")  
draws=fit$draws(format="matrix")

#-----
# 1: Producing Figure 5: Results of the simulation experiment. Estimated versus
# true values for f*pAm and v. Each point corresponds to the values of an individual 
# fish. The CI interval denotes 95% of the posterior distribution. Simulation 
# settings emulate the actual scenario for S. aurata in terms of between-individual 
# variability, number of fish, and number of repeated measures per fish. The line 
# represents perfect prediction.
#-----

# Nedeed data from out_sim.RDS to create the graphic fpAm:

pAm.hat = data.frame(as_draws_df(fit$draws("pAm")))[,1:N]

data.fpAm = data.frame(
  True = sim.par[,1],
  Estimated = apply(pAm.hat, 2, median),
  Lower = apply(pAm.hat, 2, quantile, 0.225),
  Upper = apply(pAm.hat, 2, quantile, 0.975)
)

# Graphic for f*pAm:

fpAm=ggplot(data.fpAm, aes(x = True, y = Estimated)) +
  geom_point(shape = 19, size = 0.75,color="#333399") +
  geom_abline(intercept = 0, slope = 1,color="darkgrey") +
  geom_segment(aes(xend = True, yend = Lower, x = True, y = Upper), linetype = "solid",color="#333399") +
  labs(x = "True", y = "Estimated",title=bquote(bold(symbol("\042")))) +
  ylim(range(data.fpAm$Lower, data.fpAm$Upper)) +
  theme(
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.title = element_text(hjust = 0.5)
  )


# Nedeed data from out_sim.RDS to create the graphic v:

v.hat = data.frame(as_draws_df(fit$draws("v")))[,1:N]

data.v = data.frame(
  True = sim.par[,2],
  Estimated = apply(v.hat, 2, median),
  Lower = apply(v.hat, 2, quantile, 0.225),
  Upper = apply(v.hat, 2, quantile, 0.975)
)

# Graphic for v:

v=ggplot(data.v, aes(x = True, y = Estimated)) +
  geom_point(shape = 19, size = 0.75,color="#333399") +
  geom_abline(intercept = 0, slope = 1,color="darkgrey") +
  geom_segment(aes(xend = True, yend = Lower, x = True, y = Upper), linetype = "solid",color="#333399") +
  labs(x = "True", y = "Estimated",title=bquote(dot(nu))) +
  ylim(range(data.v$Lower, data.v$Upper)) +
  theme(
    axis.title.x=element_text(size=10),
    axis.title.y=element_text(size=10),
    axis.text = element_text(size = 8),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="transparent"),
    plot.title = element_text(hjust = 0.5)
  )


ggdraw (
  plot_grid (fpAm,
             v,
             ncol=2,align='h',rel_widths = c(1, 1)))


# ggsave(filename=paste("Fig5.png"),path=paste0(getwd(),"/Out",sep=""),
#        width=25,
#        height=10,
#        units="cm")
