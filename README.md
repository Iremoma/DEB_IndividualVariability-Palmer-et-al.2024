# DEB_IIV_STAN  Palmer-et-al.-2024

This repository contains all the code related to the manuscript entitled **Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications**, by Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga, Marine Herlin, Orestis Stavrakidis-Zachoua and Andrea Campos-Candela.

This repository is licensed under a xxxxx. 

# Requirements
Install: 
- [R](https://www.r-project.org/about.html) (4.1.3 or higher)
- CmdStanR (Command Stan R) in R: CmdStanR is a lightweight interface to Stan for R users. If you are using Stan in R for the first time, or you have no experience with this package, we encourage you to visit [Getting started with CmdStanR]( https://mc-stan.org/cmdstanr/articles/cmdstanr.html)  and follow the recommended steps to install the package.

Recommended properties of the computer, and those that characterized the computer used for running these models, are: 
- Intel(R) Xeon(R)  CPU E5-2620 v4 @ 2.10GHz (40 cores and 64 GB RAM)

# Contents and working flow
We provide all the R scripts necessary to reproduce the data exploration and analysis, simulations, and results presented on the paper Palmer _et al._, 2024. **Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps.** All the code is organized in the following folders: 

## 1. [Analyses from real data](./1_REALDATA)
In this section you will find the STAN routine developed for estimating four DEB parameters at the individual level for 69 cultured Gilt-head breams (_Sparus aurata_) for which up to 11 measures of length and wet weight were available. We provide R code for exploring the empirical data, implementeing the STAN model, and exploring the results of the model. 

## 2. [Simulations](./2_SIMULATIONS)

## 3. [Supplementary materials](./3_SUPPLEMENTARY_MATERIALS) 

# Citation
If you use this code, please consider citing our work:

Palmer, M.; Moro-Martinez, I.; Tomas-Ferrer, J.; Grau,A.; Lopez-Belluga, M.D.; Herlin, M.; Stavrakidis-Zachoua, O. and Campos-Candela, A. (2024) Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications.


Any questions: palmer@imedea.uib-csic.es
