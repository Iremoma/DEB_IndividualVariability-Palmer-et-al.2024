# DEB_IndividualVariability Palmer _et-al._ 2024

This repository contains all the code related to the manuscript entitled **Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications**, by Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga, Marine Herlin, Orestis Stavrakidis-Zachoua and Andrea Campos-Candela.

The work presented here is a contribution towards gaining greater representation of the between-individual variability in Dynamic Energy Budget (**DEB**) models. [More about DEB theory [here](https://debportal.debtheory.org/docs/)]. We present our analytical strategy based on Bayesian inference and hierarchical models for inferring variability in DEB parameters from data across individuals (also known as between-individual variability). The analytical strategy is developed in **STAN**, which is establishing itself as a general-purpose statistical software suitable for the process of fitting complex mechanistic and dynamic models. [More about STAN [here](https://mc-stan.org/)].

This repository is licensed under a [MIT License Copyright (c) 2024 Miquel Palmer (et.al. 2024)](./LICENSE). 

# Requirements
Install: 
- [R](https://www.r-project.org/about.html) (4.1.3 or higher)
- CmdStanR (Command Stan R) in R: CmdStanR is a lightweight interface to Stan for R users. If you are using Stan in R for the first time, or you have no experience with this package, we encourage you to visit [Getting started with CmdStanR]( https://mc-stan.org/cmdstanr/articles/cmdstanr.html)  and follow the recommended steps to install the package.

Recommended properties of the computer, and those that characterized the computer used for running these models, are: 
- Intel(R) Xeon(R)  CPU E5-2620 v4 @ 2.10GHz (40 cores and 64 GB RAM)

# Contents and working flow
We provide all the R scripts necessary to reproduce the data exploration and analysis, simulations, and results presented on the paper Palmer _et al._, 2024. **Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps.** By running the code you will be guided in the process of thinking and exploration we performed in this work, and you will be able to reproduce the estimation of parameters and explore the results by reproducing (in a basic/raw format) all the figures in the main text and supplementeray materials. 

All the code is organized in the following folders: 

## 1. [Analyses from real data](./1_REALDATA)
In this section you will find the STAN routine developed for estimating four DEB parameters at the individual level for 69 cultured Gilt-head breams (_Sparus aurata_) for which up to 11 measures of length and wet weight were available. We provide R code for exploring the empirical data, implementeing the STAN model, and exploring the results of the model. 

## 2. [Simulations](./2_SIMULATIONS)
Data-simulation experiments were completed for demonstrating opportunities and limitations of our analytical strategy followed in the section above for estimating DEB parameters at the individual level. 

## 3. [Supplementary materials](./3_SUPPLEMENTARY_MATERIALS) 
The codes of the supplementary materials of article Palmer _et al._ 2024 show a battery of analyzes and/or simulations carried out to explore complementary parts to the general analytical strategy presented in the sections [1. Analyses from real data](./1_REALDATA) and [2. Simulations](./2_SIMULATIONS). We advise users to browse these folders with the texts of the article's supplementary materials in order to have all the information and complete justification for each analysis.

# Citation
If you use this code, please consider citing our work:

Palmer, M.; Moro-Martinez, I.; Tomas-Ferrer, J.; Grau,A.; Lopez-Belluga, M.D.; Herlin, M.; Stavrakidis-Zachoua, O. and Campos-Candela, A. (2024) Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications. _Submitted in a journal and awaiting for evaluation and review from the journal_.


Any questions: palmer@imedea.uib-csic.es
