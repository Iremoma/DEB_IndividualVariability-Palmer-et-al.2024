# DEB_IIV_STAN  Palmer-et-al.-2024

This repository contains all the code related to the manuscript entitled **Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications**. 
Authors: Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga, Marine Herlin, Orestis Stavrakidis-Zachoua and Andrea Campos-Candela.

This repository is licensed under a xxxxx. 

# Contents
In this repository we provide all the scripts necessary to reproduce the data exploration and analysis, simulations, and results presented on the paper Palmer _et al._, 2024. **Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications**.

We provide the STAN routine that we have developed for estimating four DEB parameters at the individual level for 69 cultured Gilt-head breams (_Sparus aurata_) for which up to 11 measures of length and wet weight were available. 

For this, all the code is organized in the following folders: 
1. [Analyses from real data](#Analyses_from_real_data)
2. [Simulations](./2_SIMULATIONS)
3. [Supplementary materials](./3_SUPPLEMENTARY_MATERIALS)
 
# Requirements
Install: 
- R (4.1.3)
- CmdStanR (Command Stan R) in R: CmdStanR is a lightweight interface to Stan for R users. If you are using Stan in R for the first time, or you have no experience with this package, we encourage you to visit [Getting started with CmdStanR]( https://mc-stan.org/cmdstanr/articles/cmdstanr.html)  and follow the recommended steps to install the package.

Recommended properties of the computer, and those that characterized the computer used for running these models, are: 
- Intel(R) Xeon(R)  CPU E5-2620 v4 @ 2.10GHz (40 cores and 64 GB RAM)

# Working flow

## 1. [Analyses from real data](./1_REALDATA/README.md)
Here we provide the main code for the paper Palmer _et al._ 2024. By running [Analysis.R](./REALDATA/Analysis.R), you will estimate four DEB parameters from the provided observed data in [_input.Rdata_](./REALDATA/input.Rdata) and obtain the output (which is also provided in the file called _out.RDS_) from which all the main content of the paper and figures have been obtained [_ExploringResult.R_](./REALDATA/ExploringResult.R). Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps. 




# Citation
If you use this code, please consider citing our work:
Palmer, M., Moro-Martínez, I., Tomás-Ferrer, J., Orestis, Marilo, Marine, Amalia, Lika, Campos-Candela, A. (2024) Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications.


Any questions: .....@......
