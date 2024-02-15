# DEB_IIV_STAN  Palmer-et-al.-2024

This repository contains all the code related to the manuscript entitled **Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications**. 
Authors: Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga, Marine Herlin, Orestis Stavrakidis-Zachoua and Andrea Campos-Candela.

This repository is licensed under a xxxxx. 

# Contents
In this repository we provide all the scripts necessary to reproduce the data exploration and analysis, simulations, and results presented on the paper Palmer _et al._, 2024. **Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications**.

We provide the STAN routine that we have developed for estimating four DEB parameters at the individual level for 69 cultured Gilt-head breams (_Sparus aurata_) for which up to 11 measures of length and wet weight were available. 

For this, all the code is organized in the following folders: 
1. [Analyses from real data](./1_REALDATA)
2. [Simulations](./2_SIMULATIONS)
3. [Supplementary materials](./3_SUPPLEMENTARY_MATERIALS)
 
# Requirements
Install: 
- R (4.1.3)
- CmdStanR (Command Stan R) in R: CmdStanR is a lightweight interface to Stan for R users. If you are using Stan in R for the first time, or you have no experience with this package, we encourage you to visit [Getting started with CmdStanR]( https://mc-stan.org/cmdstanr/articles/cmdstanr.html)  and follow the recommended steps to install the package.

Recommended properties of the computer, and those that characterized the computer used for running these models, are: 
- Intel(R) Xeon(R)  CPU E5-2620 v4 @ 2.10GHz (40 cores and 64 GB RAM)

# Working flow

## Analyses from real data
Here we provide the main code for the paper Palmer _et al._ 2024. By running [Analysis.R](./REALDATA/Analysis.R), you will estimate four DEB parameters from the provided observed data in [_input.Rdata_](./REALDATA/input.Rdata) and obtain the output (which is also provided in the file called _out.RDS_) from which all the main content of the paper and figures have been obtained [_ExploringResult.R_](./REALDATA/ExploringResult.R). Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps. 

Here you find the following files: 

### Core script: 
1. [Analysis.R](./REALDATA/Analysis.R) This script containing the code for estimating four DEB parameters from the observed data contained in [_input.Rdata_](./REALDATA/input.Rdata). This script has served as the skeleton of the work presented in the main text in Palmer _et al._ 2024. It can be run in multicore or in local version. As the script is presented here, it is prepared for be run in a multicore version (the  local version is on the script but currently locked). Note however that if you do not have multicore, running this script may take a considerable amount of time (in the order of several days). For this reason, we provide the output of the model as well.
   - *Input*: _input.RData_
   - *Outputs*: _out.RDS_; _out_dignosi.txt_;_out_model.R_

 #### Input data: 
- [_input.Rdata_](./REALDATA/input.Rdata) This file contains the experimental data used for estimating the parameters of the DEB model for our case study (more datils in Palmer _et al._ 2024). This is the input data for
  
 #### Output(s): 
- _out.RDS_ This is the expected output from running _[Analysis.R](./REALDATA/Analysis.R)_ and using the input data provided. Reading this file requires to have installed CmdStanR (see section Requirements).
- _out_dignosi.txt_ When running _1._ _Analysis.R_, you will obtain a .txt file that will contain several quality control descriptors for the model. 
- _out_model.R_ When running _1._ _Analysis.R_, you will obtain a .R file with all the implemented DEB model for this analysis.
  
### Exploration: 
2. [_ExploringResult.R_](./REALDATA/ExploringResult.R) This is a simplified script to show graphically the exploration for the output of the model. This script will help you to navigate and interpret the _out.RDS_ file. 
   - *Inputs*: _input.RData_; _out.RDS_
   - *Outputs*: Table 3 in main test; Figures 1, 2, 3 and 4 in main test; Tabla S12 

# Citation
If you use this code, please consider citing our work:
Palmer, M., Moro-Martínez, I., Tomás-Ferrer, J., Orestis, Marilo, Marine, Amalia, Lika, Campos-Candela, A. (2024) Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications.


Any questions: .....@......
