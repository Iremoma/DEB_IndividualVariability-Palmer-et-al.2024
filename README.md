# DEB_IIV_STAN  Palmer-et-al.-2024
This is the code related to the publication: Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications.
Authors: Miquel palmer ....
Prior to use the code:

# Contents
In this repository we provide all the scripts necessary to reproduce the data exploration and analysis, simulations, and results presented on the paper Palmer et al., 2024. Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications.

We provide the STAN routine that we have developed for estimating four DEB parameters at the individual level for 69 cultured Gilt-head breams (Sparus aurata) for which up to 11 measures   of length and wet weight were available. 
For this, all the code is organized in the following folders: 
1. Analyses
2. Simulation

SEGUIR LLISTA 
# Requirements
Install: 
- R (4.1.3)
- CmdStanR (Command Stan R) in R: CmdStanR is a lightweight interface to Stan for R users. If you are using Stan in R for the first time, or you have no experience with this package, we encourage you to visit [Getting started with CmdStanR]( https://mc-stan.org/cmdstanr/articles/cmdstanr.html)  and follow the recommended steps to install the package.

Recommended properties of the computer, and those that characterized the computer used for running these models, are: 
- Intel(R) Xeon(R)  CPU E5-2620 v4 @ 2.10GHz (40 cores and 64 GB RAM)
# Working flow

## Analyses
Here we provide the main script of the paper Palmer _et al._ 2024. By running it (_1_Analysis.R_), you will estimate four DEB parameters from the provided observed data (_input.Rdata_) and obtain the output (which is also provided in _out.RDS_) from which all the main content of the paper and figures have been obtained (_2_ExploringResult.R_). Note that the script is self-explanatory and has a lot of notations to guide the user on the different steps. Here you find the following files: 
### Input data: 
_input.Rdata_ This file contains the experimental data used for estimating the parameters of the DEB model for our case study (more datils in Palmer _et al._ 2024). This is the input data for _1_Analysis.R_.
### Core script: 
_1_Analysis.R_ This is the main script, containing the code for estimating four DEB parameters from the observed data contained in input.RData. This script has served as the skeleton of the project in Palmer et al. 2024. It can be run in multicore or local version. As the script is presented here, it is prepared for be run in a multicore version (the  local version is on the script but currently locked). Note however that if you do not have multicore, running this script may take a considerable amount of time (in the order of several days). For this reason, we provide the output of the model as well.
### Output(s): 
_out.RDS_ This is the expected output from running _1_Analysis.R_ and using the input data provided. Reading this file requires to have installed CmdStanR (see requirements above).
_out_dignosi.txt_ When running _1_Analysis.R_, you will obtain a .txt file that will contain several quality control descriptors for the model. 
_out_model.R_ When running _1_Analysis.R_, you will obtain a .R file with all the implemented DEB model for this analysis.
### Out Exploration: 
_2_ExploringResult.R_  This is a simplified script to show graphically the exploration for the output of the model. This script will help you to navigate and interpret the _out.RDS_ file. 

REPETIRIA EL MATEIX PER LA RESTA DE CARPETES PENDENTS AMB CADASCUN DELS ANALISIS REALITZATS.

# Citation
If you use this code, please consider citing our work:
Palmer, M., Moro-Martínez, I., Tomás-Ferrer, J., Orestis, Marilo, Marine, Amalia, Lika, Campos-Candela, A. (2024) Assessing between-individual variability in bioenergetics modelling: opportunities, challenges, and potential applications.


Any questions: .....@......
