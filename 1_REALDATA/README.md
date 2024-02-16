# Analyses from real data
Here we provide the code for obtaining the main content in the paper Palmer _et al._ 2024. 
By running [Analysis.R](../1_REALDATA/Analysis.R), you will estimate four DEB parameters from the provided observed data in [_input.Rdata_](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v1.0/input.RData) and obtain the output (which is also provided in the file called _out.RDS_). Afterwards, you will be able to explore the output of the model by running [_ExploringResult.R_](../1_REALDATA/ExploringResult.R). Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps.

In this section you find the following files: 

## Work flow: 
1. [Analysis.R](../1_REALDATA/Analysis.R) This script containing the **STAN code** for estimating four DEB parameters from the observed data contained in _input.Rdata_. This script has served as the skeleton of the work presented in the main text in Palmer _et al._ 2024. It can be run in multicore or in local version. As the script is presented here, it is prepared for be run in a multicore version (the  local version is on the script but currently locked). Note however that if you do not have multicore, running this script may take a considerable amount of time (in the order of several days). For this reason, we provide the output of the model as well.
   - *Input*:
      - [_input.Rdata_](../1_REALDATA/input.RData) This file contains the experimental data used for estimating the parameters of the DEB model for our case study (more datils in Palmer _et al._ 2024). 
   - *Outputs*:
      - _out.RDS_ This is the expected output from running _Analysis.R_ and using the input data provided. Reading this file requires to have installed CmdStanR (see section Requirements). Note: the size of this file exceeds the limits for github repositoru, you can downlaod from releases [here](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v1.1/out.RDS).
      - _out_dignosi.txt_ When running _Analysis.R_, you will obtain a .txt file that will contain several quality control descriptors for the model.
      - _out_model.R_ When running _Analysis.R_, you will obtain a .R file with all the implemented DEB model in **STAN code** for this analysis.
   
2. [_ExploringResult.R_](../1_REALDATA/ExploringResult.R) This is a simplified script to show graphically the exploration for the output of the model. This script will help you to navigate and interpret the _out.RDS_ file. 
   - *Inputs*:
     - _input.RData_
     - _out.RDS_ Note: the size of this file exceeds the limits for github. Note: the size of this file exceeds the limits for github repositoru, you can downlaod from releases [here](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v1.1/out.RDS).
   - *Outputs*:
        - Table 3 in main test
        - Figures 1, 2, 3 and 4 in main test
        - Tabla S12 in Supplementary materials.
