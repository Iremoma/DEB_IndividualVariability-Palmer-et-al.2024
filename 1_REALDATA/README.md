# Analyses from real data
Here we provide the main code for the paper Palmer _et al._ 2024. By running [Analysis.R](./REALDATA/Analysis.R), you will estimate four DEB parameters from the provided observed data in [_input.Rdata_](./REALDATA/input.Rdata) and obtain the output (which is also provided in the file called _out.RDS_) from which all the main content of the paper and figures have been obtained by running [_ExploringResult.R_](./REALDATA/ExploringResult.R). 
Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps. 
In this section you find the following files: 

## Core script: 
1. [Analysis.R](./REALDATA/Analysis.R) This script containing the code for estimating four DEB parameters from the observed data contained in _input.Rdata_. This script has served as the skeleton of the work presented in the main text in Palmer _et al._ 2024. It can be run in multicore or in local version. As the script is presented here, it is prepared for be run in a multicore version (the  local version is on the script but currently locked). Note however that if you do not have multicore, running this script may take a considerable amount of time (in the order of several days). For this reason, we provide the output of the model as well.
   - *Input*: _input.RData_
   - *Outputs*: _out.RDS_; _out_dignosi.txt_;_out_model.R_
 ### Input: 
 - [_input.Rdata_](./REALDATA/input.Rdata) This file contains the experimental data used for estimating the parameters of the DEB model for our case study (more datils in Palmer _et al._ 2024). 
 ### Output(s): 
- _out.RDS_ This is the expected output from running _Analysis.R_ and using the input data provided. Reading this file requires to have installed CmdStanR (see section [Requirements](./)).
- _out_dignosi.txt_ When running _1._ _Analysis.R_, you will obtain a .txt file that will contain several quality control descriptors for the model. 
- _out_model.R_ When running _1._ _Analysis.R_, you will obtain a .R file with all the implemented DEB model for this analysis.
  
## Exploration: 
2. [_ExploringResult.R_](./REALDATA/ExploringResult.R) This is a simplified script to show graphically the exploration for the output of the model. This script will help you to navigate and interpret the _out.RDS_ file. 
   - *Inputs*: _input.RData_; _out.RDS_
   - *Outputs*: Table 3 in main test; Figures 1, 2, 3 and 4 in main test; Tabla S12 in Supplementary materials.