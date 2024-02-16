# Simulations
Data-simulation experiments were completed for demonstrating opportunities and limitations of our analytical strategy for estimating DEB parameters at the individual level. We designed simulation experiments to emulate the observed patterns in our empirical data aiming to assess the accuracy and precision of the proposed analytical framework.

In this section you find the following files:

## Working flow

1. [**Simulation.R**](../2_SIMULATIONS/Simulation.R)
  - *Inputs*:
      - [_input.Rdata_](../2_SIMULATIONS/input.RData)  You can also download it directly from releases [clicking here](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v1.0/input.RData).
      - _out.RDS_ This is the output from running _Analysis.R_ and using the input data provided. Reading this file requires to have installed CmdStanR (see section Requirements). Note: the size of this file exceeds the limits for github repository, you can download it directly from releases [clicking here](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v1.1/out.RDS).
      - [_DEB.R_](../2_SIMULATIONS/DEB.R) Here you will find the DEB model function we used to simulate data. The code is written in STAN.
  - *Output/Input*:
      -  [_input_sim.RData_](../2_SIMULATIONS/input_sim.RData) This is the simulated data we obtained from the DEB model in _DEB.R_ and with "known" values for the DEB parameters at the individual level (this data is simulated in the first half of the script _Simulations.R_). DEB "known" values were based on the predictions obtained in the [first section](../1_REALDATA). Simulated data will be used to estimate DEB parameters again (in the second half of the script _Simulations.R_). You can also download it directly from releases [clicking here](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v2.0/input_sim.RData)
  - *Output*:
      -  _out_sim.RDS_ This is the output of our routing for estimating DEB parameters at the individual level from the simulated data. Reading this file requires to have installed CmdStanR (see section Requirements). Note: the size of this file exceeds the limits for github repository, you can download it directly from releases [clicking here](https://github.com/Iremoma/DEB_IndividualVariability-Palmer-et-al.2024/releases/download/v2.1/out_sim.RDS)


2. [**Visual exploration of results**](../2_SIMULATIONS/Fig5.R)
   - *Inputs*:
      -  _input_sim.RData_ Simulated data. More details above. 
      -  _out_sim.RDS_ Predictions for DEB parameters at the individual level from the simulated data. More details above.    
    - *Outputs*:
      - Fig 5 in main test in Palmer _et al._ 2024 This is a visualization to compare the prediction of the model with the simulated data and "known" parameters from which the simulated data was obtained.



