# Simulations
Data-simulation experiments were completed for demonstrating opportunities and limitations of our analytical strategy for estimating DEB parameters at the individual level. We designed simulation experiments to emulate the observed patterns in our empirical data aiming to assess the accuracy and precision of the proposed analytical framework.

In this section you find the following files:

## Working flow

1. [**Simulation.R**](../2_SIMULATIONS/Simulation.R)
  - *Inputs*:
      - [_input.Rdata_](../2_SIMULATIONS/input.RData)
      - _out.RDS_ Reading this file requires to have installed CmdStanR (see section Requirements). Note: the size of this file exceeds the limits for github. We are exploring a way to make it available.
      - [_DEB.R_](../2_SIMULATIONS/DEB.R)
  - *Outputs*:
      -  _out_sim.RDS_ Reading this file requires to have installed CmdStanR (see section Requirements). Note: the size of this file exceeds the limits for github. We are exploring a way to make it available.
      -  [_input_sim.RData_](../2_SIMULATIONS/input_sim.RData)

2. [**Visual exploration of results**](../2_SIMULATIONS/Fig5.R)
   - *Inputs*:
      -  _out_sim.RDS_ Reading this file requires to have installed CmdStanR (see section Requirements). Note: the size of this file exceeds the limits for github. We are exploring a way to make it available.
      -  [_input_sim.RData_](../2_SIMULATIONS/input_sim.RData)
    - *Outputs*:
      - Fig 5 in main test in Palmer _et al._ 2024



