#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 7
#---------------------------------------
# Last update: 09 july 2024
# NON PARALLELIZED VERSION

#-----
# 1: Loading data and libraries
#-----

remove(list=ls())
library(cmdstanr)
load("input.RData")

# N                 # Number of fish  
# n                 # Number of sampling events
# obs               # Observations: (3D array) obs[fish,[length,weight],repeated measures]
# obs0              # Observations at t0: (2D array) obs0[fish,[length,weight]]
# t0                # Age from which integrals are calculated
# age               # Ages at which the system is observed (sampling events)
# birth.correction  # Number of days between birthday and Jan 1st (phase for sinusoidal function of temperature)
# length_method     # Device used for measuring length (1: ictiometer; 2:digital ruler)
# ID.label          # fish ID
# length_precision  # Precision (depends on device and sampling event
# weight_precision  # Precision (note that precision for the first sampling event are not included becasue they are the same than for the second sampling event)
# parms             # DEB parameters


#-------------------------------------
# "pars_init_sparus_aurata.m" (DEB parameters file at add-my-pet; accessed 21 mar 2023)
#-------------------------------------
# %% Core primary parameters: 
# par.z = 1.253;        free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
# par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
# par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
# par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
# par.v = 0.04372;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
# par.kap = 0.9542;     free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
# par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
# par.p_M = 19.17;      free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
# par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
# par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
# par.E_G = 5229;       free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
# par.E_Hb = 0.04368;   free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
# par.E_Hj = 290.8;     free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
# par.E_Hp = 145400;    free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
# par.h_a = 1.885e-09;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
# par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 
# 
# %% Other parameters: 
# par.E_Hh = 0.01242;   free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
# par.T_A = 8000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
# par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 
# par.del_M = 0.1917;   free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient'; 
# par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
# par.f_tL = 1.2;       free.f_tL  = 0;   units.f_tL = '-';         label.f_tL = 'scaled functional response for 1-var data'; 

#-------------------------------------
# DEB paramters (parms)
#-------------------------------------
# parms: DEB parameters
# parms$Kelvin = 273.15        # Kelvin transformation
# parms$w_E = 23.9             # (g/mol) specific dry weight of water free reserve (assuming C:H1.8:O0.5:N0.15)
# parms$mu_E = 550000          # (J/mol) chemical potential of reserve energy
# parms$d_V = 0.2              # (gr/cm^3) specific density of reserve (reference is water => d_V can be used to transform dry weight into wet weight)
# parms$T1 = 273.15 + 20       # reference temperature (20 C)
# parms$TA = 8000              # Arrhenius temperature 
# parms$f = 1                  # functional response (-) 1: feed at libitum
# parms$E_G = 5229             # (J/cm^3) specific cost for structure 
# parms$k_J = 0.002            # (1/d) maturity maint rate coefficient 
# parms$kap = 0.9542           # (-) allocation fraction to soma 
# parms$pAm = 25.1729          # (J/d/cm^2) maximum assimilation flux (BEFORE acceleration)
# parms$p_M = 19.17            # (J/d/cm^3) vol-spec somatic maint
# parms$v = 0.04372            # (cm/d) energy conductance (BEFORE acceleration)
# parms$rho = 0.1917           # (-) shape coef
# parms$s_M = 18.701 	         # (-) acceleration factor at f=1 and 17.5 degrees
# parms$pAm_after = parms$pAm*parms$s_M
# parms$v_after = parms$v*parms$s_M

# # Parameters of the sinusoidal function for temperature (Supplementary materials 1):
# # Temp_phi have been estimated assuming t0 = "20XX-01-01",
# # thus, (birthday - "20XX-01-01") days must be added to t for proper
# # estimation of the fish-specific temperature profile
# parms$Temp_mean = 18.72
# parms$Temp_amp = 1.568635
# parms$Temp_phi = 4.18562
# parms$pi2f = 1.721421e-02

#-----       
# 2: STAN model
#-----
sink("DEB.stan")
cat(" // first line
functions {
  // DEB model (ODE)
    vector DEB(
    real t,          // time
    vector y,        // state variables
    // core parameters
    real f,
    real E_G,
    real kap,    
    real pAm, 
    real p_M,
    real v,
    
    // temperature parameters
    real Kelvin,
    real Temp_mean,
    real Temp_amp,
    real pi2f,
    real Temp_phi,
    real TA,
    real T1,
    real birth_correction,
    
    // crowding parameters
    real W0,
    real W1
    ){
      // Derivatives
      vector[2] dydt;  // 1: dEdt; 2: dVdt;

      // Auxiliary variables
      // fluxes (j/day)
      real pA; 
      real pS;
      real pC;
      real pG;
      
      // auxiliary variables and parameters
      real Temp;         // Temperature
      real cT;           // Arrhenius temperature correction
      real sT1;
      real sT;
      real pAm_T;        // temperature corrected parameters
      real v_T;
      real p_M_T;
      real biomass;      // Biomass; f(t)
      real W;            // Biomass correction [0,1]
      
      // temperature correction
      Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(t + birth_correction) + Temp_phi);
      cT = exp(TA/T1 - TA/Temp);
      
      // crowding (W) [0,1]
      biomass = 3.173425e+03-3.720773e+02*t+2.362952e+00*t^2-1.431395e-03*t^3;
      // biomass f(t) (parameters estimated from the actual biomass data) 
      // W0 is the actual biomass at t0
      // W1 (esimated) is the biomass at which W=0 (no assimilation)
      W = 1*(biomass<=W0)+
          (1+(-W0/(W0-W1))+(biomass/(W0-W1)))*(biomass>W0)*(biomass<W1)+
          0*(biomass>=W1);
          
      // Corrected parameters  (crowding effects on pAm only)    
      pAm_T = cT*pAm*W;             
      v_T = cT*v;
      p_M_T = cT*p_M;
      
      // fluxes (j/day)
      pA = f*pAm_T*y[2]^(2.0/3.0);                                                // assimilation rate    
      pS = p_M_T*y[2];			                                                      // somatic maintenance rate;			                                           // maturity maintenance rate
      pC = (y[1]/y[2])*((v_T*E_G*pow(y[2],2.0/3.0)+pS)/(E_G + kap*y[1]/y[2]));    // mobilization rate
      pG = (kap*pC - pS);                                                         // growth rate (priority for maintenance)

      // derivatives
      dydt[1]=pA-pC;
      dydt[2]=pG/E_G;
      
      return dydt;
    }
}

data {
  int n;                             // number of observations per fish
  int N;                             // number of fish
  array [N,2,n] real obs;            // observations; fishes,c(length,weight), repeated measures (1:n)
  array [N,2] real obs0;             // observations at t0; fishes,c(length,weight)
  real t0;                           // initial time
  array [n] real age;                // times at which the system is observed
  array [n] real weight_precision;   // sampling event-specific precision for weigth
  array [n] real length_precision;   // sampling event-specific precision for length
  array [n] int length_method;       // code for length methods (1: ictiometer; 2 electronic device)
  
  // DEB parameters asumed be be known (not estimated) 
  real f;           
  real E_G;         
  real kap;     
  real p_M;

  // other parameters
  real del_M;        
  real w_E;        
  real mu_E;        
  real d_V;
  real W0;
  
  // temperature parameters
  real Kelvin;
  real Temp_mean;
  real Temp_amp;
  real pi2f;
  real Temp_phi;
  real TA;
  real T1;
  real birth_correction;
  
  // priors
  array [2] real prior_E0_mu;
  array [2] real prior_E0_sd;
  array [2] real prior_V0_mu;
  array [2] real prior_V0_sd;
  array [2] real prior_pAm;
  array [2] real prior_pAm_sd;
  array [2] real prior_v;
  array [2] real prior_v_sd;
  array [2] real prior_sd_length;
  array [2] real prior_sd_weight;
  array [2] real prior_W1;
  
  // for parallelization
  int grainsize;    // parameter for optimizing multicore parametrization
}

transformed data {
}

parameters {  // DEB parameters
  real <lower=0> pAm_mu;
  real <lower=0> pAm_sd;
  array [N] real <lower=0> pAm;

  real <lower=0> v_mu;
  real <lower=0> v_sd;
  array [N] real <lower=0> v;
  
  // initial state
  real <lower=0> E0_mu;
  real <lower=0> E0_sd;
  real <lower=0> V0_mu;
  real <lower=0> V0_sd;
  array [N] vector <lower=0> [2] y0;
  
  // observation error
  array[2] real <lower=0> sd_length;
  real <lower=0> sd_weight;
   
  // effects of crowding
  real <lower=W0> W1;
}

transformed parameters {
  // setting sampling event specific error
  // Accounting for differnet methods and different precisions
  array[n]  real<lower=0> sd_length_n;
  array[n]  real<lower=0> sd_weight_n;
  for (i in 1:n){
    sd_length_n[i] = sqrt (sd_length[length_method[i]]^2 + 1.0/12.0*length_precision[i]^2);
    sd_weight_n[i] = sqrt (sd_weight^2 + 1.0/12.0*weight_precision[i]^2);
  }
}   

model {
  array [N,n] vector[2] y;  // template for the numerical integration results

  // pAm
  pAm_mu ~ normal(prior_pAm[1] , prior_pAm[2]);
  pAm_sd ~ cauchy(prior_pAm_sd[1] , prior_pAm_sd[2]);
  pAm ~ normal(pAm_mu , pAm_sd);
  
  // v 
  v_mu ~ normal(prior_v[1] , prior_v[2]);
  v_sd ~ normal(prior_v_sd[1] , prior_v_sd[2]);
  v ~ normal(v_mu , v_sd);
  
  // error
  sd_length ~ cauchy(prior_sd_length[1],prior_sd_length[2]);
  sd_weight ~ cauchy(prior_sd_weight[1],prior_sd_weight[2]);
  
  // crowding
  W1 ~ normal(prior_W1[1],prior_W1[2]);
    
  // Initial values for state variables
  E0_mu ~ normal(prior_E0_mu[1], prior_E0_mu[2]);
  E0_sd ~ normal(prior_E0_sd[1], prior_E0_sd[2]);
  V0_mu ~ normal(prior_V0_mu[1], prior_V0_mu[2]);
  V0_sd ~ normal(prior_V0_sd[1], prior_V0_sd[2]);
  for (i in 1:N){  // fish level
    y0[i,1] ~ normal(E0_mu, E0_sd);
    y0[i,2] ~ normal(V0_mu, V0_sd);
    obs0[i,1] ~ normal( ((y0[i,2]^(1.0/3.0))/del_M), sd_length_n[1]);        // length at t0
    obs0[i,2] ~ normal( (y0[i,2] + y0[i,1]*w_E/(mu_E*d_V)), sd_weight_n[1]); // weight at t0
  }
  
  // Likelihood
  for (i in 1:N){  // fish level
  // Expected values of the state variables
    y[i,1:n] = ode_rk45(DEB, y0[i], t0 , age,
                // DEB main parameters
                f,E_G,kap,pAm[i],p_M,v[i],
                // temperature-related paramters
                Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction,
                //crowding
                W0,W1
    );
    for (j in 1:n){ // Observation level
      target +=  // likelihood
         normal_lupdf( obs[i,1,j] | pow(y[i,j,2],(1.0/3.0))/del_M , sd_length_n[j] )
         +
         normal_lupdf( obs[i,2,j] | y[i,j,2] + w_E/(mu_E*d_V)*y[i,j,1] , sd_weight_n[j] )
         ;
    }
  }
}

generated quantities {
  array [N,n] vector[2] y_hat;
  array [N,2,n] real obs_hat;   //model prediction
  array [N,2] real obs0_hat;    //model prediction

  for(i in 1:N){
    y_hat[i,1:n] = ode_rk45(DEB, y0[i], t0 , age,
          // DEB main parameters
          f,E_G,kap,pAm[i],p_M,v[i],
          // temperature-related paramters
          Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction,
          // crowding
          W0,W1
          );
    for (j in 1:n){
      obs_hat[i,1,j] =  pow(y_hat[i,j,2],(1.0/3.0))/del_M;
      obs_hat[i,2,j] =  y_hat[i,j,2] + w_E/(mu_E*d_V)*y_hat[i,j,1];
    }
    obs0_hat[i,1] =  pow(y0[i,2],(1.0/3.0))/del_M;
    obs0_hat[i,2] =  y0[i,2] + y0[i,1]*w_E/(mu_E*d_V);
  }
}

" # End of model code
,fill = TRUE)
sink()

#-----
# 3: Compiling model
#-----
#--- Local version ---

mod = cmdstan_model("DEB.stan")

#--- Multicore version ---

#mod = cmdstan_model("DEB.stan",cpp_options = list(stan_threads = TRUE))

#-----
# 4: Initializing the model
#-----
# Setting priors (Supplementary materials 9)
prior_E0_mu=c(50000,500000)# normal
prior_E0_sd=c(20000,200000)# normal
prior_V0_mu=c(30,100)# normal
prior_V0_sd=c(10,100)# normal
prior_pAm = c(parms$pAm_after , parms$pAm_after) # normal
prior_pAm_sd = c(50,100) # cauchy
prior_v = c(0, 2) # cauchy
prior_v_sd = c(0, 1) # cauchy
prior_sd_length=c(0,10) # cauchy
prior_sd_weight=c(0,100) # cauchy

# Crowding effect (Supplementary materials 4)
t=seq(t0,age[n])
biomass = 3.173425e+03-3.720773e+02*t+2.362952e+00*t^2-1.431395e-03*t^3; 
prior_W1=c( mean(biomass), 10*mean(biomass)) # normal

# Initial values
n.chains = 4
initializer = function() list(
  "pAm_mu"=parms$pAm_after,
  "pAm_sd"=50,
  "pAm"=rep(parms$pAm_after,N),
  "v_mu"=parms$v_after,
  "v_sd"=0.01,
  "v"=rep(parms$v_after,N),
  "E0_mu"=prior_E0_mu[1],
  "E0_sd"=prior_E0_sd[1],
  "V0_mu"=prior_V0_mu[1],
  "V0_sd"=prior_V0_sd[1],
  "y0"=cbind(rep(prior_E0_mu[1],N),rep(prior_V0_mu[1],N)),
  "sd_length"=c(1,1),
  "sd_weight"=10,
  "W1"=as.numeric(quantile(biomass,0.9))
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()

#-----
# 5: Running
#-----

fit = mod$sample(
  data =list (
    N = N,                                
    n = n,                                
    obs = obs,                            
    obs0 = obs0,
    t0 = t0,                              
    age = age,                            
    length_method = length_method,        
    length_precision=length_precision,
    weight_precision=length_precision,
    
    # DEB main parameters
    f = parms$f, 
    k_J = parms$k_J,
    kap_R = parms$kap_R,
    E_G = parms$E_G,
    kap = parms$kap,
    p_M = parms$p_M,
    
    # DEB secondary parameters
    del_M = 0.2805034,   #Supplementary materials 6
    w_E = parms$w_E,        
    mu_E = parms$mu_E,
    d_V = parms$d_V,
    
    # Temperature-related parameters (Supplementary materials 1)
    Kelvin = parms$Kelvin,
    Temp_mean = parms$Temp_mean,
    Temp_amp = parms$Temp_amp,
    pi2f = parms$pi2f,
    Temp_phi = parms$Temp_phi,
    TA = parms$TA,
    T1 = parms$T1,
    birth_correction = birth_correction,
    
    # Crowding (Supplementary materials 4)
    W0=27652.82, # Biomass at t0 (gr) 
    
    # Priors
    prior_E0_mu=prior_E0_mu,
    prior_E0_sd=prior_E0_sd,
    prior_V0_mu=prior_V0_mu,
    prior_V0_sd=prior_V0_sd,
    prior_pAm=prior_pAm,
    prior_pAm_sd=prior_pAm_sd,
    prior_v=prior_v,
    prior_v_sd=prior_v_sd,
    prior_sd_length=prior_sd_length,
    prior_sd_weight=prior_sd_weight,
    prior_W1=prior_W1,
    grainsize=1
  ),
  chains = n.chains,
  parallel_chains = n.chains,
  #threads_per_chain = 18,
  iter_warmup = 2000,
  iter_sampling = 2000,
  init = inits,
  max_treedepth = 12,
  adapt_delta = 0.99
)

#-----
# 6: Saving
#-----
fit$save_object(file = "out.RDS")
#sink(file = "out_dignosi.txt")
#fit$cmdstan_diagnose()
#sink()
#file.rename("DEB.stan","out_model.R")

