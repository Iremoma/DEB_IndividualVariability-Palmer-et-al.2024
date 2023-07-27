#---------------------------------------
# Estimating DEB parameters for gilt-head bream (Sparus aurata)
#---------------------------------------
# Miquel Palmer
# Last update: 17-07-2022
# version with threads_per_chain (more than one core per chain)

#-----
# 1: Loading data and libraries
#-----
remove(list=ls())
library(cmdstanr)  # Installing cmdstanr at: https://mc-stan.org/cmdstanr/

load("input.RData")
# N                 # number of fish  
# n                 # number of sampling events
# obs               # Observations: (3D array) obs[repeated measures,fish,[length,weight]]
# t0                # t from which integrals are calculated
# ts                # times at which the system is observed (sampling events)
# day.correction    # number of days between birthday and Jan 1st (phase for sinusoidal function of temperature)
# length_method     # device used for measuring length (1: ichtiometer; 2:digital ruler))
# ID.label          # ABSA fish ID

# CAUTION! the dimension to be sliced (fish) must be the first
dim(obs)
#[1] 11 69  2
new.obs=array(NA,dim=c(N,2,n))
for(i in 1:n){
  for (j in 1:2){
    new.obs[,j,i]=obs[i,,j]
  }
}
dim(new.obs)


load("parms.RData")
# parms=NULL
# parms$Kelvin = 273.15       # Kelvin transformation
# parms$cw = 5                # water fraction of the organism (-); note that wet weight density is dV*cW = 1 gr/cm3 (assumption)
# parms$wE = 23.9             # wE=23.9 (g/mol) specific dry weight of water free reserve (assuming C:H1.8:O0.5:N0.15)
# parms$muE = 550000          # 550000  (J/mol) chemical potential of reserve energy (add-my-pet)
# parms$dv = 0.2              # (g cm-3) specific density of structure (dry weight)  
# parms$dVw = 1               # dv*cw; (g cm-3) specific density of structure (wet weight) 
# parms$T1 = 273.15 + 20      # reference temperature (20 C)
# parms$TA = 8000             # Arrhenius temperature (add-my-pet)
# parms$f = 1                 # functional response (-) Caution! (add-my-pet) 1: feed at libitum
# parms$EG = 5229             # (J/cm^3) specific cost for structure (add-my-pet)
# parms$kJ = 0.002            # (1/d) maturity maint rate coefficient (add-my-pet)
# parms$kappa = 0.9542        # (-) allocation fraction to soma (add-my-pet)
# parms$pAm =  NULL           # (J/d.cm^2) maximum assimilation flux; To be estimated
# parms$kM = 19.17            # (J/d.cm^3) vol-spec somatic maint (add-my-pet)
# parms$v = NULL              # (cm/d) energy conductance; To be estimated
# parms$rho = 0.1917          # shape coef; add-my-pet
# 
# parms$acceleration_factor= 	18.701#parms$Lj/parms$Lb # add-my-pet=6.82314
# parms$acceleration=F
# parms$pAm = 25.1729
# parms$v = 0.04372
# parms$pAm_after = parms$pAm*parms$acceleration_factor
# parms$v_after = parms$v*parms$acceleration_factor
# 
# parms$acceleration = F      # acceleration?
# parms$biomass = F           # biomass effects?
# parms$spawning = F          # reproductive dynamics?
# # parameters of the sinusoidal function for temperature
# # CAUTION! Temp_phi have been estimated assuming that t0 = "2017-01-01",
# # thus, (birthday - "2017-01-01") days must be added to t for proper
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
    // ODE Parameters
    real Kelvin,
    real Temp_mean,
    real Temp_amp,
    real pi2f,
    real Temp_phi,
    real TA,
    real T1,
    real f,
    real EG,
    real kappa,    
    real pAm, 
    real pM,
    real v, 
    real birth_correction,
    real W0,
    real W1
    ){
      // Derivatives
      // real dEdt;  // Reserve Energy (j)
      // real dVdt;  // Structural length (cm)
      vector[2] dydt;

      // Auxiliary variables
      real Temp;         // Temperature; f(t)
      real cT;           // Arrhenius temperature correction [0,1]
      real biomass;      // Biomass; f(t)
      real W;            // Biomass correction [0,1]
      real pAm_T;        // corrected (tempertaure and biomass) parameters
      real v_T;
      real pM_T;
      //real kJ_T;
      
      // fluxes (j/day)
      real pA; 
      real pS;
      real pC;
      real pG;
      
      // Temperature f(t) (parameters estimated from actual temperature data)
      Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(t + birth_correction) + Temp_phi); 
      // temperature correction (cT)
      cT = exp(TA/T1 - TA/Temp); // simple monotonic version
      
      // biomass f(t) (parameters estimated from actual biomass data)
      biomass = 3.173425e+03-3.720773e+02*t+2.362952e+00*t^2-1.431395e-03*t^3; 
      // biomass correction (W) [0,1]
      // W0 is the actual biomass at t0
      // W1 (esimated) is the biomass at which W=0 (neither assimilation nor mobilization)
      W = 1*(biomass<=W0)+
          (1+(-W0/(W0-W1))+(biomass/(W0-W1)))*(biomass>W0)*(biomass<W1)+
          0*(biomass>=W1);
      
      // corrected parameters      
      pAm_T=cT*pAm*W;             
      v_T=cT*v*W;
      pM_T=cT*pM*W;
		  
      // fluxes
      pA = f*pAm_T*pow(y[2],2.0/3.0);                                           // assimilation rate    
      pS = pM_T*y[2];			                                                      // somatic maintenance rate;			                                           // maturity maintenance rate
      pC = (y[1]/y[2])*((v_T*EG*pow(y[2],2.0/3.0)+pS)/(EG + kappa*y[1]/y[2]));  // mobilization rate
      pG = (kappa*pC - pS);//*(kappa*pC>pS);                                    // growth rate (priority for maintenance)

      // derivatives
     dydt[1]=pA-pC;
     dydt[2]=pG/EG;

     return dydt;
    }
                        
  // Function for reduce_sum (for multicore optimization)
  real partial_sum_lpdf(array [,,] real slice_obs,  // 0bservations
                        int start, int end,         // reduce_sum  paramters,
                        array [] vector y0,         // initial values
                        real t0,                    // initial time (fish specific)
                        array []real ts,            // observation times
                        int n,                      // number of observations
                        int N,                      // number of fish
                        // ODE parameters 
                        real Kelvin,
                        real Temp_mean,
                        real Temp_amp,
                        real pi2f,
                        real Temp_phi,
                        real TA,
                        real T1,
                        real f,
                        real EG,
                        real pM,
                        real birth_correction,
                        array[] real pAm,
                        real kappa,
                        array [] real v,
                        real rho,
                        array[] real sd_length_n,
                        real cw,              
                        real wE,
                        real muE,
                        real sd_weight,
                        real W0,
                        real W1
                         ) {

    array [n] vector[2] y;  // numerical integration results
    real likelihood=0;      // cummulated likelihood across fish
    for (i in 1:end-start+1){  //fish i
      y = ode_rk45_tol(DEB, y0[start+i-1], t0 , ts[],
        1e-6,1e-6,10000,// accuracy 
        // ODE parameters
        Kelvin, Temp_mean, Temp_amp, pi2f, Temp_phi, 
        TA,
        T1,
        f, EG, 
        kappa,
        pAm[start+i-1],
        pM,
        v[start+i-1],
        birth_correction,
        W0, W1
      ); 
     likelihood +=
       lognormal_lupdf(to_vector(slice_obs[i,1,1:n]) | log(pow(to_vector(y[1:n,2]),(1/3.0))/rho) ,   to_vector(sd_length_n[1:n]))
        +
       lognormal_lupdf(to_vector(slice_obs[i,2,1:n]) | log(to_vector(y[1:n,2]) + cw*wE/muE*to_vector(y[1:n,1])) , sd_weight)
       ;
    }
    return
      likelihood;
  }
}

data {
  int<lower=0> n;            // number of observations pr fish
  int<lower=0> N;            // number of fish
  array [N,2,n] real obs;    // observations; fishes,c(length,weight), repeated measures (1:n)
  array [n] int length_method;  // method used for measure length (three methods)
  real t0;                   // initial time (1 day before the first sampling date)
  array [n] real ts;         // times at which the system is observed
  real birth_correction;     // number of days between jan 1th and birthday
  int<lower=1> grainsize;    // parameter for optimizing multicore parametrization

  // paramters for the integrate function (ODE) 
  real Kelvin;      
  real Temp_mean;  
  real Temp_amp;    
  real pi2f;        
  real Temp_phi;    
  real TA;
  real T1;          
  real f;           
  real EG;         
  real kappa;     
  real pM;

  // other parameters
  real rho;        
  real wE;        
  real muE;        
  real cw;
  real W0;
  
  // priors
  real pAm_prior;
  real v_prior;
}

transformed data {
}

parameters {
  // DEB parameters
  real <lower=0> pAm_mu;
  real <lower=0> pAm_sd;
  array[N] real <lower=0> pAm;
  real <lower=0> v_mu;
  real <lower=0> v_sd;
  array[N] real<lower=0> v;
  
  // initial values for state variables
  real <lower=0>y0_1mu;
  real <lower=0>y0_1sd;
  real <lower=0>y0_2mu;
  real <lower=0>y0_2sd;
  array [N] vector<lower=0>[2] y0; 
  
  // observation error
  array[3] real<lower=0> sd_length;
  real <lower=0> sd_weight;
  
  // effects of crowding
  real <lower=W0> W1;
}

transformed parameters {
// setting sd_length for each sampling date
  array[n]  real<lower=0> sd_length_n;
  for (i in 1:n){
    sd_length_n[i]= sd_length[length_method[i]];
  }
}   

model {
  // priors
  // pAm
  pAm_mu ~ normal(pAm_prior,pAm_prior);
  pAm_sd ~ normal(pAm_prior*0.2,pAm_prior*0.2);
  pAm ~ normal(pAm_mu,pAm_sd);
  
  // v
  v_mu ~ normal(v_prior,v_prior);
  v_sd ~ normal(v_prior*0.2,v_prior*0.2);
  v ~ normal(v_mu,v_sd);
  
  // initial values for state variables
  y0_1mu ~ normal(1.0e6,1.0e6);     // reserve energy
  y0_1sd ~ normal(1.0e5,1.0e5);
  y0_2mu ~ normal(9.49186,9.49186); // structural volume 
  y0_2sd ~ normal(2.22539,2.22539);
  for (i in 1:N){
    y0[i,1] ~ normal(y0_1mu,y0_1sd);
    y0[i,2] ~ normal(y0_2mu,y0_2sd);
  }
  
  // water quality
  W1 ~ normal(500000.0,300000.0);
  
  // error
  sd_length ~ lognormal(-3,1);  //several values
  sd_weight ~ lognormal(-3,1);
  
  // likelihood
  target += reduce_sum(partial_sum_lupdf,
                        obs[,,],   // array to be slice (multi threads_per_chain)
                        grainsize, // slicing grainsize (multi threads_per_chain)
                        // arguments for ode_rk45 (DEB model)
                        y0[], t0, ts[], n, N, // state variables at t0, observation dates, number of obsrvations per fish, number of fish
                        Kelvin, Temp_mean, Temp_amp, pi2f, Temp_phi,  // temperature parameters
                        TA,T1, // Monototic temperature effect
                        f, EG,pM,birth_correction,
                        pAm[],  // one value per fish
                        kappa,
                        v[],    // one value per fish
                        rho,
                        sd_length_n[], // three methods
                        cw, wE, muE, sd_weight,W0,W1
                        );

}

generated quantities {
  array [n,N] vector[2] y_hat;
  array [n,N] vector[2] obs_hat;   //model prediction
  
  for(i in 1:N){
   y_hat [1:n,i] = ode_rk45_tol(DEB, y0[i], t0 , ts[1:n],
      1e-6,1e-6,10000, // accuracy numerical integration
      // ODE parameters
      Kelvin, Temp_mean, Temp_amp, pi2f, Temp_phi, TA, T1,
      f, EG, kappa, pAm[i], pM, v[i], birth_correction, W0,W1
    );
    for (j in 1:n){
      obs_hat[j,i,1] =  pow(y_hat[j,i,2],(1.0/3.0))/rho;
      obs_hat[j,i,2] =  y_hat[j,i,2] + cw*wE/muE*y_hat[j,i,1];
    }
  }
}

" # end of model code
,fill = TRUE)
sink()

#-----
# 4: compiling model
#-----
mod = cmdstan_model("DEB.stan",cpp_options = list(stan_threads = TRUE))

#-----
# 5: initializating the model
#-----
n.chains = 3
initializer = function() list(
        "pAm_mu"=parms$pAm_after,
        "pAm_sd"=10.0,
        "pAm"=rep(parms$pAm_after,N),
        "v_mu"=parms$v_after,
        "v_sd"=0.001,
        "v"=rep(0.05,N),
        "y0_1mu"=1e6,
        "y0_1sd"=1e6,
        "y0_2mu"=10,
        "y0_2sd"=10,
        "y0"=t(c(1e6,10)%*%t(rep(1,N))),
        "sd_length"=c(0.01,0.01,0.01),
        "sd_weight"=0.05,
        "W1"=500000
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()

#-----
# 6: running
#-----
fit = mod$sample(
        data =list (
                N = N,                                # number of fish
                n = n,                                # number of observations per fish (number sampling events)
                obs = new.obs,                            # observations (dim=74*13*2)
                t0 = ts[1]-1,                              # initial time
                ts = ts,                              # times at which the system is observed (sampling dates)
                birth_correction = day.correction,    # days between 1 jan and birthday
                grainsize = 1, #28                       # technical parameters for multiple cores
                #paramters
                Kelvin = parms$Kelvin,                # Kelvin zero
                Temp_mean = parms$Temp_mean,          # mean temperature (sinusoidal fit)
                Temp_amp = parms$Temp_amp,            # amplitude temperature (sinusoidal fit)
                pi2f = parms$pi2f,                    # 2pi/365
                Temp_phi = parms$Temp_phi,            # phase temperature (sinusoidal fit)
                TA = parms$TA,                        # Arhhenius temperaure
                T1 = parms$T1,                        # reference temperture
                f = parms$f,                          # functional response
                EG = parms$EG,                        # specific cost for structure
                kappa = parms$kappa,                  # allocation fraction to soma
                pM = parms$kM,                        # vol-spec somatic maint
                rho = parms$rho,                      # shape coefficient
                wE = parms$wE,                        # specific dry weight of reserve (note that density of wet structure is assumed to be 1.0 but it is not included)
                muE = parms$muE,                      # chemical potential of reserve
                cw = parms$cw,                        # water fraction
                pAm_prior=parms$pAm_after,            # mean prior for pAm
                v_prior=parms$v_after,                # mean prior for v
                length_method = length_method,        # measurement method for length
                W0=27652.82                           # biomass at t0 (gr)
                ),
        chains = n.chains,
        parallel_chains = n.chains,
        threads_per_chain = 6,
        iter_warmup = 3000,
        iter_sampling = 3000,
        init = inits,
        max_treedepth = 15,
        adapt_delta = 0.99
)

#-----
# 8: saving
#-----
fit$cmdstan_diagnose()
#fit$save_object(file = "parallel.RDS")


