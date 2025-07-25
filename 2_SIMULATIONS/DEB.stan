 // first line
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
    real birth_correction//,
    
    // crowding parameters
    //real W0,
    //real W1
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
      //real biomass;      // Biomass; f(t)
      //real W;            // Biomass correction [0,1]
      
      // temperature correction
      Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(t + birth_correction) + Temp_phi);
      cT = exp(TA/T1 - TA/Temp);
      
      // crowding (W) [0,1]
      //biomass = 3.173425e+03-3.720773e+02*t+2.362952e+00*t^2-1.431395e-03*t^3;
      // biomass f(t) (parameters estimated from the actual biomass data) 
      // W0 is the actual biomass at t0
      // W1 (esimated) is the biomass at which W=0 (no assimilation)
      //W = 1*(biomass<=W0)+
      //    (1+(-W0/(W0-W1))+(biomass/(W0-W1)))*(biomass>W0)*(biomass<W1)+
      //    0*(biomass>=W1);
          
      // Corrected parameters  (crowding effects on pAm only)    
      pAm_T = cT*pAm;//*W;             
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
    
  // Function for reduce_sum (for multicore optimization)
  real partial_sum_lpdf(array [,,] real slice_obs,  // 0bservations
                        int start, int end,         // reduce_sum  paramters,
                        // arguments for ode_rk45 (DEB model)
                        // state variables at t0, observation dates, number of obsrvations per fish, number of fish
                        array [] vector y0, real t0, array [] real age, int n, int N, 
                        // DEB main parameters
                        real f, real E_G, real kap, array [] real pAm, real p_M, array [] real v,
                        // temperature-related paramters
                        real Kelvin, real Temp_mean, real Temp_amp, real pi2f, real Temp_phi, 
                        real TA,real T1, real birth_correction,
                        // linking state variables to obervable variable
                        real del_M, real d_V, real w_E, real mu_E,
                        // error
                        //array[] real sd_length_n,
                        //array[] real sd_weight_n,
                        real sd_length,
                        real sd_weight
                        // crowding effects
                        //real W0,real W1
                        ) {

    array [n] vector[2] y;  // numerical integration results
    real likelihood=0;      // cummulated likelihood across fish

    for (i in 1:end-start+1){  // (fish i; start-to-end defines the especific fish managed by a given core)
      y = ode_rk45(DEB, y0[start+i-1], t0 , age,
                  // DEB main parameters
                  f,E_G,kap,pAm[start+i-1],p_M,v[start+i-1],
                  // temperature-related paramters
                  Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction//,
                  //crowding
                  //W0,W1
      );
      
     
     likelihood +=
       //normal_lupdf(to_vector(slice_obs[i,1,1:n]) | (pow(to_vector(y[1:n,2]),(1.0/3.0))/del_M) , to_vector(sd_length_n[1:n]) )
       normal_lupdf(to_vector(slice_obs[i,1,1:n]) | (pow(to_vector(y[1:n,2]),(1.0/3.0))/del_M) , sd_length )
       +
       //normal_lupdf(to_vector(slice_obs[i,2,1:n]) | (to_vector(y[1:n,2]) + w_E/(mu_E*d_V)*to_vector(y[1:n,1])) , to_vector(sd_weight_n[1:n]) )
       normal_lupdf(to_vector(slice_obs[i,2,1:n]) | (to_vector(y[1:n,2]) + w_E/(mu_E*d_V)*to_vector(y[1:n,1])) , sd_weight )
       ;
    }
    return
      likelihood;
  }
}

data {
  int n;                             // number of observations per fish
  int N;                             // number of fish
  array [N,2,n] real obs;            // observations; fishes,c(length,weight), repeated measures (1:n)
  array [N,2] real obs0;             // observations at t0; fishes,c(length,weight)
  real t0;                           // initial time
  array [n] real age;                // times at which the system is observed
  //array [n] real weight_precision;   // sampling event-specific precision for weigth
  //array [n] real length_precision;   // sampling event-specific precision for length
  //array [n] int length_method;       // code for length methods (1: ictiometer; 2 electronic device)
  
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
  //real W0;
  
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
  //array [2] real prior_W1;
  
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
  //array[2] real <lower=0> sd_length;
  real <lower=0> sd_length;
  real <lower=0> sd_weight;
   
  // effects of crowding
  //real <lower=W0> W1;
}

transformed parameters {
  // setting sampling event specific error
  // Accounting for differnet methods and different precisions
  //array[n]  real<lower=0> sd_length_n;
  //array[n]  real<lower=0> sd_weight_n;
  //for (i in 1:n){
  //  sd_length_n[i] = sqrt (sd_length[length_method[i]]^2 + 1.0/12.0*length_precision[i]^2);
  //  sd_weight_n[i] = sqrt (sd_weight^2 + 1.0/12.0*weight_precision[i]^2);
  //}
}   

model {
  // pAm
  pAm_mu ~ normal(prior_pAm[1] , prior_pAm[2]);
  pAm_sd ~ cauchy(prior_pAm_sd[1] , prior_pAm_sd[2]);
  pAm ~ normal(pAm_mu , pAm_sd);
  
  // v 
  v_mu ~ cauchy(prior_v[1] , prior_v[2]);
  v_sd ~ cauchy(prior_v_sd[1] , prior_v_sd[2]);
  //v ~ gamma(v_mu^2/v_sd^2 , v_mu/v_sd^2);
  v ~ normal(v_mu , v_sd);
  
  // error
  sd_length ~ cauchy(prior_sd_length[1],prior_sd_length[2]);
  sd_weight ~ cauchy(prior_sd_weight[1],prior_sd_weight[2]);
  
  // crowding
  //W1 ~ normal(prior_W1[1],prior_W1[2]);
    
  // Initial values for state variables
  E0_mu ~ normal(prior_E0_mu[1], prior_E0_mu[2]);
  E0_sd ~ normal(prior_E0_sd[1], prior_E0_sd[2]);
  V0_mu ~ normal(prior_V0_mu[1], prior_V0_mu[2]);
  V0_sd ~ normal(prior_V0_sd[1], prior_V0_sd[2]);
  for (i in 1:N){
    y0[i,1] ~ normal(E0_mu, E0_sd);
    y0[i,2] ~ normal(V0_mu, V0_sd);
    //obs0[i,1] ~ normal( ((y0[i,2]^(1.0/3.0))/del_M), sd_length_n[1]);        // length at t0
    //obs0[i,2] ~ normal( (y0[i,2] + y0[i,1]*w_E/(mu_E*d_V)), sd_weight_n[1]); // weight at t0
    obs0[i,1] ~ normal( ((y0[i,2]^(1.0/3.0))/del_M), sd_length);        // length at t0
    obs0[i,2] ~ normal( (y0[i,2] + y0[i,1]*w_E/(mu_E*d_V)), sd_weight); // weight at t0
  }

  // likelihood
  target += reduce_sum(partial_sum_lupdf,
                        obs[,,],   // array to be slice (multi threads_per_chain)
                        grainsize, // slicing grainsize (multi threads_per_chain)
                        // arguments for ode_rk45 (DEB model)
                        y0[], t0, age[], n, N, // state variables at t0, observation dates, number of obsrvations per fish, number of fish
                        // DEB main parameters
                        f,E_G,kap,pAm[],p_M,v[],
                        // temperature-related paramters
                        Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction,
                        // state variables to obervable variables 
                        del_M,d_V, w_E, mu_E,
                        // error
                        //sd_length_n[],
                        //sd_weight_n[],
                        sd_length,
                        sd_weight                  
                        // crowding
                        //W0,W1
                        );
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
          Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction//,
          // crowding
          //W0,W1
          );
    for (j in 1:n){
      obs_hat[i,1,j] =  pow(y_hat[i,j,2],(1.0/3.0))/del_M;
      obs_hat[i,2,j] =  y_hat[i,j,2] + w_E/(mu_E*d_V)*y_hat[i,j,1];
    }
    obs0_hat[i,1] =  pow(y0[i,2],(1.0/3.0))/del_M;
    obs0_hat[i,2] =  y0[i,2] + y0[i,1]*w_E/(mu_E*d_V);
  }
}


