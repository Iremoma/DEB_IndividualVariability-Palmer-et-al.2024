#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 8
#---------------------------------------
# Last update: 30 January 2024


# -----
# 1: Loading data and libraries
# -----
remove(list=ls())
load("input.RData") # observed data and standard DEB parameters  
library(deSolve)   
library(bayesplot) 
library(posterior) 
library(ggplot2)   
fit = readRDS("out.RDS") # results of the data analysis as in the main text
draws=fit$draws(format="matrix") 

source("DEB.R") # DEB version for simulation new data
# (1) Temperature correction; (2) without crowding;
# (3) Observable variables modeled:
# Length, weight, O2 consumption, food consumed and feaces

# Additional DEB parameters needed for data simulation
parms$d_V
parms$w_V=23.9
parms$mu_E
parms$kap_X=0.68
parms$mu_X=478000  
parms$kap_P
parms$mu_P=480000
parms$w_P=23.9
parms$w_X=23.9
parms$d_X=0.9

#-----
# 2: data simulation
#-----
# parameters from the observed data
pAm = median(data.frame(as_draws_df(fit$draws("pAm_mu")))[,1])
v = median(data.frame(as_draws_df(fit$draws("v_mu")))[,1])
E0 = median(data.frame(as_draws_df(fit$draws("E0_mu")))[,1])
V0 =  median(data.frame(as_draws_df(fit$draws("V0_mu")))[,1])

# simulation settings
t0=236 # initial age 
tf=300 # final age
n=30   # number of observations
times=round(seq(t0,tf,length.out=n+1))

E_H0 = 27382.78 # Supp mat S5
state=NULL
state["E"] = E0
state["V"] = V0
state["E_H"] = E_H0
parms2=parms
parms2$del_M = 0.2805034 # Supp mat S6
parms2$pAm = pAm
parms2$v = v
parms2$Hp=35889.07 # Supp mat S5
parms2$birth_correction=as.numeric(birth_correction)

res = ode(y=state,
          times=times,
          func=DEBmodel,
          parms=parms2,
          method="ode45")
res=data.frame(res)


age=times[2:(n+1)]

# setting observatio errors (to be estimated)
sd_length=0.1
sd_weight=1
sd_energy=10
sd_ox=0.005
sd_food=0.005
sd_feaces=0.001

# templates for the simulated data
obs=array(NA,dim=c(n,6))
obs[,1]=rnorm(n,res$total_length[2:(n+1)],sd_length)
obs[,2]=rnorm(n,res$wet_weight[2:(n+1)],sd_weight)
obs[,3]=rnorm(n,res$E[2:(n+1)]/res$V[2:(n+1)],sd_energy)
obs[,4]=rnorm(n,res$ox_consumption[2:(n+1)],sd_ox)
obs[,5]=rnorm(n,res$food[2:(n+1)],sd_food)
obs[,6]=rnorm(n,res$feaces[2:(n+1)],sd_feaces)
obs0=rep(NA,6)
obs0[1]=rnorm(1,res$total_length[1],sd_length)
obs0[2]=rnorm(1,res$wet_weight[1],sd_weight)
obs0[3]=rnorm(1,res$E[1]/res$V[1],sd_energy)
obs0[4]=rnorm(1,res$ox_consumption[1],sd_ox)
obs0[5]=rnorm(1,res$food[1],sd_food)
obs0[6]=rnorm(1,res$feaces[1],sd_feaces)


# save(t0,n,age,obs,obs0,res,state,
#       sd_length,sd_weight,sd_energy,sd_ox,sd_food,sd_feaces,
#      state,parms2,file="input_S8.RData")


#-----
# 3: Analysis
#-----
#---- 
# 3.1: Re-loading input (simulated data) and libraries
#----
remove(list=ls())
library(cmdstanr)
load("input_S8.RData")

#---- 
# 3.2: STAN model
#----
sink("DEB.stan")
cat(" // first line
functions {
  // two functions: DEB function and observable rates function
  //////////////////////////////
  // DEB
  vector DEB(
    real t,          // time
    vector y,        // state variables
    // core parameters
    real f,
    real E_G,
    real kap,    
    real pAm, 
    real p_M,
    real k_J,
    real v,
    real kap_R,
    real Hp,
    
    // temperature parameters
    real Kelvin,
    real Temp_mean,
    real Temp_amp,
    real pi2f,
    real Temp_phi,
    real TA,
    real T1,
    real birth_correction
    ){
      // Derivatives
      vector[3] dydt;
      // 1: dEdt;
      // 2: dVdt;
      // 3: dE_Hdt;

      // Auxiliary variables
      // fluxes (j/day)
      real pA; 
      real pS;
      real pC;
      real pG;
      real pJ;
      real pR;
      real pD;
      
      // auxiliary variables and parameters
      real Temp;         // Temperature
      real cT;           // Arrhenius temperature correction

      real pAm_T;        // temperature corrected parameters
      real v_T;
      real p_M_T;
      real k_J_T;
      
      // temperature correction
      Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(t + birth_correction) + Temp_phi);
      cT = exp(TA/T1 - TA/Temp);
      //cT=1.0;
      
      // Temperature-corrected parameters      
      pAm_T = cT*pAm;             
      v_T = cT*v;
      p_M_T = cT*p_M;
      k_J_T = cT*k_J;
      
      // fluxes (j/day)
      pA = f*pAm_T*y[2]^(2.0/3.0);                                                // assimilation rate    
      pS = p_M_T*y[2];			                                                      // somatic maintenance rate;			                                           // maturity maintenance rate
      pC = (y[1]/y[2])*((v_T*E_G*pow(y[2],2.0/3.0)+pS)/(E_G + kap*y[1]/y[2]));    // mobilization rate
      pG = (kap*pC - pS);                                                         // growth rate (priority for maintenance)
      pJ = k_J_T*y[3];                                                            // maturity maintenance
      pR = (1-kap)*pC-pJ;                                                         // reproduction
      pD = pS + pJ + (1 - kap_R) * pR;                                            // disipation
      
      // derivatives
      dydt[1]=pA-pC;
      dydt[2]=pG/E_G;
      dydt[3] = pR; //*(y[3]<Hp);
      
      return dydt;
    }
  //////////////////////////////
  // observable rates
  array [,] real observable_rates(
    // number of samples
    int n,
    // state variables
    array [] vector y,  //E,V and E_H      
    // core parameters
    real f,
    real E_G,
    real kap,    
    real pAm, 
    real p_M,
    real k_J,
    real v,
    real kap_R,
    real kap_X,
    real Hp,
    
    // temperature parameters
    real Kelvin,
    real Temp_mean,
    real Temp_amp,
    real pi2f,
    real Temp_phi,
    real TA,
    real T1,
    real birth_correction,
    array [] real age,
    
    // auxiliary paramters
    //array [,] real eta_O,
    array [,] real n_M_n_O,
    real d_V,
    real w_V,
    real mu_E,
    real mu_X,
    real kap_P,
    real mu_P,
    real w_P,
    real w_X,
    real d_X

    ){
      array [3,n] real rates; // expected rates (O2/day,food/day,feaces/day)

      // Auxiliary variables
      // fluxes (j/day)
      array[n] real pA; 
      array[n] real pS;
      array[n] real pC;
      array[n] real pG;
      array[n] real pJ;
      array[n] real pR;
      array[n] real pD;
      
      // Arrhenius temperature correction
      array[n] real Temp;           // Temperature
      array[n] real cT;             // correction
      // temperature corrected parameters
      array[n] real pAm_T;        
      array[n] real v_T;
      array[n] real p_M_T;
      array[n] real k_J_T;
      // yields
      real M_V;
      real y_V_E;
      real y_E_X;
      real y_X_E;
      real y_P_X;
      real y_P_E;
      // Mass-power couplers
      real eta_XA;
      real eta_PA;
      real eta_VG;
      array [4,3] real eta_O;
      // 
      array [3,n] real p_ADG;
      array [4,n] real J_O;
      array [4,n] real J_M;
      
      for (i in 1:n){
        // temperature correction
        Temp[i] = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(age[i] + birth_correction) + Temp_phi);
        cT[i] = exp(TA/T1 - TA/Temp[i]);
        //cT[i]=1.0;
        
        // Temperature-corrected parameters      
        pAm_T[i] = cT[i]*pAm;             
        v_T[i] = cT[i]*v;
        p_M_T[i] = cT[i]*p_M;
        k_J_T[i] = cT[i]*k_J;
        
        // fluxes
        pA[i] = f*pAm_T[i]*y[i,2]^(2.0/3.0);                                        
        pS[i] = p_M_T[i]*y[i,2];			                                                     
        pC[i] = (y[i,1]/y[i,2])*((v_T[i]*E_G*y[i,2]^(2.0/3.0)+pS[i])/(E_G + kap*y[i,1]/y[i,2]));
        pG[i] = (kap*pC[i] - pS[i]);                                                  
        pJ[i] = k_J_T[i]*y[i,3];
        pR[i] = (1-kap)*pC[i]-pJ[i];
        pD[i] = pS[i] + pJ[i] + (1 - kap_R) * pR[i];
        
        // Oxygen 
        p_ADG[1,i] = pA[i];
        p_ADG[2,i] = pD[i];
        p_ADG[3,i] = pG[i];
      }
      
      // eta_O, J_O and J_M
      M_V = d_V/ w_V;
      y_V_E = mu_E * M_V/ E_G; 
      y_E_X  = kap_X * mu_X/ mu_E;
      y_X_E  = 1/ y_E_X;
      y_P_X  = kap_P *mu_X/ mu_P;
      y_P_E  = y_P_X/ y_E_X;   
      eta_XA = y_X_E/ mu_E;
      eta_PA = y_P_E/ mu_E;
      eta_VG = y_V_E/ mu_E;
      eta_O[1,1]=-eta_XA;
      eta_O[1,2]=0;
      eta_O[1,3]=0;
      eta_O[2,1]=0;
      eta_O[2,2]=0;
      eta_O[2,3]=eta_VG;
      eta_O[3,1]=1.0/mu_E;
      eta_O[3,2]=-1.0/mu_E;
      eta_O[3,3]=-1.0/mu_E;
      eta_O[4,1]=eta_PA;
      eta_O[4,2]=0;
      eta_O[4,3]=0;
      J_O = to_array_2d(to_matrix(eta_O[1:4,1:3]) * to_matrix(p_ADG[1:3,1:n])); 
      J_M = to_array_2d(- to_matrix(n_M_n_O[1:4,1:4]) * to_matrix(J_O[1:4,1:n]));
           
      for (i in 1:n){  
        rates[1,i] = -32.0 * J_M[3,i]; // (gr/day) O2 consumption per day
        rates[2,i] = - w_X * J_O[1,i] / d_X;  // (g/day), ingested food wet mass
        rates[3,i] = w_P * J_O[4,i];  // (g/day), feaces
      }
      return rates;
    }
  
}

data {
  int n;                           // number of observations per fish
  array [n,6] real obs;            // observations; repeated measures (1:n),c(length,weight,E/V, O2,food,feaces)
  real t0;                         // initial time (1 days before the firts sampling date)
  array [n] real age;              // times at which the system is observed
  array [6] real obs0;             // observations at t0; c(length,weight,E/V, O2,food,feaces)

  // parameters for the integrate function (ODE) 
  real f;           
  //real E_G;         
  //real kap;     
  //real p_M;
  real k_J;
  real kap_R;
  real Hp;
    
  // other parameters
  real del_M;        
  real w_E;        
  //real mu_E;        
  //real d_V;
  
  // temperature parameters
  real Kelvin;
  real Temp_mean;
  real Temp_amp;
  real pi2f;
  real Temp_phi;
  real TA;
  real T1;
  real birth_correction;
  
  // rates parameters
  //array [4,3] real eta_O;
  array [4,4] real n_M_n_O;
  real d_V;
  real w_V;
  real mu_E;
  real mu_X;
  real kap_P;
  real mu_P;
  real w_P;
  real w_X;
  real d_X;
    
  // priors
  array [2] real prior_V0;
  array [2] real prior_E0;
  array [2] real prior_E_H0;
  array [2] real prior_pAm;        
  array [2] real prior_v;
  array [2] real prior_p_M;
  array [2] real prior_E_G;         
  array [2] real prior_kap;
  array [2] real prior_kap_X;
  array [2] real prior_sd_length;
  array [2] real prior_sd_weight;
  array [2] real prior_sd_energy;
  array [2] real prior_sd_ox;
  array [2] real prior_sd_food;
  array [2] real prior_sd_feaces;
}

transformed data {
}

parameters {
  // DEB parameters
  real <lower=0> pAm;
  real <lower=0> v;
  real <lower=0> p_M;
  real <lower=0> E_G;         
  real <lower=0, upper=1> kap;
  real <lower=0, upper=1> kap_X;

  // initial state
  vector <lower=0> [3] y0; 
  
  // observation error
  real <lower=0> sd_length;
  real <lower=0> sd_weight;
  real <lower=0> sd_energy;
  real <lower=0> sd_ox;
  real <lower=0> sd_food;
  real <lower=0> sd_feaces;
}

transformed parameters {
}   

model {
  array [n] vector[3] y;
  array [3,n] real rates;  //oxygen and food

  pAm ~ normal(prior_pAm[1] , prior_pAm[2]);
  v ~ normal(prior_v[1] , prior_v[2]);
  p_M ~ normal(prior_p_M[1] , prior_p_M[2]);
  E_G ~ normal(prior_E_G[1] , prior_E_G[2]);
  kap ~ beta(prior_kap[1] , prior_kap[2]);
  kap_X ~ beta(prior_kap_X[1] , prior_kap_X[2]);

  y0[1] ~ normal(prior_E0[1], prior_E0[2]); // reserve energy
  y0[2] ~ normal(prior_V0[1], prior_V0[2]); // structural volume
  y0[3] ~ normal(prior_E_H0[1], prior_E_H0[2]); // maturation energy
  obs0[1] ~ normal( ((y0[2]^(1.0/3.0))/del_M), sd_length);        // length at t0
  obs0[2] ~ normal( (y0[2] + y0[1]*w_E/(mu_E*d_V)), sd_weight);   // weight at t0
  // no likelihood for rates at t0

  sd_length ~ normal(prior_sd_length[1],prior_sd_length[2]);
  sd_weight ~ normal(prior_sd_weight[1],prior_sd_weight[2]);
  sd_energy ~ normal(prior_sd_energy[1],prior_sd_energy[2]);
  sd_ox ~ normal(prior_sd_ox[1],prior_sd_ox[2]);
  sd_food ~ normal(prior_sd_food[1],prior_sd_food[2]);
  sd_feaces ~ normal(prior_sd_feaces[1],prior_sd_feaces[2]);
    
  y[1:n] = ode_rk45(DEB, y0[], t0 , age[],
            // DEB main parameters
            f,E_G,kap,pAm,p_M,k_J,v,kap_R,Hp,
            // temperature-related paramters
            Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction
        );
  rates[1:3,1:n] = observable_rates(n,y,
                 // DEB main parameters
                 f,E_G,kap,pAm,p_M,k_J,v,kap_R,kap_X,Hp,
                 // temperature-related paramters
                 Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction,
                 age,
                 // rates
                 // eta_O[,],
                 n_M_n_O[,],d_V,w_V,mu_E,mu_X,kap_P,mu_P,w_P,w_X,d_X
                 );
  
  // likelihood
  for (i in 1:n){
      obs[i,1] ~ normal(pow(y[i,2],(1.0/3.0))/del_M,sd_length);
      obs[i,2] ~ normal(y[i,2] + y[i,1]*w_E/mu_E/d_V,sd_weight);
      obs[i,3] ~ normal(y[i,1]/y[i,2],sd_energy);
      obs[i,4] ~ normal(rates[1,i],sd_ox);
      obs[i,5] ~ normal(rates[2,i],sd_food);
      obs[i,6] ~ normal(rates[3,i],sd_feaces);
  }
}

generated quantities {
  array [n] vector[3] y_hat;
  array [3,n] real rates_hat;
  array [n,6] real obs_hat;

  y_hat[1:n] = ode_rk45(DEB, y0[], t0 , age[],
            // DEB main parameters
            f,E_G,kap,pAm,p_M,k_J,v,kap_R,Hp,
            // temperature-related paramters
            Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction
        );
  rates_hat[1:3,1:n] = observable_rates(n,y_hat,
                 // DEB main parameters
                 f,E_G,kap,pAm,p_M,k_J,v,kap_R,kap_X,Hp,
                 // temperature-related paramters
                 Kelvin,Temp_mean,Temp_amp,pi2f,Temp_phi,TA,T1,birth_correction,
                 age,
                 // rates
                 //eta_O[,],
                 n_M_n_O[,],d_V,w_V,mu_E,mu_X,kap_P,mu_P,w_P,w_X,d_X
                 );
  for (i in 1:n){
    obs_hat[i,1] =  pow(y_hat[i,2],(1.0/3.0))/del_M;
    obs_hat[i,2] =  y_hat[i,2] + y_hat[i,1]*w_E/mu_E/d_V;
    obs_hat[i,3] =  y_hat[i,1]/y_hat[i,2];
    obs_hat[i,4] =  rates_hat[1,i];
    obs_hat[i,5] =  rates_hat[2,i];
    obs_hat[i,6] =  rates_hat[3,i];
  }
}

" # end of model code
,fill = TRUE)
sink()

#-----
# 3.3: compiling model
#-----
mod = cmdstan_model("DEB.stan")
#mod = cmdstan_model("DEB.stan",cpp_options = list(stan_threads = TRUE))

#-----
# 3.4: initializating model
#-----
n.chains = 3
initializer = function() list(
  "pAm"=parms2$pAm,
  "v"=parms2$v,
  "p_M"=parms2$p_M,
  "E_G"=parms2$E_G,
  "kap"=parms2$kap,
  "kap_X"=parms2$kap_X,
  #"del_M"=parms$del_M,
  
  "y0"=c((obs0[2]-(obs0[1]*parms2$del_M)^3)/(parms2$w_E/parms2$mu_E/parms2$d_V),
         (obs0[1]*parms2$del_M)^3,
         state["E_H"]),
  
  "sd_length"=sd_length,
  "sd_weight"=sd_weight,
  "sd_energy"=sd_energy,
  "sd_ox"=sd_ox,
  "sd_food"=sd_food,
  "sd_feaces"=sd_feaces
  
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()


#------
# 3.5: priors
#------
prior_V0=c((obs0[1]*parms2$del_M)^3,
           0.5*(obs0[1]*parms2$del_M)^3
)
prior_E0=c((obs0[2]-(obs0[1]*parms2$del_M)^3)/(parms2$w_E/parms2$mu_E/parms2$d_V),
           0.5*(obs0[2]-(obs0[1]*parms2$del_M)^3)/(parms2$w_E/parms2$mu_E/parms2$d_V)
)
prior_E_H0=c(state["E_H"],0.5*state["E_H"])
prior_pAm=c(parms2$pAm,0.5*parms2$pAm)
prior_v=c(parms2$v,0.5*parms2$v)
prior_p_M=c(parms2$p_M,0.5*parms2$p_M)
prior_E_G=c(parms2$E_G,0.5*parms2$E_G)
# f.beta=function(p){
#   temp=rbeta(100000,p[1],p[2])
#   temp=temp[which(temp>0 & temp<1)]
#   q=quantile(temp,c(0.025,0.975))
#   #(q[1]-0.8)^2+(q[2]-0.99)^2
#   (q[1]-0.5)^2+(q[2]-0.8)^2
# }
# #temp=optim(c(10,1),f.beta)
# temp=optim(c(25,12),f.beta)
# quantile(rbeta(10000,temp$par[1],temp$par[2]),c(0.025,0.975))
# temp$par
# #[1]  5.952732 1.255371   # kap   [0.5 : 0.99]
# #[1]  2.099596 2.098306    # kap_X  [0.1 : 0.9]
prior_kap=c( 5.952732, 1.255371)
prior_kap_X=c( 2.099596, 2.098306)
prior_sd_length=c(sd_length,2*sd_length)
prior_sd_weight=c(sd_weight,2*sd_weight)
prior_sd_energy=c(sd_energy,2*sd_energy)
prior_sd_ox=c(sd_ox,2*sd_ox)
prior_sd_food=c(sd_food,2*sd_food)
prior_sd_feaces=c(sd_feaces,2*sd_feaces)

#-----
# 3.6: running
#-----

fit = mod$sample(
  data =list (
    n = n,                                # number of observations per fish (number sampling events)
    obs = obs,                            # observations
    obs0 = obs0,
    t0 = t0,                              # initial time
    age = age,                             # times at which the system is observed (sampling dates)
    
    # DEB main paramteters
    f = parms2$f,                          # functional response
    k_J = parms2$k_J,
    kap_R = parms2$kap_R,
    
    # Temperature
    Kelvin=parms2$Kelvin,
    Temp_mean=parms2$Temp_mean,
    Temp_amp=parms2$Temp_amp,
    pi2f=parms2$pi2f,
    Temp_phi=parms2$Temp_phi,
    TA=parms2$TA,
    T1=parms2$T1,
    birth_correction=parms2$birth_correction,
    
    # DEB secondary parameters
    del_M = parms2$del_M,                      # shape coefficient
    w_E = parms2$w_E,       
    n_M_n_O = parms2$n_M_n_O,
    d_V=parms2$d_V,
    w_V=parms2$w_V,
    mu_E=parms2$mu_E,
    mu_X=parms2$mu_X,
    kap_P=parms2$kap_P,
    mu_P=parms2$mu_P,
    w_P=parms2$w_P,
    w_X=parms2$w_X,
    d_X=parms2$d_X,
    Hp=parms2$Hp,
    
    #priors
    prior_V0=prior_V0,
    prior_E0=prior_E0,
    prior_E_H0=prior_E_H0,
    prior_pAm=prior_pAm,       
    prior_v=prior_v,
    prior_p_M=prior_p_M,
    prior_E_G=prior_E_G,         
    prior_kap=prior_kap,
    prior_kap_X=prior_kap_X,
    prior_sd_length=prior_sd_length,
    prior_sd_weight=prior_sd_weight,
    prior_sd_energy=prior_sd_energy,
    prior_sd_ox=prior_sd_ox,
    prior_sd_food=prior_sd_food,
    prior_sd_feaces=prior_sd_feaces
    
  ),
  chains = n.chains,
  parallel_chains = n.chains,
  #threads_per_chain = 8,
  iter_warmup = 2000,
  iter_sampling = 5000,
  init = inits,
  max_treedepth = 15,
  adapt_delta = 0.99
)

#-----
# 3.7: saving
#-----
# fit$save_object(file = "S8_rich.RDS")
# sink(file = "dignosi_S8_rich.txt")
# fit$cmdstan_diagnose()
# sink()
# file.rename("DEB.stan","model_S8_rich.R")

#-----
# 4: Exploring results
#-----

#-----
# 4.1: Re-loading data and libraries
#-----
remove(list=ls())
library(cmdstanr)
library(bayesplot) 
library(posterior) 
library(ggplot2)
load("input_S8.RData")
fit = readRDS("S8_rich.RDS")
draws=fit$draws(format="matrix")

#-----
# 4.2: Convergence quality control (parameter specific checking)
#-----
mcmc_trace(fit$draws(c("pAm")))

#-----
# 4.3: Tables
#-----
table=fit$summary(c("pAm","v","p_M","E_G","kap","kap_X","y0"))
table=data.frame(table)[,c(1,6,3,7,8,9)]
table[,c(2,3,4,5)]=round(table[,c(2,3,4,5)],3)
table[,6]=round(table[,6],0)
table=cbind(table,round(c(parms2$pAm,
                          parms2$v,
                          parms2$p_M,
                          parms2$E_G,
                          parms2$kap,
                          parms2$kap_X,
                          state["E"],
                          state["V"],
                          state["E_H"]),3)
)
names(table)[7]="True value"
table

#-----
# 4.4: Re-setting priors (as defined for the analysis)
#-----
iter=length(data.frame(as_draws_df(fit$draws("pAm")))[,1])
prior_V0=c((obs0[1]*parms2$del_M)^3,
           0.5*(obs0[1]*parms2$del_M)^3
)
prior_E0=c((obs0[2]-(obs0[1]*parms2$del_M)^3)/(parms2$w_E/parms2$mu_E/parms2$d_V),
           0.5*(obs0[2]-(obs0[1]*parms2$del_M)^3)/(parms2$w_E/parms2$mu_E/parms2$d_V)
)
prior_E_H0=c(state["E_H"],0.5*state["E_H"])
prior_pAm=c(parms2$pAm,0.5*parms2$pAm)
prior_v=c(parms2$v,0.5*parms2$v)
prior_p_M=c(parms2$p_M,0.5*parms2$p_M)
prior_E_G=c(parms2$E_G,0.5*parms2$E_G)
prior_kap=c( 5.952732, 1.255371)
prior_kap_X=c( 2.099596, 2.098306)
prior_sd_length=c(sd_length,2*sd_length)
prior_sd_weight=c(sd_weight,2*sd_weight)
prior_sd_energy=c(sd_energy,2*sd_energy)
prior_sd_ox=c(sd_ox,2*sd_ox)
prior_sd_food=c(sd_food,2*sd_food)
prior_sd_feaces=c(sd_feaces,2*sd_feaces)



#-----
# Figure S9: Prior (blue), posterior distribution (red) and true values 
# (vertical dashed line) for the DEB parameters estimated in a data rich scenario 
# (length, weight, energy density, oxygen consumption, food consumed, and fecal 
# matter produced). Note that in this case, precision of the estimates (posterior 
# distribution) is clearly narrower than the prior distribution in all the cases, 
# excepting for E_H0 (maturity energy at t0).
#----- 


par(mfrow=c(3,3))

# f*pAm:

temp=rnorm(iter,prior_pAm[1],prior_pAm[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("pAm")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(bold(symbol("\042"))),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=parms2$pAm,lty=2)

# v:

temp=rnorm(iter,prior_v[1],prior_v[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("v")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(dot(nu)),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=parms2$v,lty=2)

#p_M: 

temp=rnorm(iter,prior_p_M[1],prior_p_M[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("p_M")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(dot(p)*S),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=parms2$p_M,lty=2)

#E_G:

temp=rnorm(iter,prior_E_G[1],prior_E_G[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("E_G")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main = bquote(italic('[' ~ EG ~ ']'))
,xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=parms2$E_G,lty=2)

# kap:

temp=rbeta(iter,prior_kap[1],prior_kap[2])
temp=temp[which(temp>0 & temp<1)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("kap")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(kappa),xlab="",
     xlim=c(0,1),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=parms2$kap,lty=2)

#kap_X:

temp=rbeta(iter,prior_kap_X[1],prior_kap_X[2])
temp=temp[which(temp>0 & temp<1)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("kap_X")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(kappa[~X]),xlab="",
     xlim=c(0,1),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=parms2$kap_X,lty=2)

#E0:

temp=rnorm(iter,prior_E0[1],prior_E0[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("y0")))[,1]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(E[~0]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=state["E"],lty=2)

#V0:

temp=rnorm(iter,prior_V0[1],prior_V0[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("y0")))[,2]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(V[~0]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=state["V"],lty=2)

#E_H0:

temp=rnorm(iter,prior_E_H0[1],prior_E_H0[2])
temp=temp[which(temp>0)]
density_prior=density(temp)
temp=data.frame(as_draws_df(fit$draws("y0")))[,3]
density_post=density(temp)
plot(density_prior$x,density_prior$y,type="l",col="blue",ylab="Frequency",main=bquote(E[~H0]),xlab="",
     xlim=c(0,max(density_prior$x)),
     ylim=range(c(density_post$y,density_prior$y)))
lines(density_post$x,density_post$y,col="red")
abline(v=state["E_H"],lty=2)





