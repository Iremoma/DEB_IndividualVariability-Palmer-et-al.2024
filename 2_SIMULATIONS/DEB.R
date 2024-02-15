#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer,Irene Moro-Martinez, Joaquim Tomàs-Ferrer, Amalia Grau, María Dolores López-Belluga,
#  Marine Herlin, Orestis Stavrakidis-Zachoua, Andrea Campos-Candela.
#
# More information related to the Script can be find in Supplementary materials 7
#---------------------------------------
# Last update: 25 January 2024

#-----
# DEB model function
#-----


DEBmodel=function(time, state, parms) { 
  with(as.list(c(state, parms)), {	
    
    #-------
    # Temperature correction
    #-------
    
    Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(time + birth_correction) + Temp_phi) 
    cT = exp(TA/T1 - TA/Temp)

    
    W0=27652.82
    W1=1118730
    biomass = (3.173425e+03-3.720773e+02*time+2.362952e+00*time^2-1.431395e-03*time^3)
    W = 1*(biomass<=W0)+
      (1+(-W0/(W0-W1))+(biomass/(W0-W1)))*(biomass>W0)*(biomass<W1)+
      0*(biomass>=W1)
    
    pAm_T = pAm*cT
    v_T = v*cT
    p_M_T = p_M*cT  
    
    #-------
    # Fluxes
    #-------
    
    pA=f*pAm_T*V^(2/3)                            # Assimilation rate
    pS=p_M_T*V			                              # Somatic maintenance rate
    pC=(E/V)*((E_G*v_T*V^(2/3)+pS)/(kap*E/V+E_G)) # Mobilization rate
    pG = kap*pC - pS                              # Growth rate
    #pJ = k_J*E_H			                            # Maturity maintenance rate
    #pR = (1-kap)*pC-pJ                           # Reproduction rate; priority for maturity maintenance
    pD = pS #+ pJ + (1 - kap_R) * pR ;            # Dissipation
    
    #-------
    # Dynamic equations
    #-------
    
    dEdt=pA-pC	                                  # Reserve energy dynamics			             
    dVdt=pG/E_G                                   # Structural volume dynamics
    #dE_Hdt = pR#*(E_H<Hp)                        # Maturity dynamics
    #dRdt = kappaR*pR*(E_H>=Hp)                   # Reproduction dynamics
    
    #-------
    # Observable variables
    #-------
    
    # 1) Physical length (cm):
    
    total_length=V^(1/3)/del_M
    
    # 2) Wet weight (gr):
    
    wet_weight = V + E*w_E/mu_E/d_V

    # 3) Food consumed:
    
    p_ADG = matrix(c(pA, pD, pG), nrow = 1, ncol = 3, byrow = T); # fluxes affecting O2 (assimilation, dipispation and growth 3-element vector)
    J_O = eta_O %*% t(p_ADG);                                     # mol/d: fluxes of organic compounds:   J_X, J_V, J_E+JER, J_P in rows (4-element vector)
    J_M = - n_M_n_O %*% J_O;                                      # mol/d: fluxes of "mineral" compounds: J_C, J_H, J_O, J_N in rows (4-element vector)
    w_X=23.9  # gr/mol
    water=9 # water content(%); from the food producer
    food = - w_X * J_O[1] * (1+water/100);  # g/d, ingested food wet mass

    # 4) Feaces:
    
    w_P=23.9 # gr/mol
    feaces =  w_P * J_O[4]; # g/d, feaces dry weight
    
    #-------
    # Output
    #-------
    
    list(c(dEdt,dVdt),#,dE_Hdt),
         total_length=total_length,
         wet_weight=wet_weight,
         food=food,
         feaces=feaces
         )
  })
}
