DEBmodel=function(time, state, parms) { # t (=time, in Miquel's script; state=y, in Miquel's script) 
  with(as.list(c(state, parms)), {	# params and state are both introduced as a vector.
    
    #-------
    # Temperature and biomass correction
    #-------
    Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(time + birth_correction) + Temp_phi) 
    cT = exp(TA/T1 - TA/Temp)
    
    biomass = (3.173425e+03-3.720773e+02*time+2.362952e+00*time^2-1.431395e-03*time^3)
    W = 1*(biomass<=W0)+
        (1+(-W0/(W0-W1))+(biomass/(W0-W1)))*(biomass>W0)*(biomass<W1)+
        0*(biomass>=W1)
    
    
    pAm_T = pAm*cT*W
    v_T = v*cT
    p_M_T = p_M*cT
    k_J_T = k_J*cT
    
    #-------
    # spawining
    #-------
    #t2=time-(floor(time/365.0)*365.0)
    #spawn=max.spawn*(exp(-(t2-330.0)^2.0/(2.0*30.0^2)) + exp(-(t2+365.0-330.0)^2/(2.0*30.0^2.0)))
    spawn=0  # no spawning dynamics
    
    #-------
    # Fluxes
    #-------
    pA=f*pAm_T*V^(2/3)                              # assimilation rate
    pS=p_M_T*V			                                # somatic maintenance rate
    pC=(E/V)*((E_G*v_T*V^(2/3)+pS)/(kap*E/V+E_G))   # mobilization rate
    pG = kap*pC - pS                                # growth rate
    pJ = k_J*H			                                # maturity maintenance rate
    pR = (1-kap)*pC-pJ                              # reproduction rate
    pD = pS + pJ + (1 - kap_R) * pR                 # dissipation
    
    #-------
    # Dynamic equations
    #-------
    dEdt=pA-pC	                                  # reserve energy dynamics			             
    dVdt=pG/E_G                                   # structural volume dynamics
    dHdt=pR*(H<Hp)                                # maturity dynamics
    dRdt=(kap_R*pR-spawn*R)*(H>=Hp) #-UR*(spaw_event>364)*(UR>0)# reproduction dynamics
    
    #-------
    # Observable variables
    #-------
    
    # 1) Physical length (cm)
    total_length=V^(1/3)/del_M
    
    # 2) Wet weight (gr)
    wet_weight = V + E*w_E/mu_E/d_V  # no spawning dynamics
    #wet_weight = V + (E+R)*w_E/mu_E/d_V  # spawning dynamics

    # 3) food consumed
    p_ADG = matrix(c(pA, pD, pG), nrow = 1, ncol = 3, byrow = T); # fluxes affecting O2 (assimilation, dipispation and growth 3-element vector)
    J_O = eta_O %*% t(p_ADG);                                     # mol/d: fluxes of organic compounds:   J_X, J_V, J_E+JER, J_P in rows (4-element vector)
    J_M = - n_M_n_O %*% J_O;                                      # mol/d: fluxes of "mineral" compounds: J_C, J_H, J_O, J_N in rows (4-element vector)
    w_X=23.9  # gr/mol
    water = 9 # water content(%); from the food producer
    food = - w_X * J_O[1] * (1+water/100);  # g/d, ingested food wet mass

    # 4) feaces
    w_P=23.9 # gr/mol
    feaces =  w_P * J_O[4]; # g/d, feaces dry weight
    
    #-------
    # Output
    #-------
    list(c(dEdt,dVdt,dHdt,dRdt),#,dE_Hdt),
         total_length=total_length,
         wet_weight=wet_weight,
         food=food,
         feaces=feaces
         )
  })
}
