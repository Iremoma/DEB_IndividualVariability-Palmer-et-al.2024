DEBmodel=function(time, state, parms) {
  with(as.list(c(state, parms)), {
    
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
    t2=time-(floor(time/365.0)*365.0)
    spawn=max.spawn*(exp(-(t2-peak.spawn)^2.0/(2.0*spread.spawn^2)) + exp(-(t2+365.0-peak.spawn)^2/(2.0*spread.spawn^2.0)))
    #spawn=0  # no spawning dynamics
    gam_rel_weight=spawn*R*w_E/mu_E/d_V
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
    dRdt=(kap_R*pR-spawn*R)*(H>=Hp)               # reproduction dynamics
    
    #-------
    # Observable variables
    #-------
    
    # 1) Physical length (cm)
    total_length=V^(1/3)/del_M
    
    # 2) Wet weight (gr)
    #wet_weight = V + E*w_E/mu_E/d_V  # no spawning dynamics
    wet_weight = V + (E+R)*w_E/mu_E/d_V  # spawning dynamics
    gonad_weight = R*w_E/mu_E/d_V 
    
    #-------
    # Output
    #-------
    list(c(dEdt,dVdt,dHdt,dRdt),
         total_length=total_length,
         wet_weight=wet_weight,
         gonad_weight=gonad_weight,
         gam_rel_weight=gam_rel_weight
         )
  })
}
