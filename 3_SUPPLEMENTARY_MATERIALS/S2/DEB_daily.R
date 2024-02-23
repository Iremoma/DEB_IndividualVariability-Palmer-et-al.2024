DEBmodel=function(time, state, parms) { # t (=time, in Miquel's script; state=y, in Miquel's script) 
  with(as.list(c(state, parms)), {	# params and state are both introduced as a vector.
    
    #-------
    # Temperature correction
    #-------
    #Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(time + birth_correction) + Temp_phi) 
    Temp = 273.15 + forcing.var(time)   # daily temperatures
    
    cT = exp(TA/T1 - TA/Temp)
    pAm_T = pAm*cT
    v_T = v*cT
    p_M_T = p_M*cT
    
    #-------
    # Fluxes
    #-------
    pA=f*pAm_T*V^(2/3)                              # assimilation rate
    pS=p_M_T*V			                                # somatic maintenance rate
    pC=(E/V)*((E_G*v_T*V^(2/3)+pS)/(kap*E/V+E_G))   # mobilization rate
    pG = kap*pC - pS                                # growth rate
    
    #-------
    # Dynamic equations
    #-------
    dEdt=pA-pC	                                  # reserve energy dynamics			             
    dVdt=pG/E_G                                   # structural volume dynamics
    
    #-------
    # Observable variables
    #-------
    
    # 1) Physical length (cm)
    total_length=V^(1/3)/del_M
    
    # 2) Wet weight (gr)
    wet_weight = V + E*w_E/mu_E/d_V

    #-------
    # Output
    #-------
    list(c(dEdt,dVdt),#,dE_Hdt),
         total_length=total_length,
         wet_weight=wet_weight
         )
  })
}
