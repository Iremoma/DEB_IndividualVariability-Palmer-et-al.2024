DEBmodel=function(time, state, parms) { 
  with(as.list(c(state, parms)), {
    
    Temp = Kelvin + Temp_mean + Temp_amp*sin(pi2f*(time + birth_correction) + Temp_phi) 
    cT = exp(TA/T1 - TA/Temp)
    pAm_T = pAm*cT
    v_T = v*cT
    p_M_T = p_M*cT
    k_J_T = cT*k_J
    
    #-------
    # Fluxes
    #-------
    
    pA = f*pAm_T*V^(2/3)                              # assimilation rate
    pS = p_M_T*V			                                # somatic maintenance rate
    pC = (E/V)*((E_G*v_T*V^(2/3)+pS)/(kap*E/V+E_G))   # mobilization rate
    pG = kap*pC - pS                                  # growth rate
    pJ = k_J_T*E_H			                              # maturity maintenance rate
    pR = (1-kap)*pC-pJ                                # reproduction rate; priority for maturity maintenance
    pD = pS + pJ + (1 - kap_R) * pR                   # CAUTION! dissipation
    
    #-------
    # Dynamic equations
    #-------
    dEdt=pA-pC	                                  # reserve energy dynamics			             
    dVdt=pG/E_G                                   # structural volume dynamics
    dE_Hdt = pR#*(E_H<Hp)                          # maturity dynamics
    #dRdt = kappaR*pR*(E_H>=Hp)                   # reproduction dynamics
    
    #-------
    # Observable variables
    #-------
    # eta_O
    
    # Yields
    M_V     = d_V/ w_V;                   # mol/cm^3, [M_V], volume-specific mass of structure
    y_V_E   = mu_E * M_V/ E_G;            # mol/mol, yield of structure on reserve
    y_E_X  = kap_X * mu_X/ mu_E;          # mol/mol, yield of reserve on food
    y_X_E  = 1/ y_E_X;                    # mol/mol, yield of food on reserve
    y_P_X  = kap_P *mu_X/ mu_P;           # mol/mol, yield of faeces on food 
    y_P_E  = y_P_X/ y_E_X;                # mol/mol, yield of faeces on reserve            # mol/mol, yield of faeces on reserve
    
    #  Mass-power couplers
    eta_XA = y_X_E/ mu_E;          # mol/J, food-assim energy coupler
    eta_PA = y_P_E/ mu_E;          # mol/J, faeces-assim energy coupler
    eta_VG = y_V_E/ mu_E;          # mol/J, struct-growth energy coupler
    eta_O  = matrix(c(-eta_XA,     0,         0,         # mol/J, mass-energy coupler
                      0,        0,       eta_VG,       # used in: J_O = eta_O * p
                      1/mu_E, -1/mu_E,  -1/mu_E,
                      eta_PA,    0,         0),
                      nrow = 4, ncol = 3, byrow = T); 
    
    # 1) Physical length (cm)
    total_length=V^(1/3)/del_M
    
    # 2) Wet weight (gr)
    wet_weight = V + E*w_E/mu_E/d_V

    # 3) O2 consumption 
    p_ADG = matrix(c(pA, pD, pG), nrow = 1, ncol = 3, byrow = T); # fluxes affecting O2 (assimilation, dipispation and growth 3-element vector)
    J_O = eta_O %*% t(p_ADG);                                     # mol/d: fluxes of organic compounds:   J_X, J_V, J_E+JER, J_P in rows (4-element vector)
    J_M = - n_M_n_O %*% J_O;                                      # mol/d: fluxes of "mineral" compounds: J_C, J_H, J_O, J_N in rows (4-element vector)
    ox_consumption =  - 2* 16 * J_M[3];                           # g/d, O2 consumption 
    
    # 4) food consumed
    #w_X=23.9
    #d_X=0.9
    food = - w_X * J_O[1] / d_X;  # g/d, ingested food wet mass
    
    # 5) feaces
    feaces =  w_P * J_O[4]; # g/d, feaces dry weight
    
    #-------
    # Output
    #-------
    list(c(dEdt,dVdt,dE_Hdt),total_length=total_length,wet_weight=wet_weight,ox_consumption=ox_consumption,food=food,feaces=feaces)
  })
}
