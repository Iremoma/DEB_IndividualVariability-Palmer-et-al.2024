#---------------------------------------
# Assessing between-individual variability in bioenergetics modelling: 
#  opportunities, challenges, and potential applications.
#
# Miquel Palmer, Irene Moro-Martinez, Joaquim Tomas-Ferrer, Amalia Grau, Maria Dolores Lopez-Belluga,
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
    
    #-------
    # Dynamic equations
    #-------
    
    dEdt=pA-pC	                                  # Reserve energy dynamics			             
    dVdt=pG/E_G                                   # Structural volume dynamics
    
    #-------
    # Observable variables
    #-------
    
    # 1) Physical length (cm):
    
    total_length=V^(1/3)/del_M
    
    # 2) Wet weight (gr):
    
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
