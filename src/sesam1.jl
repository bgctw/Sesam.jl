function sesam1(;name)
    @parameters t [unit = uT]
    D = Differential(t)

    @parameters ϵ ϵ_tvr κ_E a_E [unit = uT^-1] m [unit = uT^-1] τ [unit = uT^-1] 
    @parameters k_L [unit = uT^-1] k_R [unit = uT^-1] k_mN [unit = uQ*uT^-1] 
    @parameters α_R0

    # state variables
    @variables B(t) [unit = uQ] L(t) [unit = uQ]  R(t) [unit = uQ] α_R(t) cumresp(t) [unit = uQ]

    @variables α_L(t) tmpNoUnit(t)
    # fluxes
    @variables (begin
        syn_Enz(t),[unit = uQ/uT], r_M(t),[unit = uQ/uT], tvr_B(t),[unit = uQ/uT],
        dec_LPot(t),[unit = uQ/uT], dec_L(t),[unit = uQ/uT], dec_RPot(t),[unit = uQ/uT],
        dec_R(t),[unit = uQ/uT], u_C(t),[unit = uQ/uT],
        C_synBCt(t),[unit = uQ/uT], C_synBC(t),[unit = uQ/uT], r_G(t),[unit = uQ/uT],
        r_tvr(t),[unit = uQ/uT], syn_B(t),[unit = uQ/uT],
        r_B(t),[unit = uQ/uT], r_GEnz(t),[unit = uQ/uT], r_O(t),[unit = uQ/uT]
    end)

    #------------- N ---------------------
    @variables L_N(t) [unit = uQ] R_N(t) [unit = uQ] I_N(t) [unit = uQ]
    @variables (begin
        Φ_N(t),[unit = uQ/uT], Φ_Nu(t),[unit = uQ/uT], Φ_NB(t),[unit = uQ/uT], 
        u_PlantN(t),[unit = uQ/uT], u_NOM(t),[unit = uQ/uT], u_immNPot(t),[unit = uQ/uT], 
        u_N(t),[unit = uQ/uT], N_synBN(t),[unit = uQ/uT], M_ImbN(t),[unit = uQ/uT], 
        β_NL(t), β_NR(t) 
    end)
    @parameters β_NE β_NB l_N [unit = uT^-1] ν_N i_BN [unit = uT^-1] 

    # need to be specified by coupled system:
    @variables (begin
        i_L(t),[unit = uQ/uT], β_Ni(t), i_IN(t),[unit = uQ/uT],
        u_PlantNmax(t),[unit = uQ/uT], k_PlantN(t),[unit = uT^-1]
    end) 

    eqs = [
        #tmpNoUnit ~ ustrip_num(B),
        D(B) ~ syn_B - tvr_B,
        D(L) ~ -dec_L + i_L,
        D(R) ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
        D(α_R) ~ 0.0,
        # C cycle
        syn_Enz ~ a_E*B,
        r_M ~ m*B,
        tvr_B ~ τ*B,
        dec_LPot ~ k_L * L,
        dec_L ~ dec_LPot*(α_L * syn_Enz)/(k_mN + α_L*syn_Enz),
        dec_RPot ~ k_R * R,
        dec_R ~ dec_RPot*(α_R * syn_Enz)/(k_mN + α_R*syn_Enz),
        u_C ~ dec_L + dec_R + κ_E*syn_Enz,
        C_synBCt ~ u_C - syn_Enz/ϵ - r_M,
        # TODO make ifelse work with Unitful
        #tmpNoUnit ~ 1.0u"g*m^-2*yr^-1" > 0.0u"g*m^-2*yr^-1",
        #tmpNoUnit ~ C_synBCt.val > 0.0u"g*m^-2*yr^-1",
        #C_synBC ~ IfElse.ifelse(C_synBCt > 0.0u"g*m^-2*yr^-1", ϵ*C_synBCt, C_synBCt), 
        C_synBC ~ ϵ*C_synBCt,
        # N cycle
        β_NL ~ L/L_N,
        β_NR ~ R/R_N,
        D(L_N) ~ -dec_L/β_NL + i_L/β_Ni,
        D(R_N) ~ -dec_R/β_NR + ϵ_tvr*tvr_B/β_NB + (1-κ_E)*syn_Enz/β_NE,
        u_PlantN ~ min(u_PlantNmax, k_PlantN*I_N), # TODO specify as Pin 
        D(I_N) ~ i_IN - u_PlantN - l_N*I_N + Φ_N,
        Φ_N ~ Φ_Nu + Φ_NB + r_tvr/β_NB,
        Φ_Nu ~ (1-ν_N) * u_NOM,
        u_N ~ ν_N * u_NOM + u_immNPot,
        u_NOM ~ dec_L/β_NL + dec_R/β_NR + κ_E*syn_Enz/β_NE,
        u_immNPot ~ i_BN * I_N,
        N_synBN ~ u_N - syn_Enz/β_NE,
        M_ImbN ~ u_N - (syn_B/β_NB + syn_Enz/β_NE),
        Φ_NB ~ M_ImbN - u_immNPot,
        # co-limitation
        #syn_B ~ C_synBC,
        syn_B ~ min(C_synBC, β_NB*N_synBN), # TODO add P limitation
        α_L ~ (1-α_R),
        # excess fluxes
        r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
        # auxilaries
        #D(cumresp) ~ r_B + r_tvr,
        r_tvr ~ (1-ϵ_tvr)*tvr_B,
        r_B ~ r_GEnz + r_G + r_M + r_O,
        r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
        # TODO make ifelse work with Unitful
        # r_G ~ IfElse.ifelse(syn_B.val > 0.0, (1-ϵ)/ϵ * syn_B, 0.0), 
        r_G ~ (1-ϵ)/ϵ * syn_B, 
        ]
    ODESystem(eqs; name)    
end










