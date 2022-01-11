function sesam1(;name)
    @parameters t 
    D = Differential(t)

    @parameters ϵ_tvr κ_E a_E m τ k_L k_R k_mN ϵ
    @parameters α_R0
    @variables B(t) L(t) R(t) α_R(t) cumresp(t)
    @variables syn_Enz(t) r_M(t) tvr_B(t) dec_LPot(t) dec_L(t) dec_RPot(t) dec_R(t) u_C(t)
    @variables C_synBCt(t) C_synBC(t) r_G(t) r_tvr(t) syn_B(t)
    @variables C_tot(t) r_B(t) r_GEnz(t) r_O(t)
    @variables α_L(t)

    @variables L_N(t) R_N(t) I_N(t)
    @variables Φ_N(t) Φ_Nu(t) Φ_NB(t) u_PlantN(t) u_NOM(t) u_immNPot(t) 
    @variables u_N(t) N_synBN(t) M_ImbN(t) β_NL(t) β_NR(t) 
    @parameters β_NE β_NB l_N ν_N i_BN 

    # need to be specified by coupled system:
    @variables i_L(t) β_Ni(t) i_IN(t)
    @variables u_PlantNmax(t) k_PlantN(t)

    eqs = [
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
        C_synBC ~ IfElse.ifelse(C_synBCt > 0.0, ϵ*C_synBCt, C_synBCt), 
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
        # excell fluxes
        r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
        # auxilaries
        #D(cumresp) ~ r_B + r_tvr,
        r_tvr ~ (1-ϵ_tvr)*tvr_B,
        r_B ~ r_GEnz + r_G + r_M + r_O,
        r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
        r_G ~ IfElse.ifelse(syn_B > 0, (1-ϵ)/ϵ * syn_B, 0.0), 
        ]
    ODESystem(eqs; name)    
end










