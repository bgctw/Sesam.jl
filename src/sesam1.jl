function sesam1C(;name)
    # @parameters t [unit = uT]
    # D = Differential(t)

    # @parameters ϵ ϵ_tvr κ_E a_E [unit = uT^-1] m [unit = uT^-1] τ [unit = uT^-1] 
    # @parameters k_L [unit = uT^-1] k_R [unit = uT^-1] k_mN [unit = uQ*uT^-1] 
    # @parameters α_R0

    # # state variables
    # @variables B(t) [unit = uQ] L(t) [unit = uQ]  R(t) [unit = uQ] α_R(t) cumresp(t) [unit = uQ]

    # @variables α_L(t) tmpNoUnit(t)
    # # fluxes
    # @variables (begin
    #     syn_Enz(t),[unit = uQ/uT], r_M(t),[unit = uQ/uT], tvr_B(t),[unit = uQ/uT],
    #     dec_LPot(t),[unit = uQ/uT], dec_L(t),[unit = uQ/uT], dec_RPot(t),[unit = uQ/uT],
    #     dec_R(t),[unit = uQ/uT], u_C(t),[unit = uQ/uT],
    #     C_synBCt(t),[unit = uQ/uT], C_synBC(t),[unit = uQ/uT], r_G(t),[unit = uQ/uT],
    #     r_tvr(t),[unit = uQ/uT], syn_B(t),[unit = uQ/uT],
    #     r_B(t),[unit = uQ/uT], r_GEnz(t),[unit = uQ/uT], r_O(t),[unit = uQ/uT]
    # end)

    # #------------- N ---------------------
    # @variables L_N(t) [unit = uQ] R_N(t) [unit = uQ] I_N(t) [unit = uQ]
    # @variables (begin
    #     Φ_N(t),[unit = uQ/uT], Φ_Nu(t),[unit = uQ/uT], Φ_NB(t),[unit = uQ/uT], 
    #     u_PlantN(t),[unit = uQ/uT], u_NOM(t),[unit = uQ/uT], u_immNPot(t),[unit = uQ/uT], 
    #     u_N(t),[unit = uQ/uT], N_synBN(t),[unit = uQ/uT], M_ImbN(t),[unit = uQ/uT], 
    #     β_NL(t), β_NR(t) 
    # end)
    # @parameters β_NE β_NB l_N [unit = uT^-1] ν_N i_BN [unit = uT^-1] 

    # # need to be specified by coupled system:
    # @variables (begin
    #     i_L(t),[unit = uQ/uT], β_Ni(t), i_IN(t),[unit = uQ/uT],
    #     u_PlantNmax(t),[unit = uQ/uT], k_PlantN(t),[unit = uT^-1]
    # end) 

    # replace patttern to remove units: ,?\[unit[^\]]+\]
    @parameters t 
    D = Differential(t)

    @parameters ϵ ϵ_tvr κ_E a_E  m  τ  
    @parameters k_L  k_R  k_mN  
    @parameters α_R0

    # state variables
    @variables B(t)  L(t)   R(t)  α_R(t) cumresp(t) 
    @variables dB(t) dL(t) dR(t) dα_R(t) r_tot(t)

    @variables α_L(t) tmpNoUnit(t)
    # fluxes
    @variables (begin
        syn_Enz(t), r_M(t), tvr_B(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t), 
        r_B(t), r_GEnz(t), r_O(t)
    end)

    # need to be specified by coupled system:
    @variables i_L(t) syn_B(t)

    eqs = [
        D(B) ~ dB,
        dB ~ syn_B - tvr_B,
        D(L) ~ dL,
        dL ~ -dec_L + i_L,
        D(R) ~ dR,
        dR ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
        D(α_R) ~ dα_R,
        dα_R ~ 0.0,
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
        # co-limitation
        # define syn_B in coupled model
        #syn_B ~ min(C_synBC, β_NB*N_synBN), # TODO add P limitation
        α_L ~ (1-α_R),
        # excess fluxes
        r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
        # auxilaries
        #D(cumresp) ~ r_B + r_tvr,
        r_tvr ~ (1-ϵ_tvr)*tvr_B,
        r_B ~ r_GEnz + r_G + r_M + r_O,
        r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
        # TODO make ifelse work with Unitful
        r_G ~ IfElse.ifelse(syn_B.val > 0.0, (1-ϵ)/ϵ * syn_B, 0.0), 
        #r_G ~ (1-ϵ)/ϵ * syn_B, 
        r_tot ~ r_B + r_tvr,
        ]
    ODESystem(eqs; name)    
end

function sesam1N(;name)
    @parameters t 
    D = Differential(t)

    sts = @variables (begin
        L_N(t),  R_N(t),  I_N(t), 
        #
        dL_N(t),  dR_N(t),  dI_N(t),
        Φ_N(t), Φ_Nu(t), Φ_NB(t), 
        u_PlantN(t), u_NOM(t), u_immNPot(t), 
        u_N(t), N_synBN(t), M_ImbN(t), 
        β_NL(t), β_NR(t),
        leach_N(t),
        # need to be specified by coupled system:
        β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t)
    end)
    ps = @parameters β_NE β_NB l_N  ν_N i_BN  

    @named sC = sesam1C()
    @unpack L, R, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, syn_B, syn_Enz, r_tvr = sC
    @unpack κ_E = sC

    eqs = [
        # N cycle
        β_NL ~ L/L_N,
        β_NR ~ R/R_N,
        D(L_N) ~ dL_N,
        dL_N ~ -dec_L/β_NL + i_L/β_Ni,
        D(R_N) ~ dR_N,
        dR_N ~ -dec_R/β_NR + ϵ_tvr*tvr_B/β_NB + (1-κ_E)*syn_Enz/β_NE,
        D(I_N) ~ dI_N,
        u_PlantN ~ min(u_PlantNmax, k_PlantN*I_N), # TODO specify as Pin 
        dI_N ~ i_IN - u_PlantN - leach_N + Φ_N,
        leach_N ~ l_N*I_N,
        Φ_N ~ Φ_Nu + Φ_NB + r_tvr/β_NB,
        Φ_Nu ~ (1-ν_N) * u_NOM,
        u_N ~ ν_N * u_NOM + u_immNPot,
        u_NOM ~ dec_L/β_NL + dec_R/β_NR + κ_E*syn_Enz/β_NE,
        u_immNPot ~ i_BN * I_N,
        N_synBN ~ u_N - syn_Enz/β_NE,
        M_ImbN ~ u_N - (syn_B/β_NB + syn_Enz/β_NE),
        Φ_NB ~ M_ImbN - u_immNPot,
        ]
    extend(ODESystem(eqs, t, sts, ps; name), sC)
end

sesam1(args...;kwargs...) = sesam1N(args...;kwargs...)










