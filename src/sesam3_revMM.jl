# N-version used in seam, not affacted by further developments of sesam3?
function sesam3C_revMM(;name, k_N=60.0)
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
    #     u_NPot(t),[unit = uQ/uT], N_synBN(t),[unit = uQ/uT], M_ImbN(t),[unit = uQ/uT], 
    #     β_NL(t), β_NR(t) 
    # end)
    # @parameters β_NEnz β_NB l_N [unit = uT^-1] ν_N i_BN [unit = uT^-1] 

    # # need to be specified by coupled system:
    # @variables (begin
    #     i_L(t),[unit = uQ/uT], β_Ni(t), i_IN(t),[unit = uQ/uT],
    #     u_PlantNmax(t),[unit = uQ/uT], k_PlantN(t),[unit = uT^-1]
    # end) 

    # replace patttern to remove units: ,?\[unit[^\]]+\]
    @parameters t 
    D = Differential(t)

    ps = @parameters(
        ϵ, ϵ_tvr, κ_E, a_E,  m,  τ,  
        k_L,  k_R,  k_mN_L, k_mN_R, k_N=k_N,
        #ρ_CBtvr = 0.0, # proportion of carbon resorption on turnover
    )    

    sts = @variables (begin
        B(t),  L(t),   R(t),  cumresp(t),
        dB(t), dL(t), dR(t), r_tot(t),
        syn_Enz(t), tvr_Enz(t), r_M(t), tvr_B(t),
        E_L(t), E_R(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t), 
        r_B(t), r_GEnz(t), r_O(t),
        α_L(t), α_R(t),
        # equations need to be defined by component across all elements
        syn_B(t), sum_w(t),
        ω_Enz(t), ω_L(t), ω_R(t),
        dα_R(t),        
        # need to be specified by coupled system:
        i_L(t)
    end)

    eqs = [
        D(B) ~ dB, dB ~ syn_B - tvr_B,
        D(L) ~ dL, dL ~ -dec_L + i_L,
        D(R) ~ dR, dR ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
        syn_Enz ~ a_E*B, tvr_Enz ~ syn_Enz,
        r_M ~ m*B,
        tvr_B ~ τ*B,
        dec_LPot ~ k_L * L,
        dec_L ~ dec_LPot*(α_L * syn_Enz)/(k_mN_L + α_L*syn_Enz),
        dec_RPot ~ k_R * R,
        dec_R ~ dec_RPot*(α_R * syn_Enz)/(k_mN_R + α_R*syn_Enz),
        u_C ~ dec_L + dec_R + κ_E*syn_Enz,
        C_synBCt ~ u_C - syn_Enz/ϵ - r_M,
        #C_synBC ~ IfElse.ifelse(C_synBCt > 0.0, ϵ*C_synBCt, C_synBCt), 
        C_synBC ~ C_synBCt - ((C_synBCt > 0.0) * (1-ϵ)*C_synBCt), 
        #syn_B ~ C_synBC # define in all-elements model
        r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
        #D(cumresp) ~ r_B + r_tvr,
        r_tvr ~ (1-ϵ_tvr)*tvr_B,
        r_B ~ r_GEnz + r_G + r_M + r_O,
        r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
        #r_G ~ IfElse.ifelse(syn_B.val > 0.0, (1-ϵ)/ϵ * syn_B, 0.0), 
        r_G ~ (syn_B.val > 0.0) * (1-ϵ)/ϵ * syn_B, 
        r_tot ~ r_B + r_tvr,
        E_L ~ (α_L * syn_Enz)/k_N,
        E_R ~ (α_R * syn_Enz)/k_N,
        ]
    ODESystem(eqs, t, sts, ps; name)    
end

function sesam3N_revMM(;name, sC = sesam3C_revMM(name=:sC))
    @parameters t 
    D = Differential(t)

    sts = @variables (begin
        L_N(t),  R_N(t),  I_N(t), 
        dL_N(t),  dR_N(t),  dI_N(t),
        Φ_N(t), Φ_Nu(t), Φ_NB(t), Φ_Ntvr(t),
        u_PlantN(t), u_NOM(t), u_immNPot(t), 
        u_NPot(t), N_synBN(t), M_ImbN(t), 
        β_NL(t), β_NR(t),
        leach_N(t),
        p_uNmic(t), ν_TN(t),
        # need to be specified by coupled system:
        β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t)
    end)
    ps = @parameters β_NEnz β_NB l_N  ν_N i_BN  

    @unpack L, R, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, syn_B, syn_Enz, tvr_Enz, r_tvr, κ_E = sC

    eqs = [
        β_NL ~ L/L_N, β_NR ~ R/R_N,
        D(L_N) ~ dL_N, dL_N ~ -dec_L/β_NL + i_L/β_Ni,
        D(R_N) ~ dR_N, dR_N ~ -dec_R/β_NR + ϵ_tvr*tvr_B/β_NB + (1-κ_E)*tvr_Enz/β_NEnz,
        D(I_N) ~ dI_N,
        u_PlantN ~ min(u_PlantNmax, k_PlantN*I_N), 
        dI_N ~ i_IN - u_PlantN - leach_N + Φ_N,
        leach_N ~ l_N*I_N,
        Φ_N ~ Φ_Nu + Φ_NB + Φ_Ntvr,
        Φ_Ntvr ~ r_tvr/β_NB,
        Φ_Nu ~ (1-ν_N) * u_NOM,
        u_NPot ~ ν_N * u_NOM + u_immNPot,
        u_NOM ~ dec_L/β_NL + dec_R/β_NR + κ_E*tvr_Enz/β_NEnz,
        u_immNPot ~ i_BN * I_N,
        N_synBN ~ u_NPot - syn_Enz/β_NEnz,
        M_ImbN ~ u_NPot - (syn_B/β_NB + syn_Enz/β_NEnz),
        Φ_NB ~ M_ImbN - u_immNPot,
        p_uNmic ~ u_immNPot/(u_immNPot + u_PlantN),
        ν_TN ~ ν_N+(1-ν_N)*p_uNmic,        
        ]
    extend(ODESystem(eqs, t, sts, ps; name), sC)
end

function sesam3CN_revMM(;name, δ=40.0, max_w=12, use_seam_revenue=false, sN=sesam3N_revMM(name=:sN))
    @parameters t 
    D = Differential(t)
    @unpack α_L, α_R, syn_B, B, C_synBC, β_NB, N_synBN, tvr_B, τ, ϵ = sN
    @unpack β_NB, β_NEnz, β_NL, β_NR = sN
    @unpack ω_Enz, ω_L, ω_R = sN
    @unpack u_immNPot, u_PlantN, ν_N, ν_TN = sN    
    sts = @variables (begin
        C_synBN(t), 
        C_synBmC(t), C_synBmN(t),
        dα_L(t), dα_R(t),
        w_C(t), w_N(t), lim_C(t), lim_N(t)
    end)
    ps = @parameters δ=δ
    eqs_rev, sts_rev = use_seam_revenue ? 
        get_dα_eq_seam(sN) : get_dα_eq_sesam3CN(sN)
    @variables α_LT(t) α_RT(t)
    lim_E = SA[lim_C, lim_N]
    β_B = SA[1.0, β_NB]
    ν_TZ = SA[ϵ, ν_TN] 
    eqs = [
        C_synBN ~ β_NB*N_synBN,
        syn_B ~ min(C_synBC, C_synBN), 
        C_synBmC ~ min(C_synBN), 
        C_synBmN ~ min(C_synBC), 
        # need minimum, otherwise danger of Inf and nan -> instability
        w_C ~ exp(min(max_w, -δ/tvr_B*(C_synBC - syn_B))),
        w_N ~ exp(min(max_w, -δ/tvr_B*(C_synBN - syn_B))),
        # w_C ~ min(max_w, exp(δ/tvr_B*(C_synBmC - C_synBC))),
        # w_N ~ min(max_w, exp(δ/tvr_B*(C_synBmN - C_synBN))),
        # w_C ~ exp(δ/tvr_B*(C_synBmC - C_synBC)),
        # w_N ~ exp(δ/tvr_B*(C_synBmN - C_synBN)),
        lim_C ~ w_C/(w_C + w_N), lim_N ~ w_N/(w_C + w_N), # normalized for plot
        # α_LT, α_RT by get_dα_eq_X
        ω_Enz ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NEnz], β_B),
        ω_L ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NL], β_B, ν_TZ),
        ω_R ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NR], β_B, ν_TZ),
        D(α_L) ~ dα_L, dα_L ~ (α_LT - α_L)*(τ + abs(syn_B)/B),
        D(α_R) ~ dα_R #, dα_R ~ (α_RT - α_R)*(τ + abs(syn_B)/B),
        ]
    extend(ODESystem(vcat(eqs,eqs_rev), t, vcat(sts, sts_rev), ps; name), sN)
end

sesam3_revMM(args...;kwargs...) = sesam3CN_revMM(args...;kwargs...)


function get_dα_eq_seam(sN)
    # needs dec_RPot, which is specific to revMM decomposition formulation
    @parameters t 
    @unpack dec_LPot, dec_RPot, k_mN_L, k_mN_R, syn_Enz, α_L, α_R, β_NL, β_NR = sN
    @unpack τ, syn_B, B = sN
    sts = @variables (begin
        α_LT(t), α_RT(t),
        rev_LC(t), rev_RC(t), rev_LN(t), rev_RN(t), 
        α_RC(t), α_RN(t) 
    end)
    # need to be defined in coupled component:
    @variables w_C(t), w_N(t)
    eqs = [
        rev_LC ~ dec_LPot/(k_mN_L + α_L*syn_Enz),
        rev_RC ~ dec_RPot/(k_mN_R + α_R*syn_Enz),
        rev_LN ~ rev_LC/β_NL,
        rev_RN ~ rev_RC/β_NR,
        α_RC ~ rev_RC/(rev_RC + rev_LC),
        α_RN ~ rev_RN/(rev_RN + rev_LN),
        α_RT ~ (w_C * α_RC + w_N * α_RN)/(w_C + w_N),
        α_LT ~ 1.0 - α_RT,
        dα_R ~ (α_RT - α_R)*(τ + abs(syn_B)/B),
    ]
    (;eqs, sts)
end








