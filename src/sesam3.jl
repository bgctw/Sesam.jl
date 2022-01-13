function sesam3C(;name)
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
    # @parameters β_NEnz β_NB l_N [unit = uT^-1] ν_N i_BN [unit = uT^-1] 

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

    @variables (begin
        B(t),  L(t),   R(t),  cumresp(t),
        dB(t), dL(t), dR(t), r_tot(t),
        syn_Enz(t), r_M(t), tvr_B(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t), 
        r_B(t), r_GEnz(t), r_O(t),
        # need to be defined by component across all elements
        α_L(t), α_R(t),
        # need to be specified by coupled system:
        i_L(t), syn_B(t) 
    end)

    eqs = [
        D(B) ~ dB, dB ~ syn_B - tvr_B,
        D(L) ~ dL, dL ~ -dec_L + i_L,
        D(R) ~ dR, dR ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
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
        #syn_B ~ C_synBC # define in all-elements model
        r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
        #D(cumresp) ~ r_B + r_tvr,
        r_tvr ~ (1-ϵ_tvr)*tvr_B,
        r_B ~ r_GEnz + r_G + r_M + r_O,
        r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
        r_G ~ IfElse.ifelse(syn_B.val > 0.0, (1-ϵ)/ϵ * syn_B, 0.0), 
        r_tot ~ r_B + r_tvr,
        ]
    ODESystem(eqs; name)    
end

function sesam3N(;name)
    @parameters t 
    D = Differential(t)

    sts = @variables (begin
        L_N(t),  R_N(t),  I_N(t), 
        dL_N(t),  dR_N(t),  dI_N(t),
        Φ_N(t), Φ_Nu(t), Φ_NB(t), Φ_Ntvr(t),
        u_PlantN(t), u_NOM(t), u_immNPot(t), 
        u_N(t), N_synBN(t), M_ImbN(t), 
        β_NL(t), β_NR(t),
        leach_N(t),
        # need to be specified by coupled system:
        β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t)
    end)
    ps = @parameters β_NEnz β_NB l_N  ν_N i_BN  

    @named sC = sesam3C()
    @unpack L, R, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, syn_B, syn_Enz, r_tvr, κ_E = sC

    eqs = [
        β_NL ~ L/L_N, β_NR ~ R/R_N,
        D(L_N) ~ dL_N, dL_N ~ -dec_L/β_NL + i_L/β_Ni,
        D(R_N) ~ dR_N, dR_N ~ -dec_R/β_NR + ϵ_tvr*tvr_B/β_NB + (1-κ_E)*syn_Enz/β_NEnz,
        D(I_N) ~ dI_N,
        u_PlantN ~ min(u_PlantNmax, k_PlantN*I_N), 
        dI_N ~ i_IN - u_PlantN - leach_N + Φ_N,
        leach_N ~ l_N*I_N,
        Φ_N ~ Φ_Nu + Φ_NB + Φ_Ntvr,
        Φ_Ntvr ~ r_tvr/β_NB,
        Φ_Nu ~ (1-ν_N) * u_NOM,
        u_N ~ ν_N * u_NOM + u_immNPot,
        u_NOM ~ dec_L/β_NL + dec_R/β_NR + κ_E*syn_Enz/β_NEnz,
        u_immNPot ~ i_BN * I_N,
        N_synBN ~ u_N - syn_Enz/β_NEnz,
        M_ImbN ~ u_N - (syn_B/β_NB + syn_Enz/β_NEnz),
        Φ_NB ~ M_ImbN - u_immNPot,
        ]
    extend(ODESystem(eqs, t, sts, ps; name), sC)
end


function get_revenue_eq_seam(sN)
    @parameters t 
    @unpack dec_LPot, dec_RPot, k_mN, syn_Enz, α_L, α_R, β_NL, β_NR = sN
    sts = @variables (begin
        α_LT(t), α_RT(t),
        rev_LC(t), rev_RC(t), rev_LN(t), rev_RN(t), 
        α_RC(t), α_RN(t) 
    end)
    # need to be defined in coupled component:
    @variables w_C(t), w_N(t)
    eqs = [
        rev_LC ~ dec_LPot/(k_mN + α_L*syn_Enz),
        rev_RC ~ dec_RPot/(k_mN + α_R*syn_Enz),
        rev_LN ~ rev_LC/β_NL,
        rev_RN ~ rev_RC/β_NR,
        α_RC ~ rev_RC/(rev_RC + rev_LC),
        α_RN ~ rev_RN/(rev_RN + rev_LN),
        α_RT ~ (w_C * α_RC + w_N * α_RN)/(w_C + w_N),
        α_LT ~ 1.0 - α_RT,
    ]
    (;eqs, sts)
end

function get_revenue_eq_sesam3CN(sN)
    @parameters t 
    @unpack α_L, α_R, dec_L, dec_R, β_NL, β_NR, β_NEnz, syn_Enz = sN
    sts = @variables (begin
        α_LT(t), α_RT(t),
        invest_L(t), invest_R(t), return_L(t), return_R(t), revenue_L(t), revenue_R(t),
        invest_Ln(t), invest_Rn(t), return_Ln(t), return_Rn(t), 
        revenue_sum(t)
    end)
    # need to be defined in coupled component:
    @variables w_C(t), w_N(t)
    eqs = [
        invest_L ~ α_L*syn_Enz*(w_C + w_N/β_NEnz),
        invest_R ~ α_R*syn_Enz*(w_C + w_N/β_NEnz),
        invest_Ln ~ invest_L/(invest_L + invest_R),
        invest_Rn ~ invest_R/(invest_L + invest_R),
        return_L ~ dec_L * (w_C + w_N/β_NL), 
        return_R ~ dec_R * (w_C + w_N/β_NR), 
        return_Ln ~ return_L/(return_L + return_R),
        return_Rn ~ return_R/(return_L + return_R),
        revenue_L ~ return_L / invest_L,
        revenue_R ~ return_R / invest_R,
        revenue_sum ~ revenue_L + revenue_R,
        α_LT ~ revenue_L/revenue_sum,
        α_RT ~ revenue_R/revenue_sum,
        ]
    (;eqs, sts)
end


function sesam3CN(;name, δ=20.0, max_w=1e5, use_seam_revenue=false)
    @parameters t 
    D = Differential(t)
    @named sN = sesam3N()
    @unpack α_L, α_R, syn_B, B, C_synBC, β_NB, N_synBN, tvr_B, τ = sN
    sts = @variables (begin
        C_synBN(t), 
        C_synBmC(t), C_synBmN(t),
        dα_L(t), dα_R(t),
        w_C(t), w_N(t), lim_C(t), lim_N(t)
    end)
    ps = @parameters δ=δ
    eqs_rev, sts_rev = use_seam_revenue ? 
        get_revenue_eq_seam(sN) : get_revenue_eq_sesam3CN(sN)
    @variables α_LT(t) α_RT(t)
    eqs = [
        C_synBN ~ β_NB*N_synBN,
        syn_B ~ min(C_synBC, C_synBN), 
        C_synBmC ~ min(C_synBN), 
        C_synBmN ~ min(C_synBC), 
        w_C ~ min(max_w, exp(δ/tvr_B*(C_synBmC - C_synBC))),
        w_N ~ min(max_w, exp(δ/tvr_B*(C_synBmN - C_synBN))),
        lim_C ~ w_C/(w_C + w_N), lim_N ~ w_N/(w_C + w_N),
        # α_LT, α_RT by get_revenue_eq_X
        D(α_L) ~ dα_L, dα_L ~ (α_LT - α_L)*(τ + abs(syn_B)/B),
        D(α_R) ~ dα_R, dα_R ~ (α_RT - α_R)*(τ + abs(syn_B)/B),
        ]
    extend(ODESystem(vcat(eqs,eqs_rev), t, vcat(sts, sts_rev), ps; name), sN)
end

sesam3(args...;kwargs...) = sesam3CN(args...;kwargs...)

"""
R pool is a mixture of microbial turnover and enzymes.
Here the stable C/N ratio is computed based on given parameterization.
"""
function calculate_β_NR_sesam3(p,s)
    w_REnz = p[s.a_E]*(1-p[s.κ_E]) # a_E (1-κ_E) B
    w_RB = p[s.ϵ_tvr]*p[s.τ]       # τ ϵ_tvr B
    β_NR = 1/((w_REnz/p[s.β_NEnz] + w_RB/p[s.β_NB])/(w_REnz + w_RB))
end








