function sesam3C(;name, k_N=60.0)
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

    @parameters ϵ ϵ_tvr κ_E a_E  m  τ  
    @parameters k_L  k_R  k_mN_L k_mN_R k_N=k_N
    @parameters α_R0
    @parameters ρ_CBtvr = 0.0 # proportion of carbon resorption on turnover

    @variables (begin
        B(t),  L(t),   R(t),  cumresp(t),
        dB(t), dL(t), dR(t), r_tot(t),
        syn_Enz(t), tvr_Enz(t), r_M(t), tvr_B(t),
        E_L(t), E_R(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t), 
        r_B(t), r_GEnz(t), r_O(t),
        tvr_B0(t), resorp_C(t), 
        # need to be defined by component across all elements
        α_L(t), α_R(t), 
        # need to be specified by coupled system:
        i_L(t), syn_B(t)
    end)

    eqs = [
        D(B) ~ dB, dB ~ syn_B - tvr_B,
        D(L) ~ dL, dL ~ -dec_L + i_L,
        D(R) ~ dR, dR ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
        syn_Enz ~ a_E*B, tvr_Enz ~ syn_Enz,
        r_M ~ m*B,
        # element resorption changes stoichiometry: β_EB -> β_EBtvr
        tvr_B ~ tvr_B0 - resorp_C,
        tvr_B0 ~ τ*B, # turnover before resorption 
        resorp_C ~ ρ_CBtvr * tvr_B0,
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
    ODESystem(eqs; name)    
end

function sesam3N(;name, sC = sesam3C(name=:sC))
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
        resorp_N(t), β_NBtvr(t),
        # need to be specified by coupled system:
        β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t),
        SOC(t), SON(t), β_NSOM(t), BtoSOC(t)
    end)
    ps = @parameters β_NEnz β_NB l_N  ν_N i_BN ρ_NBtvr=0.0

    @unpack L, R, B, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, tvr_B0, syn_B, syn_Enz, tvr_Enz, r_tvr, κ_E = sC

    eqs = [
        β_NL ~ L/L_N, β_NR ~ R/R_N,
        D(L_N) ~ dL_N, dL_N ~ -dec_L/β_NL + i_L/β_Ni,
        D(R_N) ~ dR_N, dR_N ~ -dec_R/β_NR + ϵ_tvr*tvr_B/β_NBtvr + (1-κ_E)*tvr_Enz/β_NEnz,
        D(I_N) ~ dI_N,
        u_PlantN ~ min(u_PlantNmax, k_PlantN*I_N), 
        dI_N ~ i_IN - u_PlantN - leach_N + Φ_N,
        leach_N ~ l_N*I_N,
        Φ_N ~ Φ_Nu + Φ_NB + Φ_Ntvr,
        Φ_Ntvr ~ r_tvr/β_NBtvr,
        Φ_Nu ~ (1-ν_N) * u_NOM,
        u_NPot ~ ν_N * u_NOM + u_immNPot + resorp_N,
        u_NOM ~ dec_L/β_NL + dec_R/β_NR + κ_E*tvr_Enz/β_NEnz,
        u_immNPot ~ i_BN * I_N,
        N_synBN ~ u_NPot - syn_Enz/β_NEnz,
        M_ImbN ~ u_NPot - (syn_B/β_NB + syn_Enz/β_NEnz),
        Φ_NB ~ M_ImbN - u_immNPot,
        # make sure ρ_NBtvr >= 0
        resorp_N ~ ρ_NBtvr * tvr_B0/β_NB,
        β_NBtvr ~ tvr_B / (tvr_B0/β_NB - resorp_N),
        # observables for diagnostic output
        SOC ~ L + R + B,
        SON ~ L_N + R_N + B/β_NB,
        β_NSOM ~ SOC/SON,
        BtoSOC ~ B/SOC,
    ]
    extend(ODESystem(eqs, t, sts, ps; name), sC)
end

function get_revenue_eq_sesam3CN(sN)
    @parameters t 
    @unpack α_L, α_R, dec_L, dec_R, β_NL, β_NR, β_NEnz, syn_Enz = sN
    @unpack u_immNPot, u_PlantN, ν_N = sN
    sts = @variables (begin
        α_LT(t), α_RT(t),
        p_uNmic(t), p_oNmic(t),         
        invest_L(t), invest_R(t), return_L(t), return_R(t), revenue_L(t), revenue_R(t),
        #invest_Ln(t), invest_Rn(t), return_Ln(t), return_Rn(t), 
        revenue_sum(t)
    end)
    # need to be defined in coupled component:
    @variables w_C(t), w_N(t)
    eqs = [
        p_uNmic ~ u_immNPot/(u_immNPot + u_PlantN),
        p_oNmic ~ ν_N+(1-ν_N)*p_uNmic,        
        invest_L ~ α_L*syn_Enz*(w_C + w_N/β_NEnz),
        invest_R ~ α_R*syn_Enz*(w_C + w_N/β_NEnz),
        # return_L ~ dec_L * (w_C + w_N/β_NL), 
        # return_R ~ dec_R * (w_C + w_N/β_NR), 
        return_L ~ dec_L * (w_C + w_N/β_NL*p_oNmic), 
        return_R ~ dec_R * (w_C + w_N/β_NR*p_oNmic),         
        revenue_L ~ return_L / invest_L,
        revenue_R ~ return_R / invest_R,
        revenue_sum ~ revenue_L + revenue_R,
        α_LT ~ revenue_L/revenue_sum,
        α_RT ~ revenue_R/revenue_sum,
        # auxilary for plotting
        # invest_Ln ~ invest_L/(invest_L + invest_R),
        # invest_Rn ~ invest_R/(invest_L + invest_R),
        # return_Ln ~ return_L/(return_L + return_R),
        # return_Rn ~ return_R/(return_L + return_R),
        ]
    (;eqs, sts)
end

function sesam3CN(;name, δ=40.0, max_w=12, use_seam_revenue=false, sN=sesam3N(name=:sN))
    @parameters t 
    D = Differential(t)
    @unpack α_L, α_R, syn_B, B, C_synBC, β_NB, N_synBN, tvr_B, τ = sN
    @unpack u_immNPot, u_PlantN, ν_N = sN    
    sts = @variables (begin
        C_synBN(t), 
        #dα_L(t), 
        dα_R(t),
        w_C(t), w_N(t), lim_C(t), lim_N(t)
    end)
    ps = @parameters δ=δ
    eqs_rev, sts_rev = use_seam_revenue ? 
        get_revenue_eq_seam(sN) : get_revenue_eq_sesam3CN(sN)
    @variables α_LT(t) α_RT(t)
    eqs = [
        C_synBN ~ β_NB*N_synBN,
        syn_B ~ min(C_synBC, C_synBN), 
        # need minimum, otherwise danger of Inf and nan -> instability
        w_C ~ exp(min(max_w, -δ/tvr_B*(C_synBC - syn_B))),
        w_N ~ exp(min(max_w, -δ/tvr_B*(C_synBN - syn_B))),
        lim_C ~ w_C/(w_C + w_N), lim_N ~ w_N/(w_C + w_N), # normalized for plot
        # α_LT, α_RT by get_revenue_eq_X
        #D(α_L) ~ dα_L, dα_L ~ (α_LT - α_L)*(τ + abs(syn_B)/B),
        α_L ~ 1 - α_R,
        D(α_R) ~ dα_R, dα_R ~ (α_RT - α_R)*(τ + abs(syn_B)/B),
        ]
    extend(ODESystem(vcat(eqs,eqs_rev), t, vcat(sts, sts_rev), ps; name), sN)
end

# sesam3 now refers to sesam3CNP
#sesam3(args...;kwargs...) = sesam3CN(args...;kwargs...)

"""
R pool is a mixture of microbial turnover and enzymes.
Here the stable C/N ratio is computed based on given parameterization.
"""
function calculate_β_NR_sesam3(p,s)
    w_REnz = p[s.a_E]*(1-p[s.κ_E]) # a_E (1-κ_E) B
    w_RBtvr = p[s.ϵ_tvr]*p[s.τ]       # τ ϵ_tvr B
    # β_NBtvr may differ from β_NB because of elemental resorption
    # if ρ is not in parameter dictionary, its default is zero
    ρ_CBtvr = (hasproperty(s, :ρ_CBtvr) && haskey(p, s.ρ_CBtvr)) ? p[s.ρ_CBtvr] : 0.0
    ρ_NBtvr = (hasproperty(s, :ρ_NBtvr) && haskey(p, s.ρ_NBtvr)) ? p[s.ρ_NBtvr] : 0.0
    β_NBtvr = p[s.β_NB] * (1-ρ_CBtvr)/(1-ρ_NBtvr)
    β_NR = 1/((w_REnz/p[s.β_NEnz] + w_RBtvr/β_NBtvr)/(w_REnz + w_RBtvr))
end

function calculate_β_NR_sesam3(p_sym)
    w_REnz = p_sym[:s₊a_E]*(1-p_sym[:s₊κ_E]) # a_E (1-κ_E) B
    w_RBtvr = p_sym[:s₊ϵ_tvr]*p_sym[:s₊τ]       # τ ϵ_tvr B
    # β_NBtvr may differ from β_NB because of elemental resorption
    # if ρ is not in parameter dictionary, its default is zero
    ρ_CBtvr = (haskey(p_sym, :s₊ρ_CBtvr)) ? p_sym[:s₊ρ_CBtvr] : 0.0
    ρ_NBtvr = (haskey(p_sym, :s₊ρ_NBtvr)) ? p_sym[:s₊ρ_NBtvr] : 0.0
    β_NBtvr = p_sym[:s₊β_NB] * (1-ρ_CBtvr)/(1-ρ_NBtvr)
    β_NR = 1/((w_REnz/p_sym[:s₊β_NEnz] + w_RBtvr/β_NBtvr)/(w_REnz + w_RBtvr))
end







