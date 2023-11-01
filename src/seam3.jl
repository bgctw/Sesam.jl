function seam3C(;name)
    @parameters t 
    D = Differential(t)

    ps = @parameters(
        ϵ, ϵ_tvr, κ_E, a_E,  m,  τ,  
        k_L,  k_R,  k_m, k_N,  
        #α_R0
         )    

    sts = @variables (begin
        B(t),  L(t),   R(t),  cumresp(t),
        E_L(t), E_R(t), tvr_Enz(t),
        dB(t), dL(t), dR(t), r_tot(t),
        dE_L(t), dE_R(t),
        syn_Enz(t), r_M(t), tvr_B(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t), 
        r_B(t), r_GEnz(t), r_O(t),
        α_L(t), α_R(t),
        # need to be defined by component across all elements
        syn_B(t), sum_w(t),
        ω_Enz(t), ω_L(t), ω_R(t),
        dα_R(t),
        # need to be specified by coupled system:
        i_L(t)      
    end)

    eqs = [
        D(B) ~ dB, dB ~ syn_B - tvr_B,
        D(L) ~ dL, dL ~ -dec_L + i_L,
        D(R) ~ dR, dR ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*tvr_Enz,
        D(E_L) ~ dE_L, dE_L ~ syn_Enz * α_L - k_N * E_L,
        D(E_R) ~ dE_R, dE_R ~ syn_Enz * α_R - k_N * E_R,
        syn_Enz ~ a_E*B,
        tvr_Enz ~ k_N * (E_L + E_R),
        r_M ~ m*B,
        tvr_B ~ τ*B,
        dec_LPot ~ k_L * L,
        dec_L ~ dec_LPot*(E_L)/(k_m + E_L),
        dec_RPot ~ k_R * R,
        dec_R ~ dec_RPot*(E_R)/(k_m + E_R),
        u_C ~ dec_L + dec_R + κ_E*tvr_Enz,
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
        ]
    ODESystem(eqs, t, sts, ps; name)    
end

# nitrogen dynamics as in sesam3N_revMM
seam3N(;name) = sesam3N_revMM(;name, sC = seam3C(name=:sC))

function seam3CN(;name, δ=40.0, max_w=12, use_seam_revenue=false)
    @parameters t 
    D = Differential(t)
    @named sN = seam3N()
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
        # need to limit, otherwise danger of Inf and nan after too long steps
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
        D(α_R) ~ dα_R, dα_R ~ (α_RT - α_R)*(τ + abs(syn_B)/B),
        ]
    extend(ODESystem(vcat(eqs,eqs_rev), t, vcat(sts, sts_rev), ps; name), sN)
end

seam3(args...;kwargs...) = seam3CN(args...;kwargs...)







