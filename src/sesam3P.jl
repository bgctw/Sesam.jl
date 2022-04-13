function sesam3P(;name, sN = sesam3N(name=:sN))
    @parameters t 
    D = Differential(t)

    sts = @variables (begin
        L_P(t),  R_P(t),  I_P(t), 
        dL_P(t),  dR_P(t),  dI_P(t),
        Φ_P(t), Φ_Pu(t), Φ_PB(t), Φ_Ptvr(t),
        u_PlantP(t), u_POM(t), u_immPPot(t), 
        u_PPot(t), P_synBP(t), M_ImbP(t), 
        β_PL(t), β_PR(t),
        leach_P(t),
        α_LP(t), 
        α_RP(t),
        dec_LP(t), dec_RP(t), dec_LPPlant(t), dec_RPPlant(t),
        # need to be specified by coupled system:
        β_Pi(t), i_IP(t),
        u_PlantPmax(t), k_PlantP(t),
        s_EP(t), pL_sEP(t) # synthesis of E_LP adn E_RP enzymes by plants
    end)
    ps = @parameters β_PEnz β_PB l_P  ν_P i_BN i_BP k_LP k_RP 

    @unpack L, R, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, syn_B, syn_Enz, tvr_Enz, r_tvr, κ_E = sN
    @unpack k_mN_L, k_mN_R = sN

    eqs = [
        β_PL ~ L/L_P, β_PR ~ R/R_P,
        D(L_P) ~ dL_P, dL_P ~ -dec_L/β_PL - dec_LP + i_L/β_Pi,
        D(R_P) ~ dR_P, dR_P ~ -dec_R/β_PR - dec_RP + ϵ_tvr*tvr_B/β_PB + (1-κ_E)*tvr_Enz/β_PEnz,
        D(I_P) ~ dI_P,
        # assume that proportion of s_EP to biomineralizing L is 
        # proportional to potential biomineralization flux
        pL_sEP ~ k_LP*L_P/(k_LP*L_P + k_RP*R_P),
        dec_LP ~ k_LP * L_P * (α_LP*syn_Enz + pL_sEP*s_EP)/(k_mN_L + α_LP*syn_Enz + pL_sEP*s_EP),
        dec_LPPlant ~ k_LP * L_P * (pL_sEP*s_EP)/(k_mN_L + pL_sEP*s_EP), 
        # while_RP is decomp due to the sum of enzmyes of mic and plant (assumed to be the same with same k_M)
        # dec_RPPlant is the decomp that would be by enzymes levels by solely plant production, used in revenue calculation
        dec_RP ~ k_RP * R_P * (α_RP*syn_Enz + (1-pL_sEP)*s_EP)/(k_mN_R + α_RP*syn_Enz + (1-pL_sEP)*s_EP),
        dec_RPPlant ~ k_RP * R_P * ((1-pL_sEP)*s_EP)/(k_mN_R + (1-pL_sEP)*s_EP), 
        u_PlantP ~ min(u_PlantPmax, k_PlantP*I_P), 
        dI_P ~ i_IP - u_PlantP - leach_P  + dec_RP + dec_LP + Φ_P,
        leach_P ~ l_P*I_P,
        Φ_P ~ Φ_Pu + Φ_PB + Φ_Ptvr,
        Φ_Ptvr ~ r_tvr/β_PB,
        Φ_Pu ~ (1-ν_P) * u_POM,
        # dec_RP is biomineralization - cleaved P is all mineralized
        # u_POM is from depolymerizing - only a fraction ν_P is mineralized
        u_PPot ~ ν_P * u_POM + u_immPPot,
        u_POM ~ dec_L/β_PL + dec_R/β_PR + κ_E*tvr_Enz/β_PEnz,
        u_immPPot ~ i_BP * I_P,
        P_synBP ~ u_PPot - syn_Enz/β_PEnz,
        M_ImbP ~ u_PPot - (syn_B/β_PB + syn_Enz/β_PEnz),
        Φ_PB ~ M_ImbP - u_immPPot,
        ]
    extend(ODESystem(eqs, t, sts, ps; name), sN)
end

function get_revenue_eq_sesam3CNP(sP)
    @parameters t 
    @unpack α_L, α_R, dec_L, dec_R, β_NL, β_NR, β_NEnz, syn_Enz = sP
    @unpack β_PL, β_PR, β_PEnz = sP
    @unpack α_LP, α_RP, dec_LP, dec_RP, dec_RPPlant, dec_LPPlant = sP
    @unpack u_immPPot, u_PlantP = sP
    sts = @variables (begin
        α_LT(t), α_RT(t),
        α_LPT(t), α_RPT(t),
        syn_Enz_w(t), 
        p_uPmic(t),
        return_L(t), return_R(t), revenue_L(t), revenue_R(t),
        return_LP(t), return_RP(t), revenue_LP(t), revenue_RP(t),
        #invest_Ln(t), invest_Rn(t), return_Ln(t), return_Rn(t), 
        revenue_sum(t)
    end)
    # need to be defined in coupled component:
    @variables w_C(t), w_N(t), w_P(t)
    eqs = [
        syn_Enz_w ~ syn_Enz*(w_C + w_N/β_NEnz + w_P/β_PEnz),
        return_L ~ dec_L * (w_C + w_N/β_NL + w_P/β_PL), 
        return_R ~ dec_R * (w_C + w_N/β_NR + w_P/β_PR), 
        # return of enzymes produced in addition to that of plants
        # only proportion of the biomineralization flux ends up in microbes: part it taken up by plant
        # dec_LP is already in P units, no need to devide by β_P
        p_uPmic ~ u_immPPot/(u_immPPot + u_PlantP),
        return_LP ~ (dec_LP - dec_LPPlant) * p_uPmic * w_P, 
        return_RP ~ (dec_RP - dec_RPPlant) * p_uPmic * w_P, 
        revenue_L ~ return_L / (α_L * syn_Enz_w),
        revenue_R ~ return_R / (α_R * syn_Enz_w),
        revenue_LP ~ return_LP / (α_LP * syn_Enz_w),
        revenue_RP ~ return_RP / (α_RP * syn_Enz_w),
        revenue_sum ~ revenue_L + revenue_R + revenue_LP + revenue_RP,
        α_LT ~ revenue_L/revenue_sum,
        α_RT ~ revenue_R/revenue_sum,
        α_LPT ~ revenue_LP/revenue_sum,
        α_RPT ~ revenue_RP/revenue_sum,
        # auxilary for plotting
        # invest_Ln ~ invest_L/(invest_L + invest_R),
        # invest_Rn ~ invest_R/(invest_L + invest_R),
        # return_Ln ~ return_L/(return_L + return_R),
        # return_Rn ~ return_R/(return_L + return_R),
        ]
    (;eqs, sts)
end


function sesam3CNP(;name, δ=20.0, max_w=1e5, use_seam_revenue=false, sP=sesam3P(name=:sP))
    @parameters t 
    D = Differential(t)
    @unpack α_L, α_R, syn_B, B, C_synBC, β_NB, N_synBN, tvr_B, τ = sP
    @unpack β_PB, P_synBP, α_LP, α_RP = sP
    sts = @variables (begin
        C_synBmC(t), 
        C_synBN(t), C_synBmN(t),
        C_synBP(t), C_synBmP(t),
        #dα_L(t), 
        dα_R(t),
        dα_LP(t), dα_RP(t),
        w_C(t), w_N(t), lim_C(t), lim_N(t),
        sum_w(t),
        w_P(t), lim_P(t)
    end)
    ps = @parameters δ=δ
    eqs_rev, sts_rev = get_revenue_eq_sesam3CNP(sP)
    @variables α_LT(t) α_RT(t)
    @variables α_LPT(t) α_RPT(t)
    eqs = [
        C_synBN ~ β_NB*N_synBN,
        C_synBP ~ β_PB*P_synBP,
        syn_B ~ min(C_synBC, C_synBN, C_synBP), 
        C_synBmC ~ min(C_synBN, C_synBP), 
        C_synBmN ~ min(C_synBC, C_synBP), 
        C_synBmP ~ min(C_synBC, C_synBN), 
        # need minimum, otherwise danger of Inf and nan -> instability
        w_C ~ min(max_w, exp(δ/tvr_B*(C_synBmC - C_synBC))),
        w_N ~ min(max_w, exp(δ/tvr_B*(C_synBmN - C_synBN))),
        w_P ~ min(max_w, exp(δ/tvr_B*(C_synBmP - C_synBP))),
        # w_C ~ exp(δ/tvr_B*(C_synBmC - C_synBC)),
        # w_N ~ exp(δ/tvr_B*(C_synBmN - C_synBN)),
        sum_w ~ w_C + w_N + w_P,
        lim_C ~ w_C/sum_w, lim_N ~ w_N/sum_w, lim_P ~ w_P/sum_w,# normalized for plot
        # α_LT, α_RT by get_revenue_eq_X
        #D(α_L) ~ dα_L, dα_L ~ (α_LT - α_L)*(τ + abs(syn_B)/B),
        α_L ~ 1 - (α_R + α_LP + α_RP),
        D(α_R) ~ dα_R, dα_R ~ (α_RT - α_R)*(τ + abs(syn_B)/B),
        D(α_LP) ~ dα_LP, dα_LP ~ (α_LPT - α_LP)*(τ + abs(syn_B)/B),
        D(α_RP) ~ dα_RP, dα_RP ~ (α_RPT - α_RP)*(τ + abs(syn_B)/B),
        ]
    extend(ODESystem(vcat(eqs,eqs_rev), t, vcat(sts, sts_rev), ps; name), sP)
end

sesam3(args...;kwargs...) = sesam3CNP(args...;kwargs...)

"""
R pool is a mixture of microbial turnover and enzymes.
Here the stable C/N ratio is computed based on given parameterization.
"""
function calculate_β_PR_sesam3(p,s)
    w_REnz = p[s.a_E]*(1-p[s.κ_E]) # a_E (1-κ_E) B
    w_RB = p[s.ϵ_tvr]*p[s.τ]       # τ ϵ_tvr B
    β_PR = 1/((w_REnz/p[s.β_PEnz] + w_RB/p[s.β_PB])/(w_REnz + w_RB))
end








