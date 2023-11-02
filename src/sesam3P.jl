function sesam3P(; name, δ = 40.0, max_w = 12, sN = sesam3N(; name = :sN, δ, max_w))
    @parameters t
    D = Differential(t)

    sts = @variables (begin
        L_P(t), R_P(t), I_P(t),
        dL_P(t), dR_P(t), dI_P(t),
        i_LP(t), tvr_P(t), tvr_Enz_P(t),
        Φ_P(t), Φ_Pu(t), Φ_PB(t), Φ_Ptvr(t),
        u_PlantP(t), u_POM(t), u_immPPot(t),
        u_PPot(t), P_synBP(t), C_synBP(t), M_ImbP(t),
        β_PL(t), β_PR(t),
        leach_P(t),
        resorp_P(t), β_PBtvr(t),
        α_P(t), lim_enz_P(t), dec_LPPot(t), dec_RPPot(t), dec_PPot(t),
        dec_LP_P(t), dec_RP_P(t), dec_PPlant(t),
        lim_LP(t), lim_RP(t),
        p_uPmic(t), ν_TP(t),
        # need to be specified by coupled system:
        β_Pi(t), i_IP(t),
        u_PlantPmax(t), k_PlantP(t),
        s_EP(t), pL_sEP(t), # synthesis of E_LP adn E_RP enzymes by plants
        SOP(t), β_PSOM(t), β_NPSOM(t), β_NPB(t), β_NPi(t), β_NPL(t), β_NPR(t),
        lim_P(t), ω_P(t), ν_TP(t),
        dα_P(t), w_P(t)
    end)
    ps = @parameters(β_PEnz, β_PB, l_P, ν_P, i_BN, i_BP, k_LP, k_RP, k_mN_P, β_Pm,
        ρ_PBtvr=0.0)

    @unpack L, R, B, SOC, SON, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, syn_B, tvr_B0 = sN
    @unpack syn_Enz, tvr_Enz, r_tvr, κ_E = sN
    @unpack β_NB, β_Ni, L_N, R_N = sN
    @unpack k_mN_L, k_mN_R = sN
    @unpack sum_w = sN

    eqs = [
        β_PL ~ L / L_P, β_PR ~ R / R_P,
        D(L_P) ~ dL_P, dL_P ~ -dec_L / β_PL - dec_LP_P + i_LP,
        D(R_P) ~ dR_P, dR_P ~ -dec_R / β_PR - dec_RP_P + tvr_P + tvr_Enz_P,
        D(I_P) ~ dI_P,
        D(α_P) ~ dα_P,
        i_LP ~ i_L / β_Pi, tvr_P ~ ϵ_tvr * tvr_B / β_PBtvr,
        tvr_Enz_P ~ (1 - κ_E) * tvr_Enz / β_PEnz,
        # potential rate of biomineralization decreases with  increasing C:P ratio β_PS
        lim_LP ~ 1 / (1 + β_PL / β_Pm),
        lim_RP ~ 1 / (1 + β_PR / β_Pm),
        dec_LPPot ~ (k_LP * L_P * lim_LP),
        dec_RPPot ~ (k_RP * R_P * lim_RP),
        dec_PPot ~ dec_LPPot + dec_RPPot,
        lim_enz_P ~ (s_EP + α_P * syn_Enz) / (k_mN_P + s_EP + α_P * syn_Enz),
        dec_LP_P ~ dec_LPPot * lim_enz_P,
        dec_RP_P ~ dec_RPPot * lim_enz_P,
        # although its not linear, estimate biomineralization by plant enzyme production only
        dec_PPlant ~ dec_PPot * (s_EP) / (k_mN_L + s_EP),
        u_PlantP ~ min(u_PlantPmax, k_PlantP * I_P),
        dI_P ~ i_IP - u_PlantP - leach_P + dec_LP_P + dec_RP_P + Φ_P,
        leach_P ~ l_P * I_P,
        Φ_P ~ Φ_Pu + Φ_PB + Φ_Ptvr,
        Φ_Ptvr ~ r_tvr / β_PBtvr,
        Φ_Pu ~ (1 - ν_P) * u_POM,
        # dec_RP_P is biomineralization - cleaved P is all mineralized
        # u_POM is from depolymerizing - only a fraction ν_P is mineralized
        u_PPot ~ ν_P * u_POM + u_immPPot + resorp_P,
        u_POM ~ dec_L / β_PL + dec_R / β_PR + κ_E * tvr_Enz / β_PEnz,
        u_immPPot ~ i_BP * I_P,
        P_synBP ~ u_PPot - syn_Enz / β_PEnz,
        C_synBP ~ β_PB * P_synBP,
        M_ImbP ~ u_PPot - (syn_B / β_PB + syn_Enz / β_PEnz),
        Φ_PB ~ M_ImbP - u_immPPot,
        # make sure ρ_PBtvr >= 0
        resorp_P ~ ρ_PBtvr * tvr_B0 / β_PB,
        β_PBtvr ~ tvr_B / (tvr_B0 / β_PB - resorp_P),
        # loss by indirect pathway
        ν_TP ~ ν_P + (1 - ν_P) * p_uPmic,
        p_uPmic ~ u_immPPot / (u_immPPot + u_PlantP),
        w_P ~ exp(min(max_w, -δ / tvr_B * (C_synBP - syn_B))), lim_P ~ w_P / sum_w,
        # observables for diagnostic output
        SOP ~ L_P + R_P + B / β_PB,
        β_PSOM ~ SOC / SOP,
        β_NPSOM ~ SON / SOP,
        β_NPB ~ β_PB / β_NB, # N:P of microbial biomass
        β_NPi ~ β_Pi / β_Ni, # N:P of litter inputs
        β_NPL ~ L_N / L_P,
        β_NPR ~ R_N / R_P,
    ]
    extend(ODESystem(eqs, t, sts, ps; name), sN)
end

function get_dα_eq_sesam3CNP(sP)
    @parameters t
    @unpack α_L, α_R, dec_L, dec_R, β_NL, β_NR, β_NB, β_NEnz, syn_Enz, ϵ = sP
    @unpack β_PL, β_PR, β_PB, β_PEnz = sP
    @unpack lim_C, lim_N, lim_P, ω_Enz, ω_L, ω_R, ω_P = sP
    @unpack α_P, dec_LP_P, dec_RP_P, dec_PPlant, syn_B, B, τ = sP
    @unpack u_immPPot, u_PlantP, u_immNPot, u_PlantN, ν_TN, ν_TP = sP
    sts = @variables (begin
        α_LT(t), α_RT(t),
        α_PT(t),
        syn_Enz_w(t),
        return_L(t), return_R(t), revenue_L(t), revenue_R(t),
        return_P(t), revenue_P(t),
        #invest_Ln(t), invest_Rn(t), return_Ln(t), return_Rn(t), 
        revenue_sum(t),
        τsyn(t), dα_R(t), dα_P(t)
    end)
    lim_E = SA[lim_C, lim_N, lim_P]
    β_B = SA[1.0, β_NB, β_PB]
    ν_TZ = SA[ϵ, ν_TN, ν_TP]
    eqs = [
        syn_Enz_w ~ syn_Enz * ω_Enz,
        # return_L ~ dec_L * (w_C*ϵ + w_N/β_NL*ν_TN + w_P/β_PL*ν_TP), 
        # return_R ~ dec_R * (w_C*ϵ + w_N/β_NR*ν_TN + w_P/β_PR*ν_TP), 
        # return_L ~ dec_L * (w_C + w_N/β_NL + w_P/β_PL), 
        # return_R ~ dec_R * (w_C + w_N/β_NR + w_P/β_PR), 
        return_L ~ dec_L * ω_L,
        return_R ~ dec_R * ω_R,
        # return of enzymes produced in addition to that of plants
        # only proportion of the biomineralization flux ends up in microbes: part it taken up by plant
        # TODO dec_LP_P is already in P units, no need to devide by β_P
        #return_P ~ (dec_LP_P + dec_RP_P - dec_PPlant) * p_uPmic * w_P, 
        return_P ~ (dec_LP_P + dec_RP_P - dec_PPlant) * ω_P,
        revenue_L ~ return_L / (α_L * syn_Enz_w),
        revenue_R ~ return_R / (α_R * syn_Enz_w),
        # avoid deviding by zero. Checking > 0 fails with Derivative need to check > small_number
        revenue_P ~ IfElse.ifelse(α_P * syn_Enz_w > 1e-16,
            return_P / (α_P * syn_Enz_w),
            0.0),
        revenue_sum ~ revenue_L + revenue_R + revenue_P,
        α_LT ~ revenue_L / revenue_sum,
        α_RT ~ revenue_R / revenue_sum,
        α_PT ~ revenue_P / revenue_sum,
        # auxilary for plotting
        # invest_Ln ~ invest_L/(invest_L + invest_R),
        # invest_Rn ~ invest_R/(invest_L + invest_R),
        # return_Ln ~ return_L/(return_L + return_R),
        # return_Rn ~ return_R/(return_L + return_R),
        τsyn ~ τ + abs(syn_B) / B,
        dα_R ~ (α_RT - α_R) * τsyn,
        dα_P ~ (α_PT - α_P) * τsyn,
    ]
    (; eqs, sts)
end

function get_dα_eq_sesam3CNP_deriv(sP)
    @parameters t
    @unpack α_L, α_R, dec_LPot, dec_RPot, β_NL, β_NR, β_NB, β_NEnz, syn_Enz, ϵ = sP
    @unpack β_PL, β_PR, β_PB, β_PEnz = sP
    @unpack lim_C, lim_N, lim_P, ω_Enz, ω_L, ω_R, ω_P = sP
    @unpack α_P, dec_PPot, k_mN_L, k_mN_R, k_mN_P, s_EP, syn_B, B, τ = sP
    @unpack u_immPPot, u_PlantP, u_immNPot, u_PlantN, ν_N, ν_P, ν_TP(t), ν_TN = sP
    @unpack dα_R, dα_P = sP
    sts = @variables (begin
        syn_Enz_w(t),
        d_L(t), d_R(t), d_P(t),
        du_L(t), du_R(t), du_P(t),
        mdu(t),
        τsyn(t)
        #dα_L(t),
    end)
    eqs = [
        syn_Enz_w ~ syn_Enz * ω_Enz,
        d_L ~ dec_LPot * ω_L,
        d_R ~ dec_RPot * ω_R,
        d_P ~ dec_PPot * ω_P,
        du_L ~ syn_Enz * k_mN_L * d_L / (k_mN_L + α_L * syn_Enz)^2,
        du_R ~ syn_Enz * k_mN_R * d_R / (k_mN_R + α_R * syn_Enz)^2,
        # TODO here only k_mN_R (assume k_mN_L same)
        du_P ~ syn_Enz * k_mN_R * d_P / (s_EP + k_mN_P + α_P * syn_Enz)^2,
        mdu ~ compute_mean_du3(du_L, α_L, du_R, α_R, du_P, α_P),
        τsyn ~ τ + abs(syn_B) / B,
        #dα_L ~ τsyn * max((du_L - mdu)/mdu, -α_L),
        dα_R ~ τsyn * max((du_R - mdu) / mdu, -α_R),
        dα_P ~ τsyn * max((du_P - mdu) / mdu, -α_P),
    ]
    (; eqs, sts)
end

"""  
    compute_mean_du3(du1, α1, du2, α2, du3, α3; sum_α_others=0)

Compute mean across derivates, du, with checking exluding enzymes for 2 enzymes
Given derivates du1, du2, du3 and allocations α1, α2, α3.
(Wutzler22 Appendix D1)

sum_α_others: allocation to enzymes that should be excluded from the mean

Checks all the combinations of relative change of du being smaller than -α_Z
"""
function compute_mean_du3(du1, α1, du2, α2, du3, α3; sum_α_others = 0)
    mdu_raw = (du1 + du2 + du3) / 3
    mdu3 = 3 * mdu_raw / (3 + sum_α_others) # adjust for otherwise excluded enzymes
    mdu = IfElse.ifelse((du3 - mdu3) <= -α3 * mdu3,
        # exclude enzyme 3 and check remaining two
        compute_mean_du2(du1, α1, du2, α2; sum_α_others = sum_α_others + α3),
        IfElse.ifelse((du2 - mdu3) <= -α2 * mdu3,
            # exclude enzyme 2 and check remaining two
            compute_mean_du2(du1, α1, du3, α3; sum_α_others = sum_α_others + α2),
            IfElse.ifelse((du1 - mdu3) <= -α1 * mdu3,
                # exclude enzyme 1 and check remaining two
                compute_mean_du2(du2, α2, du3, α3; sum_α_others = sum_α_others + α1),
                # exclude no enzyme
                mdu3)))
    mdu
end

function sesam3CNP(;
    name, δ = 40.0, max_w = 12, use_proportional_revenue = false,
    sP = sesam3P(; name = :sP, δ, max_w))
    @parameters t
    D = Differential(t)
    @unpack α_L, α_R, syn_B, B, C_synBC, β_NB, N_synBN, tvr_B = sP
    @unpack lim_C, lim_N, lim_P, ω_Enz, ω_L, ω_R, ω_P = sP
    @unpack ϵ, ν_TN, ν_TP = sP
    @unpack β_PB, P_synBP, α_P = sP
    @unpack β_NL, β_NR, β_NEnz, β_PL, β_PR, β_PEnz = sP
    @unpack C_synBN, C_synBP, dα_P, w_P = sP
    @unpack w_C, w_N, w_P = sP
    @unpack sum_w = sP
    ps = @parameters δ=δ max_w=max_w
    eqs_rev, sts_rev = use_proportional_revenue ?
                       get_dα_eq_sesam3CNP(sP) : get_dα_eq_sesam3CNP_deriv(sP)
    lim_E = SA[lim_C, lim_N, lim_P]
    β_B = SA[1.0, β_NB, β_PB]
    ν_TZ = SA[ϵ, ν_TN, ν_TP]
    eqs = [
        syn_B ~ min(C_synBC, C_synBN, C_synBP),
        sum_w ~ w_C + w_N + w_P,
        ω_Enz ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NEnz, β_PEnz], β_B),
        ω_L ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NL, β_PL], β_B, ν_TZ),
        ω_R ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NR, β_PR], β_B, ν_TZ),
        ω_P ~ compute_elemental_weightfactor(SA[lim_P], SA[1.0], SA[β_PB], SA[ν_TP]),
        α_L ~ 1 - (α_R + α_P),
    ]
    extend(ODESystem(vcat(eqs, eqs_rev), t, vcat(sts_rev), ps; name), sP)
end

sesam3(args...; kwargs...) = sesam3CNP(args...; kwargs...)

"""
R pool is a mixture of microbial turnover and enzymes.
Here the stable C/N ratio is computed based on given parameterization.
"""
function calculate_β_PR_sesam3(p, s)
    w_REnz = p[s.a_E] * (1 - p[s.κ_E]) # a_E (1-κ_E) B
    w_RBtvr = p[s.ϵ_tvr] * p[s.τ]       # τ ϵ_tvr B
    # if ρ is not in parameter dictionary, its default is zero
    ρ_CBtvr = haskey(p, s.ρ_CBtvr) ? p[s.ρ_CBtvr] : 0.0
    ρ_PBtvr = haskey(p, s.ρ_PBtvr) ? p[s.ρ_PBtvr] : 0.0
    β_PBtvr = p[s.β_PB] * (1 - ρ_CBtvr) / (1 - ρ_PBtvr)
    β_PR = 1 / ((w_REnz / p[s.β_PEnz] + w_RBtvr / β_PBtvr) / (w_REnz + w_RBtvr))
end

function calculate_β_PR_sesam3(p_sym)
    # ensure with test that this is the same as with Dictionaries of Symbol.Num
    # for performance recode rather than converting Nums to Symbols each time
    w_REnz = p_sym[:s₊a_E] * (1 - p_sym[:s₊κ_E]) # a_E (1-κ_E) B
    w_RBtvr = p_sym[:s₊ϵ_tvr] * p_sym[:s₊τ]       # τ ϵ_tvr B
    # if ρ is not in parameter dictionary, its default is zero
    ρ_CBtvr = haskey(p_sym, :s₊ρ_CBtvr) ? p_sym[:s₊ρ_CBtvr] : 0.0
    ρ_PBtvr = haskey(p_sym, :s₊ρ_PBtvr) ? p_sym[:s₊ρ_PBtvr] : 0.0
    β_PBtvr = p_sym[:s₊β_PB] * (1 - ρ_CBtvr) / (1 - ρ_PBtvr)
    β_PR = 1 / ((w_REnz / p_sym[:s₊β_PEnz] + w_RBtvr / β_PBtvr) / (w_REnz + w_RBtvr))
end

"""
    get_updated_Penz_pars(pDict, s)

Create a copy of parameters and set k_mN_P, k_LP, and k_RP 
to values of the L/R enzymes.
"""
function get_updated_Penz_pars(pDict, s)
    p = copy(pDict)
    p[s.k_mN_P] = p[s.k_mN_R]
    p[s.k_RP] = p[s.k_R]
    p[s.k_LP] = p[s.k_L]
    p
end
function get_updated_Penz_pars(p_sym)
    p = copy(p_sym)
    p[:s₊k_mN_P] = p[:s₊k_mN_R]
    p[:s₊k_RP] = p[:s₊k_R]
    p[:s₊k_LP] = p[:s₊k_L]
    p
end
