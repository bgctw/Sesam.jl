function sesam3C(; name, k_N = 60.0, δ = 40.0, max_w = 12)
    @parameters t
    D = Differential(t)

    ps = @parameters(ϵ, ϵ_tvr, κ_E, a_E, m, τ,
        k_L, k_R, k_mN_L, k_mN_R, k_N=k_N,
        ρ_CBtvr=0.0,)

    sts = @variables (begin
        B(t), L(t), R(t), cumresp(t),
        dB(t), dL(t), dR(t), r_tot(t),
        syn_Enz(t), tvr_Enz(t), r_M(t), tvr_B(t),
        E_L(t), E_R(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t),
        r_B(t), r_GEnz(t), r_O(t),
        tvr_B0(t), resorp_C(t),
        w_C(t), lim_C(t),
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
        D(R) ~ dR, dR ~ -dec_R + ϵ_tvr * tvr_B + (1 - κ_E) * tvr_Enz,
        D(α_R) ~ dα_R,
        syn_Enz ~ a_E * B, tvr_Enz ~ syn_Enz,
        r_M ~ m * B,
        # element resorption changes stoichiometry: β_EB -> β_EBtvr
        tvr_B ~ tvr_B0 - resorp_C,
        tvr_B0 ~ τ * B, # turnover before resorption 
        resorp_C ~ ρ_CBtvr * tvr_B0,
        dec_LPot ~ k_L * L,
        dec_L ~ dec_LPot * (α_L * syn_Enz) / (k_mN_L + α_L * syn_Enz),
        dec_RPot ~ k_R * R,
        dec_R ~ dec_RPot * (α_R * syn_Enz) / (k_mN_R + α_R * syn_Enz),
        u_C ~ dec_L + dec_R + κ_E * syn_Enz,
        C_synBCt ~ u_C - syn_Enz / ϵ - r_M,
        #C_synBC ~ IfElse.ifelse(C_synBCt > 0.0, ϵ*C_synBCt, C_synBCt), 
        C_synBC ~ C_synBCt - ((C_synBCt > 0.0) * (1 - ϵ) * C_synBCt),
        #syn_B ~ C_synBC # define in all-elements model
        r_O ~ u_C - (syn_Enz / ϵ + syn_B + r_G + r_M),
        #D(cumresp) ~ r_B + r_tvr,
        r_tvr ~ (1 - ϵ_tvr) * tvr_B,
        r_B ~ r_GEnz + r_G + r_M + r_O,
        r_GEnz ~ (1 - ϵ) / ϵ * syn_Enz,
        #r_G ~ IfElse.ifelse(syn_B.val > 0.0, (1-ϵ)/ϵ * syn_B, 0.0), 
        r_G ~ (syn_B.val > 0.0) * (1 - ϵ) / ϵ * syn_B,
        r_tot ~ r_B + r_tvr,
        E_L ~ (α_L * syn_Enz) / k_N,
        E_R ~ (α_R * syn_Enz) / k_N,
        w_C ~ exp(min(max_w, -δ / tvr_B * (C_synBC - syn_B))), lim_C ~ w_C / sum_w,
    ]
    ODESystem(eqs, t, sts, ps; name)
end

function sesam3N(; name, δ = 40.0, max_w = 12, sC = sesam3C(; name = :sC, δ, max_w))
    @parameters t
    D = Differential(t)

    sts = @variables (begin
        L_N(t), R_N(t), I_N(t),
        dL_N(t), dR_N(t), dI_N(t),
        Φ_N(t), Φ_Nu(t), Φ_NB(t), Φ_Ntvr(t),
        u_PlantN(t), u_NOM(t), u_immNPot(t),
        u_NPot(t), N_synBN(t), C_synBN(t), M_ImbN(t),
        β_NL(t), β_NR(t),
        leach_N(t),
        resorp_N(t), β_NBtvr(t),
        p_uNmic(t), ν_TN(t),
        # need to be specified by coupled system:
        β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t),
        SOC(t), SON(t), β_NSOM(t), BtoSOC(t),
        w_N(t), lim_N(t)
    end)
    ps = @parameters β_NEnz β_NB l_N ν_N i_BN ρ_NBtvr=0.0

    @unpack L, R, B, dec_L, dec_R, i_L, ϵ_tvr, tvr_B, tvr_B0, syn_B, syn_Enz, tvr_Enz = sC
    @unpack r_tvr, κ_E, sum_w = sC

    eqs = [
        β_NL ~ L / L_N, β_NR ~ R / R_N,
        D(L_N) ~ dL_N, dL_N ~ -dec_L / β_NL + i_L / β_Ni,
        D(R_N) ~ dR_N,
        dR_N ~ -dec_R / β_NR + ϵ_tvr * tvr_B / β_NBtvr + (1 - κ_E) * tvr_Enz / β_NEnz,
        D(I_N) ~ dI_N,
        u_PlantN ~ min(u_PlantNmax, k_PlantN * I_N),
        dI_N ~ i_IN - u_PlantN - leach_N + Φ_N,
        leach_N ~ l_N * I_N,
        Φ_N ~ Φ_Nu + Φ_NB + Φ_Ntvr,
        Φ_Ntvr ~ r_tvr / β_NBtvr,
        Φ_Nu ~ (1 - ν_N) * u_NOM,
        u_NPot ~ ν_N * u_NOM + u_immNPot + resorp_N,
        u_NOM ~ dec_L / β_NL + dec_R / β_NR + κ_E * tvr_Enz / β_NEnz,
        u_immNPot ~ i_BN * I_N,
        N_synBN ~ u_NPot - syn_Enz / β_NEnz,
        C_synBN ~ β_NB * N_synBN,
        M_ImbN ~ u_NPot - (syn_B / β_NB + syn_Enz / β_NEnz),
        Φ_NB ~ M_ImbN - u_immNPot,
        # make sure ρ_NBtvr >= 0
        resorp_N ~ ρ_NBtvr * tvr_B0 / β_NB,
        β_NBtvr ~ tvr_B / (tvr_B0 / β_NB - resorp_N),
        p_uNmic ~ u_immNPot / (u_immNPot + u_PlantN),
        ν_TN ~ ν_N + (1 - ν_N) * p_uNmic,
        w_N ~ exp(min(max_w, -δ / tvr_B * (C_synBN - syn_B))),
        lim_N ~ w_N / sum_w,
        # observables for diagnostic output
        SOC ~ L + R + B,
        SON ~ L_N + R_N + B / β_NB,
        β_NSOM ~ SOC / SON,
        BtoSOC ~ B / SOC,
    ]
    extend(ODESystem(eqs, t, sts, ps; name), sC)
end

function get_dα_eq_sesam3CN(sN)
    @parameters t
    @unpack α_L, α_R, dec_L, dec_R, β_NL, β_NR, β_NB, β_NEnz, syn_Enz, τ, syn_B, B, ϵ = sN
    @unpack u_immNPot, u_PlantN, ν_N, ν_TN = sN
    @unpack ω_Enz, ω_L, ω_R = sN
    sts = @variables (begin
        α_LT(t), α_RT(t),
        invest_L(t), invest_R(t), return_L(t), return_R(t), revenue_L(t), revenue_R(t),
        #invest_Ln(t), invest_Rn(t), return_Ln(t), return_Rn(t), 
        revenue_sum(t)
    end)
    # defined for system
    @unpack dα_R = sN
    eqs = [
        invest_L ~ α_L * syn_Enz * ω_Enz,
        invest_R ~ α_R * syn_Enz * ω_Enz,
        return_L ~ dec_L * ω_L,
        return_R ~ dec_R * ω_R,
        revenue_L ~ return_L / invest_L,
        revenue_R ~ return_R / invest_R,
        revenue_sum ~ revenue_L + revenue_R,
        α_LT ~ revenue_L / revenue_sum,
        α_RT ~ revenue_R / revenue_sum,
        dα_R ~ (α_RT - α_R) * (τ + abs(syn_B) / B),
        # auxiliary for plotting
        # invest_Ln ~ invest_L/(invest_L + invest_R),
        # invest_Rn ~ invest_R/(invest_L + invest_R),
        # return_Ln ~ return_L/(return_L + return_R),
        # return_Rn ~ return_R/(return_L + return_R),
    ]
    (; eqs, sts)
end

function get_dα_eq_sesam3CN_deriv(sN)
    @parameters t
    @unpack α_L, α_R, dec_LPot, dec_RPot, β_NL, β_NR, β_NB, β_NEnz, syn_Enz, ϵ = sN
    @unpack ω_Enz, ω_L, ω_R = sN
    @unpack τ, syn_B, B, k_mN_L, k_mN_R = sN
    @unpack u_immNPot, u_PlantN, ν_N, ν_TN = sN
    sts = @variables (begin
        d_L(t), d_R(t),
        du_L(t), du_R(t),
        mdu(t)
    end)
    @unpack dα_R = sN
    eqs = [
        d_L ~ dec_LPot * ω_L,
        d_R ~ dec_RPot * ω_R,
        du_L ~ syn_Enz * k_mN_L * d_L / (k_mN_L + α_L * syn_Enz)^2,
        du_R ~ syn_Enz * k_mN_R * d_R / (k_mN_R + α_R * syn_Enz)^2,
        mdu ~ compute_mean_du2(du_L, α_L, du_R, α_R),
        dα_R ~ (τ + abs(syn_B) / B) * max((du_R - mdu) / mdu, -α_R),
    ]
    (; eqs, sts)
end

"""  
    compute_mean_du2(du1, α1, du2, α2; sum_α_others=0)

Compute mean across derivates, du, with checking excluding enzymes for 2 enzymes
Given derivates du1, du2 and allocations α1 and α2.
(Wutzler22 Appendix D1)

sum_α_others: allocation to enzymes that should be excluded from the mean

Checks all the combinations of relative change of du being smaller than -α_Z
"""
function compute_mean_du2(du1, α1, du2, α2; sum_α_others = 0)
    mdu_raw = (du1 + du2) / 2
    mdu2 = 2 * mdu_raw / (2 + sum_α_others) # both enzymes, adjusted for otherwise excluded enzymes
    mdu = IfElse.ifelse((du2 - mdu2) <= -α2 * mdu2, du1 / (1 + α2 + sum_α_others), # only first enzyme
        IfElse.ifelse((du1 - mdu2) <= -α1 * mdu2, du2 / (1 + α1 + sum_α_others), # only second enzyme
            mdu2))
    mdu
end

"""  
compute elemental-limitation weighting factor omega_Z
"""
function compute_elemental_weightfactor(lim_E, betaZ, β_B,
    ν_TZ = Fill(1.0, length(lim_E)))
    # omega_Z will be multiplied by the mineralization flux
    # for depolymerizing enzymes this is carbon-flux and beta_Z denote C:E ratios
    # for biomineralizing enzymes this is already P-flux and betaZ=1 and only limP
    #   needs to be provided
    wE = lim_E .* ν_TZ .* β_B ./ betaZ
    sum(wE)
end

function sesam3CN(; name, δ = 40.0, max_w = 12,
    use_seam_revenue = false, use_proportional_revenue = false,
    sN = sesam3N(; name = :sN, δ, max_w))
    @parameters t
    D = Differential(t)
    @unpack α_L, α_R, syn_B, B, C_synBC, N_synBN, tvr_B, τ, ϵ = sN
    @unpack β_NB, β_NEnz, β_NL, β_NR = sN
    @unpack lim_C, lim_N, ω_Enz, ω_L, ω_R = sN
    @unpack u_immNPot, u_PlantN, ν_TN = sN
    @unpack sum_w, w_C, w_N, C_synBN, dα_R = sN
    ps = @parameters δ = δ
    eqs_rev, sts_rev = use_seam_revenue ?
                       get_dα_eq_seam(sN) :
                       use_proportional_revenue ?
                       get_dα_eq_sesam3CN(sN) : get_dα_eq_sesam3CN_deriv(sN)
    @variables α_LT(t) α_RT(t)
    # two elements in weighting Enz, L and R carbon fluxes
    lim_E = SA[lim_C, lim_N]
    β_B = SA[1.0, β_NB]
    ν_TZ = SA[ϵ, ν_TN]
    eqs = [
        syn_B ~ min(C_synBC, C_synBN),
        # need minimum, otherwise danger of Inf and nan -> instability
        sum_w ~ w_C + w_N,
        ω_Enz ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NEnz], β_B),
        ω_L ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NL], β_B, ν_TZ),
        ω_R ~ compute_elemental_weightfactor(lim_E, SA[1.0, β_NR], β_B, ν_TZ),
        # α_LT, α_RT, dα_R by get_dα_eq_X
        #D(α_L) ~ dα_L, dα_L ~ (α_LT - α_L)*(τ + abs(syn_B)/B),
        α_L ~ 1 - α_R,
    ]
    extend(ODESystem(vcat(eqs, eqs_rev), t, vcat(sts_rev), ps; name), sN)
end

# sesam3 now refers to sesam3CNP
#sesam3(args...;kwargs...) = sesam3CN(args...;kwargs...)

"""
R pool is a mixture of microbial turnover and enzymes.
Here the stable C/N ratio is computed based on given parameterization.
"""
function calculate_β_NR_sesam3(p, s)
    w_REnz = p[s.a_E] * (1 - p[s.κ_E]) # a_E (1-κ_E) B
    w_RBtvr = p[s.ϵ_tvr] * p[s.τ]       # τ ϵ_tvr B
    # β_NBtvr may differ from β_NB because of elemental resorption
    # if ρ is not in parameter dictionary, its default is zero
    ρ_CBtvr = (hasproperty(s, :ρ_CBtvr) && haskey(p, s.ρ_CBtvr)) ? p[s.ρ_CBtvr] : 0.0
    ρ_NBtvr = (hasproperty(s, :ρ_NBtvr) && haskey(p, s.ρ_NBtvr)) ? p[s.ρ_NBtvr] : 0.0
    β_NBtvr = p[s.β_NB] * (1 - ρ_CBtvr) / (1 - ρ_NBtvr)
    β_NR = 1 / ((w_REnz / p[s.β_NEnz] + w_RBtvr / β_NBtvr) / (w_REnz + w_RBtvr))
end

function calculate_β_NR_sesam3(p_sym)
    pREnz = calculate_shareEnzTvr(p_sym)
    β_NBtvr = calculate_β_NBtvr(p_sym)
    β_NR = 1 / (pREnz / p_sym[:s₊β_NEnz] + (1 - pREnz) / β_NBtvr)
end

"Calculate the share of enzymes tvr contributing to total tvr in addition to biomass"
function calculate_shareEnzTvr(p_sym)
    w_REnz = p_sym[:s₊a_E] * (1 - p_sym[:s₊κ_E])    # a_E (1-κ_E) B
    w_RBtvr = p_sym[:s₊ϵ_tvr] * p_sym[:s₊τ]       # τ ϵ_tvr B
    w_REnz / (w_REnz + w_RBtvr)
end

"""
Calculate the C/N ratio of microbial biomass turnover.
It may differ from β_NB because of elemental resorption.
"""
function calculate_β_NBtvr(p_sym)
    # if ρ is not in parameter dictionary, its default is zero
    ρ_CBtvr = (haskey(p_sym, :s₊ρ_CBtvr)) ? p_sym[:s₊ρ_CBtvr] : 0.0
    ρ_NBtvr = (haskey(p_sym, :s₊ρ_NBtvr)) ? p_sym[:s₊ρ_NBtvr] : 0.0
    p_sym[:s₊β_NB] * (1 - ρ_CBtvr) / (1 - ρ_NBtvr)
end

"""
Given the C/N ratio of total SOM, litter input and input to R: tvr of B and Enz
calculate the proportion of residue pool R to total SOM
"""
function calculate_propR(p_sym, β_NSOM; β_L = p_sym[:pl₊β_Ni0])
    β_NBtvr = calculate_β_NBtvr(p_sym)
    β_NEnz = p_sym[:s₊β_NEnz]
    pREnz = calculate_shareEnzTvr(p_sym)
    # y is the inverse of C/N ratio beta: N:C
    yR = (1 - pREnz) / β_NBtvr + pREnz / β_NEnz
    yS = 1 / β_NSOM
    yL = 1 / β_L
    (yS - yL) / (yR - yL)
end
