function sesam3C_protect(;name, k_N=60.0, δ=40.0, max_w=12)
    # extends sesam3C by assuming only a portion of residues accessible to enzymatic
    # degradation. That portion is in Langmuir-equilibrium with a protected pool
    @parameters t 
    D = Differential(t)

    ps = @parameters(
        ϵ, ϵ_tvr, κ_E, a_E,  m,  τ,  
        k_L,  k_R,  k_mN_L, k_mN_R, k_N=k_N,
        ρ_NBtvr = 0.0, # used in N component, only relevant for variants implementing microbial nutrient resorption
        K_eqR=1e-3, Q_max=10_000.0,
    ) 

    sts = @variables (begin
        B(t),  L(t),   R(t),  cumresp(t), Ra(t), # accessible portion of R
        dB(t), dL(t), dR(t), r_tot(t),
        syn_Enz(t), tvr_Enz(t), r_M(t), tvr_B(t), tvr_B0(t),
        E_L(t), E_R(t),
        dec_LPot(t), dec_L(t), dec_RPot(t),
        dec_R(t), u_C(t),
        C_synBCt(t), C_synBC(t), r_G(t),
        r_tvr(t), 
        r_B(t), r_GEnz(t), r_O(t),
        w_C(t), lim_C(t),
        # new variables
        k_Rtot(t),
        # need to be defined by component across all elements
        α_L(t), α_R(t), 
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
        D(α_R) ~ dα_R, 
        syn_Enz ~ a_E*B, tvr_Enz ~ syn_Enz,
        r_M ~ m*B,
        tvr_B ~ τ*B,
        tvr_B0 ~ tvr_B, # used in N component, differs from tvr_B in microbila nutrient resorption variants
        dec_LPot ~ k_L * L,
        dec_L ~ dec_LPot*(α_L * syn_Enz)/(k_mN_L + α_L*syn_Enz),
        # the next two lines are different from sesam3
        # note the new observable Ra, k_Rtot and paramters Q_max, K_eqR
        dec_RPot ~ k_R * Ra,
        Ra ~ compute_accessible_langmuir(R, Q_max, K_eqR),
        k_Rtot ~ k_R * Ra/R,
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
        w_C ~ exp(min(max_w, -δ/tvr_B*(C_synBC - syn_B))), lim_C ~ w_C/sum_w, 
        ]
    ODESystem(eqs, t, sts, ps; name)    
end

function compute_accessible_langmuir(R, Q_max, K_eqR) 
    # positive solution of a quadratic equation 
    # see inst/pluto/R_langmuir.jl
    # provide Q_max and K_eqR as Ref-values to not interfere with symbolics when registering
    phalv = (Q_max - R + 1/K_eqR)/2
	q = -R/K_eqR
	-phalv + sqrt(phalv * phalv - q)    
end
@register_symbolic compute_accessible_langmuir(t, Q_max, K_eqR)

sesam3_protect(args...; δ=40.0, max_w=12, kwargs...) = sesam3CN(
    args...; 
    sN = sesam3N(name=:sN, sC = sesam3C_protect(;name=:sC, δ, max_w)),
    kwargs...
)






