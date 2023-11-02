using ModelingToolkit: AbstractSystem, get_eqs, get_states, get_ps, get_observed,
    get_continuous_events, get_defaults, get_systems

# moved override to MTKHelpers.override_system

function sesam_const_dLRP(dL0, dR0, dP0; name, kwargs...)
    s = sesam3(; name = :s, kwargs...)
    @unpack d_L, d_R, d_P, B = s
    @parameters t
    ps = @parameters d_L0=dL0 d_R0=dR0 d_P0=dP0
    D = Differential(t)
    eqs = [
        d_L ~ d_L0,
        d_R ~ d_R0,
        d_P ~ d_P0,
        D(B) ~ 0,
    ]
    sys_ext = override_system(eqs, s; name, ps) # extend does not overwrite
end

function sesam_const_dLRP_relative(dL0, dR0, dP0; name, kwargs...)
    s = sesam3(; name = :s, use_proportional_revenue = true, kwargs...)
    @unpack dec_LPot, dec_RPot, dec_LPPot, dec_RPPot, dec_PPlant, B = s
    @unpack w_C, w_P, w_N, lim_C, lim_P, lim_N, β_PL, β_PR, β_PB, β_NB = s
    @unpack ω_Enz, ω_L, ω_R, ω_P = s
    @unpack ν_TN, ν_TP = s
    ps = @parameters d_L0=dL0 d_R0=dR0 d_P0=dP0
    @parameters t
    D = Differential(t)
    # there is no clear correspondence to sesam_LRP_deriv.Rmd
    # has to specify something before mutlipyling with enzmye-limitation-factor
    # but after weighting elemental limitations
    # need to rework dec_LPot so that dec_L matches dL0
    eqs = [
        dec_LPPot ~ d_P0 / 2,
        dec_RPPot ~ d_P0 / 2,
        # subtract P flux - because its multiplied by weight > 1 including P flux
        # dec_LPot ~ d_L0 - dec_LPPot*w_P,#*β_PB/β_PL,  # creates implicit eq.
        # dec_RPot ~ d_R0 - dec_RPPot*w_P*β_PB/β_PR,
        # dec_LPot ~ d_L0/ω_L, # creates implicit eq.
        # dec_RPot ~ d_R0/ω_R,
        dec_LPot ~ d_L0,
        dec_RPot ~ d_R0,
        dec_PPlant ~ 0,
        D(B) ~ 0,
        # # w_C ~ 1,
        # # w_P ~ 1,
        # # w_N ~ 0,
        # # lim_C ~ 1,
        # # lim_P ~ 1,
        # # w_N ~ 0,
        ω_Enz ~ 1, ω_L ~ 1, ω_R ~ 1, ω_P ~ 1,
        ν_TN ~ 1, ν_TP ~ 1,
    ]
    sys_ext = override_system(eqs, s; name, ps) # extend does not overwrite
end

function sesam_fixed_substrates(s; name = nameof(s), kwargs...)
    # keep all the substrate pools fixed and simulate only biomass B and α
    @unpack L, R, L_N, R_N, I_N, L_P, R_P, I_P = s
    @parameters t
    D = Differential(t)
    #dR <- dL <- dRN <- dLN <- dRP <- dLP <- dIN <- dIP <- 0
    eqs = [
        D(L) ~ 0,
        D(R) ~ 0,
        D(L_N) ~ 0,
        D(R_N) ~ 0,
        D(I_N) ~ 0,
        D(L_P) ~ 0,
        D(R_P) ~ 0,
        D(I_P) ~ 0,
    ]
    sys_ext = override_system(eqs, s; name)
end
