"""
model i_L and β_Ni as as step function increased by factor fac_int at t in [t1,t2).
Superimpose a seasonal pattern on the litter input.
"""
function plant_face_fluct(; name, t1 = 0.0, t2 = 100.0, fac_inc = 1.2,
    #    autumn_start=10/12, autumn_end = 11/12, 
    autumn_start = 8.5 / 12, autumn_end = 11.5 / 12,
    share_autumn = 0.5, k_Lagr = 12 / 2, k_PlantN0v = 100, k_PlantP0v = 100)
    @parameters t
    D = Differential(t)
    @parameters i_L0 i_IN0 β_Ni0 i_IP0 β_Pi0 t1=t1 t2=t2 fac_inc=fac_inc
    @parameters u_PlantPmax0 k_PlantP0=k_PlantP0v s_EP0
    @parameters u_PlantNmax0 k_PlantN0=k_PlantN0v
    @parameters k_Lagr=k_Lagr share_autumn=share_autumn
    # in the fluctuation analysis the integrated
    # model gives slightly higher litter input
    # ad-hoc correct for providiing slightly lower 
    # amount in autumn

    # Distribution of litterfall within year (between 0 and 1)
    # avoid dependencies to DistributionFits and Optim - hardcode    
    #dist_flat = fit_mode_flat(LogitNormal, 0.3; peakedness = 3)
    dist_flat = LogitNormal{Float64}(-0.7545014913579478, 0.48165436006864837)
    d_lit_agr = autumn_start + (autumn_end - autumn_start) * dist_flat
                
    @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        i_L_anomaly(t), i_L_annual(t), β_Ni_annual(t),
        Lagr(t), i_Lagr(t), dec_Lagr(t), d_Lagr(t),
        u_PlantNmax(t), k_PlantN(t), inc_period(t),
        β_Pi_annual(t), β_Pi(t), i_IP(t),
        u_PlantPmax(t), k_PlantP(t),
        s_EP(t)
    end)
    smooth_dt = 1 / 12 / 2
    eqs = [
        D(Lagr) ~ d_Lagr, d_Lagr ~ i_Lagr - dec_Lagr,
        dec_Lagr ~ k_Lagr * Lagr,
        i_Lagr ~ share_autumn * i_L_annual * i_L_anomaly,
        i_L ~ (1 - share_autumn) * i_L_annual + dec_Lagr,
        # i_L_annual ~ IfElse.ifelse((t1 <= t) & (t < t2), fac_inc*i_L0, i_L0), 
        # β_Ni_annual ~ IfElse.ifelse((t1 <= t) & (t < t2), fac_inc*β_Ni0, β_Ni0),
        inc_period ~ smoothstep_sesam(t, t1, smooth_dt) *
                     (1 - smoothstep_sesam(t, t2, smooth_dt)),
        i_L_annual ~ i_L0 * (1 + (fac_inc - 1) * inc_period),
        β_Ni_annual ~ β_Ni0 * (1 + (fac_inc - 1) * inc_period),
        i_L_anomaly ~ get_iL_anomaly(t, d_lit_agr),
        β_Ni ~ β_Ni_annual, # does not change with season
        i_IN ~ i_IN0,
        u_PlantNmax ~ u_PlantNmax0,
        k_PlantN ~ k_PlantN0,
        #
        β_Pi_annual ~ β_Pi0 * (1 + (fac_inc - 1) * inc_period),
        β_Pi ~ β_Pi_annual,
        i_IP ~ i_IP0,
        u_PlantPmax ~ u_PlantPmax0,
        k_PlantP ~ k_PlantP0,
        s_EP ~ s_EP0,
    ]
    defaults = Pair{Num, Any}[
        # rate not constrained
        k_PlantN0 => k_PlantN0v,
        # take as much N and P as provide with litter
        u_PlantNmax0 => i_L0 / β_Ni0,
        u_PlantPmax0 => i_L0 / β_Pi0]
    continuous_events = vcat([
        t ~ t1,
        t ~ t2,
    ])
    ODESystem(eqs, t; name, defaults, continuous_events)
end

# now defined in MTKHelpers - not found in generated function
smoothstep_sesam(args...) = smoothstep(args...)
@register_symbolic smoothstep_sesam(t, t_step::Number, dt::Number)

get_iL_anomaly(t, dist) = pdf(dist, t - floor(t))
@register_symbolic get_iL_anomaly(t, dist::Distribution)

function f_interpolate9(t, sol, numref)
    sol(t, idxs = numref.x)
end
@register_symbolic f_interpolate9(t, sol::ODESolution, numref::Base.RefValue{Num})

function plant_face_fluct_fake(; name, sys, sol, t1, t2)
    @parameters t
    D = Differential(t)
    ps = @parameters t1=t1 t2=t2 # for continuous events
    sts = @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t),
        i_L_annual(t)
    end)
    eqs = [
        i_L ~ f_interpolate9(t, sol, Ref(sys.i_L)),
        β_Ni ~ f_interpolate9(t, sol, Ref(sys.β_Ni)),
        i_IN ~ f_interpolate9(t, sol, Ref(sys.i_IN)),
        u_PlantNmax ~ f_interpolate9(t, sol, Ref(sys.u_PlantNmax)),
        k_PlantN ~ f_interpolate9(t, sol, Ref(sys.k_PlantN)),
        i_L_annual ~ f_interpolate9(t, sol, Ref(sys.i_L_annual)),
    ]
    # continuous_events = vcat(
    #     [t-yr ~ autumn_start*1.05 for yr in -500:200],
    #     [t-yr ~ autumn_end*0.995 for yr in -500:200],
    #     [
    #     t ~ t1,
    #     t ~ t2,
    #     ])
    #@unpack autumn_start, autumn_end = sys # in case those become parameters
    continuous_events = first(sys.continuous_events.value).eqs
    ODESystem(eqs, t, sts, ps; name, continuous_events)
end
