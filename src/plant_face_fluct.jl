"""
model i_L and β_Ni as as step function increased by factor fac_int at t in [t1,t2).
Superimpose a seasonal pattern on the litter input.
"""
function plant_face_fluct(;name, t1=0.0, t2=100.0, fac_inc=1.2, 
    autumn_start=10/12, autumn_end = 11/12, share_autumn=0.5, k_Lagr=12/2, k_PlantN0v = 100, biasfac=1.0)
    @parameters t 
    D = Differential(t)
    @parameters i_L0  i_IN0  β_Ni0 t1=t1 t2=t2 fac_inc=fac_inc
    @parameters u_PlantNmax0  k_PlantN0=k_PlantN0v 
    @parameters k_Lagr = k_Lagr share_autumn = share_autumn
    # in the fluctuation analysis the integrated
    # model gives slightly higher litter input
    # ad-hoc correct for providiing slighly lower 
    # amount in autumn
    @parameters biasfac=biasfac
    

    @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        i_L_anomaly(t), i_L_annual(t), β_Ni_annual(t),
        Lagr(t), i_Lagr(t), dec_Lagr(t), d_Lagr(t),
        u_PlantNmax(t), k_PlantN(t), biasfac_agr(t)
    end) 
    eqs = [
        D(Lagr) ~ d_Lagr, d_Lagr ~ i_Lagr - dec_Lagr,
        dec_Lagr ~ k_Lagr*Lagr,
        # put entire bias to agr in autumn
        biasfac_agr ~ 1 + (biasfac - 1)/share_autumn,
        i_Lagr ~ share_autumn * i_L_annual * biasfac_agr * i_L_anomaly,
        i_L ~ (1-share_autumn) * i_L_annual + dec_Lagr,
        i_L_annual ~ IfElse.ifelse((t1 <= t) & (t < t2), fac_inc*i_L0, i_L0), 
        β_Ni_annual ~ IfElse.ifelse((t1 <= t) & (t < t2), fac_inc*β_Ni0, β_Ni0),
        i_L_anomaly ~ get_iL_anomaly(t, autumn_start, autumn_end),
        β_Ni ~ β_Ni_annual, # does not change with season
        i_IN ~ i_IN0, 
        u_PlantNmax ~ u_PlantNmax0,
        k_PlantN ~ k_PlantN0,
    ]
    defaults=Pair{Num, Any}[
        # rate not constrained
        k_PlantN0 => k_PlantN0v, 
        # take as much N as provide with litter
        u_PlantNmax0 => i_L0/β_Ni0, 
    ]
    continuous_events = vcat(
        [t-yr ~ autumn_start for yr in -500:200],
        [t-yr ~ autumn_end for yr in -500:200],
        [
        t ~ t1,
        t ~ t2,
        ])
    ODESystem(eqs,t; name, defaults, continuous_events)    
end

function get_iL_anomaly(t, autumn_start, autum_end)
    ty = t - floor(t)
    #dlit = shifloNormal(autumn_start, autum_end)
    dlit = Uniform(autumn_start, autum_end)
    # at the upper border return already 0
    ty == autum_end ? 0.0 :  pdf(dlit, ty) 
    #(1-share_autumn) + pdf(dlit, ty) * share_autumn
end
@register_symbolic get_iL_anomaly(t, autumn_start, autum_end)


function f_interpolate9(t, sol, numref)
	sol(t,idxs=numref.x)
end
@register f_interpolate9(t, sol::ODESolution, numref::Base.RefValue{Num})


function plant_face_fluct_fake(;name, sys, sol, t1, t2)
    @parameters t 
    D = Differential(t)
    ps = @parameters t1=t1 t2=t2 # for continuous events
    sts = @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t),
        i_L_annual(t)
    end) 
    eqs = [
        i_L ~  f_interpolate9(t, sol, Ref(sys.i_L)), 
        β_Ni ~ f_interpolate9(t, sol, Ref(sys.β_Ni)), 
        i_IN ~  f_interpolate9(t, sol, Ref(sys.i_IN)), 
        u_PlantNmax ~ f_interpolate9(t, sol, Ref(sys.u_PlantNmax)), 
        k_PlantN ~ f_interpolate9(t, sol, Ref(sys.k_PlantN)), 
        i_L_annual ~  f_interpolate9(t, sol, Ref(sys.i_L_annual)), 
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
    ODESystem(eqs,t, sts, ps; name, continuous_events)    
end


