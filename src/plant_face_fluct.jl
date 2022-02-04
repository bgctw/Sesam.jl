"""
model i_L and β_Ni as as step function increased by factor fac_int at t in [t1,t2).
Superimpose a seasonal pattern on the litter input.
"""
function plant_face_fluct(;name, t1=0.0, t2=100.0, fac_inc=1.2, 
    autumn_start=10/12, autumn_end = 11/12, share_autumn=0.5, k_Lagr=12/2, k_PlantN0v = 100)
    @parameters t 
    D = Differential(t)
    @parameters i_L0  i_IN0  β_Ni0 t1=t1 t2=t2 fac_inc=fac_inc
    @parameters u_PlantNmax0  k_PlantN0=k_PlantN0v 
    @parameters k_Lagr = k_Lagr share_autumn = share_autumn

    @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        i_L_anomaly(t), i_L_annual(t), β_Ni_annual(t),
        Lagr(t), i_Lagr(t), dec_Lagr(t), d_Lagr(t),
        u_PlantNmax(t), k_PlantN(t)
    end) 
    eqs = [
        D(Lagr) ~ d_Lagr, d_Lagr ~ i_Lagr - dec_Lagr,
        dec_Lagr ~ k_Lagr*Lagr,
        i_Lagr ~ share_autumn * i_L_annual * i_L_anomaly,
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
        [t-yr ~ autumn_start*1.05 for yr in -500:200],
        [t-yr ~ autumn_end*0.995 for yr in -500:200],
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

# function get_from_sol(sol::ODESolution, t, sym::Num) 
#     #sin(t)
#     @show t, sym
#     sol(t, idxs=sym)
# end
# @register  get_from_sol(sol::ODESolution, t, sym::Num) 

function get_from_sol(t, sol, aterm) 
    sol(t)
end
@register_symbolic get_from_sol(t, sol, term) 
#@register_symbolic (::ODESolution)(t, term) 

# function plant_face_fluct_fake(;name, sys, sol)
#     @parameters t 
#     D = Differential(t)
#     @variables (begin
#         i_L(t), β_Ni(t), i_IN(t),
#         i_L_anomaly(t), i_L_annual(t), β_Ni_annual(t),
#         Lagr(t), i_Lagr(t), dec_Lagr(t), d_Lagr(t),
#         u_PlantNmax(t), k_PlantN(t)
#     end) 
#     #@unpack i_L, β_Ni, i_IN, i_L_anomaly, i_L_annual, β_Ni_annual, Lagr, i_Lagr, dec_Lagr, d_Lagr, u_PlantNmax, k_PlantN = sys
#     eqs = [
#         D(tmp) ~ 0,
#         i_L ~ get_from_sol(sol, t, sys.i_L),
#         dec_Lagr ~ get_from_sol(sol, t, sys.dec_Lagr),
#         i_Lagr ~ get_from_sol(sol, t, sys.i_Lagr),
#         i_L ~ get_from_sol(sol, t, sys.i_L),
#         i_L_annual ~ get_from_sol(sol, t, sys.i_L_annual), 
#         β_Ni_annual ~ get_from_sol(sol, t, sys.β_Ni_annual),
#         i_L_anomaly ~ get_from_sol(sol, t, sys.i_L_anomaly),
#         β_Ni ~ get_from_sol(sol, t, sys.β_Ni), 
#         i_IN ~ get_from_sol(sol, t, sys.i_IN), 
#         u_PlantNmax ~ get_from_sol(sol, t, sys.u_PlantNmax),
#         k_PlantN ~ get_from_sol(sol, t, sys.k_PlantN),
#     ]

#     # continuous_events = vcat(
#     #     [t-yr ~ autumn_start*1.05 for yr in -500:200],
#     #     [t-yr ~ autumn_end*0.995 for yr in -500:200],
#     #     [
#     #     t ~ t1,
#     #     t ~ t2,
#     #     ])
#     ODESystem(eqs,t; name)    
# end


function f_interpolate(t, table_t, table_u)
	interpolator = LinearInterpolation(table_u, table_t)
	output = interpolator(t)
end
@register f_interpolate(t, table_t::AbstractVector, table_u::AbstractVector)


function f_interpolate2(t, interp)
	interp(t)
end
@register f_interpolate2(t, interp::DataInterpolations.AbstractInterpolation)


function f_interpolate4(t, sol, sys)
	#sol(t,idxs=sys.i_L)
	sol(t,idxs=getproperty(sys, :i_L))
end
@register f_interpolate4(t, sol::ODESolution, sys::ODESystem)

# function f_interpolate5(t, sol, term)
# 	sol(t,idxs=term)
# end
# @register f_interpolate5(t, sol::ODESolution, term::Term)

# function f_interpolate6(t, sol, termdict, termsym)
#     @variables (begin
#         i_L(t), β_Ni(t), i_IN(t),
#         u_PlantNmax(t), k_PlantN(t)
#     end) 
#     term = termdict[termsym]
#     @show term, typeof(term)
# 	sol(t,idxs=term)
# end
# @register f_interpolate6(t, sol::ODESolution, termdict::Dict, termsym::Symbol)

# function f_interpolate7(t, sol, termdict, termsym, pl)
#     @show termsym
#     @show pl.i_L
#     term = termdict[termsym]
#     @show term, typeof(term)
# 	sol(t,idxs=term)
# end
# @register f_interpolate7(t, sol::ODESolution, termdict::Dict, termsym::Symbol, pl::ODESystem)

function f_interpolate8(t, sol, sys, symref)
	#sol(t,idxs=sys.i_L)
	sol(t,idxs=getproperty(sys, symref.x))
end
@register f_interpolate8(t, sol::ODESolution, sys::ODESystem, symref::Base.RefValue{Symbol})



function plant_face_fluct_fake(;name, sys, sol)
    @parameters t 
    D = Differential(t)
    sts = @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t)
    end) 
    # @named syse = embed_system(sys)
    # terms = Symbolics.value.(getproperty.(observed(syse),:lhs))
    # term_syms = [t.f.name for t in terms]
    # term_dict = Dict(x => terms[findfirst(term_syms .== x)] for x in 
    #     (:pl₊i_L, :pl₊β_Ni, :pl₊i_IN, :pl₊u_PlantNmax, :pl₊k_PlantN))
    term_dict = Dict(x => getproperty(sys, x) for x in 
    (:i_L, :β_Ni, :i_IN, :u_PlantNmax, :k_PlantN))
    #@show term_dict[:i_L], typeof(term_dict[:i_L])
    #sol(0, idxs=first(values(term_dict)))
    #sol(0, idxs=collect(values(term_dict))[1])


    #sol(0, idxs=sys.i_L)
    
    # ts = sol.t
    # i_L_ts = sol(ts, idxs=sys.i_L)
    # interp = LinearInterpolation(i_L_ts, ts)
    # local sysc = deepcopy(sys)
    # Sesam.f_interpolate4(0, sol, sysc)
    eqs = [
        #i_L ~  f_interpolate(t, ts, i_L_ts),
        #i_L ~  f_interpolate2(t, interp),
        #i_L ~  f_interpolate4(t, sol, sys), 
        i_L ~  f_interpolate8(t, sol, sys, Ref(:i_L)), 
        #i_L ~  f_interpolate6(t, sol, term_dict, :pl₊i_L), 
        # β_Ni ~ f_interpolate6(t, sol, term_dict, :pl₊β_Ni), 
        # i_IN ~  f_interpolate6(t, sol, term_dict, :pl₊i_IN), 
        # u_PlantNmax ~ f_interpolate6(t, sol, term_dict, :pl₊u_PlantNmax), 
        # k_PlantN ~ f_interpolate6(t, sol, term_dict, :pl₊k_PlantN), 
        #i_L ~  f_interpolate7(t, sol, term_dict, :i_L, sys), 
        # β_Ni ~ 1, #f_interpolate6(t, sol, term_dict, :β_Ni), 
        # i_IN ~  1, #f_interpolate6(t, sol, term_dict, :i_IN), 
        # u_PlantNmax ~ 1, #f_interpolate6(t, sol, term_dict, :u_PlantNmax), 
        # k_PlantN ~ 1, #f_interpolate6(t, sol, term_dict, :k_PlantN), 
        β_Ni ~ f_interpolate8(t, sol, sys, Ref(:β_Ni)), 
        i_IN ~  f_interpolate8(t, sol, sys, Ref(:i_IN)), 
        u_PlantNmax ~ f_interpolate8(t, sol, sys, Ref(:u_PlantNmax)), 
        k_PlantN ~ f_interpolate8(t, sol, sys, Ref(:k_PlantN)), 

        # β_Ni ~ get_from_sol(sol, t, β_Ni), 
        # i_IN ~ get_from_sol(sol, t, i_IN), 
        # # i_L_anomaly ~ get_from_sol(sol, t, i_L_anomaly),
        # # i_L_annual ~ get_from_sol(sol, t, i_L_annual), 
        # # β_Ni_annual ~ get_from_sol(sol, t, β_Ni_annual),

        # # i_Lagr ~ get_from_sol(sol, t, i_Lagr),
        # # dec_Lagr ~ get_from_sol(sol, t, dec_Lagr),
        # # dec_Lagr ~ get_from_sol(sol, t, d_Lagr),
        # u_PlantNmax ~ get_from_sol(sol, t, u_PlantNmax),
        # k_PlantN ~ get_from_sol(sol, t, k_PlantN),
    ]
    # continuous_events = vcat(
    #     [t-yr ~ autumn_start*1.05 for yr in -500:200],
    #     [t-yr ~ autumn_end*0.995 for yr in -500:200],
    #     [
    #     t ~ t1,
    #     t ~ t2,
    #     ])
    ODESystem(eqs,t, sts, []; name)    
end


