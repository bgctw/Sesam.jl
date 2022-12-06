@named s = sesam3CN()
@named pl = plant_const()

# @named s = Sesam.sesam3CNC()
# sp = compose(ODESystem([
#     s.i_L ~ pl.i_L,
#     # colimitation
#     s.syn_B ~ s.C_synBC,
#   ], t; name=:sp), s, pl)
#sp = structural_simplify(sp)


@named sp = plant_sesam_system(s,pl)
states(sp)


sr = sesam3CN(;use_proportional_revenue=true, name=:s)
@named spr = plant_sesam_system(sr,pl)


p = pC = Dict(
    s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
        #which respire and corresponding amount of N must be mineralized)
    s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
        # added to assimilable, i.e. uptake instead of R)
    s.a_E => 0.001*365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    s.m => 0.005*365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    s.τ => 1/60*365,  ##<< biomass turnover rate (12 days)    
    s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    s.k_R => 1/(20.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    s.k_mN_L => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    s.k_mN_R => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    #s.ρ_CBtvr => 0.0,    # no carbon resorption
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    pl.β_Ni0 => 25,
    pl.i_IN0 => 0,   ##<< input of mineral N,
    #
    # P from plant model parameters not used in CN-Sesam soil model
    pl.β_Pi0 => Inf, #25*20, ## leaf litter N:P ~20(massratio Kang10)
    pl.i_IP0 => Inf, #0.65,   ##<< input of mineral P, weathering: Table3 mixed sedimentary rocks 0.65g/m2/yr Hartmann14 10.1016/j.chemgeo.2013.10.025
    pl.s_EP0 => Inf, # 0.5, # plant 1/20 of typical total microbial enzyme synthesis flux    
    pl.u_PlantPmax0 => Inf, 
    pl.k_PlantP0 => Inf,
)
pN = Dict(
    s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    #s.ρ_NBtvr => 0.0,    # no nitrogen resorption
    #s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    s.l_N => 0.0,       
    s.ν_N =>  0.9,     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    # for N uptake: take the defaults which take as much N as supplied by litter
    # pl.u_PlantNmax0 => Inf32, # only constrained by k_PlantN in min function
    # pl.k_PlantN0 => 10.57, #0.0289652*365     ##<< plant uptake rate first order of IN
)
p = merge(pC, pN)

u0 = u0C = Dict(
    s.B => 17,
    s.L => 100,
    s.R => 1100,
    #s.cumresp => 0.0,
    s.α_R => 0.1, # enzyme synthesis into L # TODO model by optimality
)
#u0C[s.α_L] = 1.0 - u0C[s.α_R]
u0N = Dict(
    s.I_N => 0.04, ##<< inorganic pool gN/m2 
    s.L_N => u0[s.L]/p[pl.β_Ni0],
    s.R_N => u0[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
    )
u0 = merge(u0C, u0N)    
#u0[s.R]/u0[s.R_N] # smaller p[s.β_NB]

tspan = (0.0,100.0)    
#tspan = (0.0,1e5)    
#tspan = (0.0,10.0)    

#prob = ODEProblem(sp,[t for t in u0], tspan, [t for t in p])
#prob = ODEProblem(sp, remove_units(u0), tspan, remove_units(p))
prob = ODEProblem(sp, u0, tspan, p)
#prob = ODEProblem(sp,u0, tspan, p, jac=true)
sol = sol_sesam3CN = solve(prob);

probr = ODEProblem(spr, u0, tspan, p)
solr = solve(probr);

i_plot = () -> begin
    #using Plots
    ts = (0,2)
    ts = extrema(sol.t)
    plot(sol)
    plot(sol, vars=[s.α_R], tspan=ts)
    plot!(solr, vars=[s.α_R], tspan=ts)
    plot(sol, vars=[s.B], tspan=ts)
    plot!(solr, vars=[s.B], tspan=ts)
    plot(sol, vars=[s.R], tspan=ts)
    plot(sol, vars=[s.du_L, s.du_R], tspan=ts)
    Plots.plot(sol, vars=[s.β_NBtvr, p[s.β_NB] * (1-p[s.ρ_CBtvr])/(1-p[s.ρ_NBtvr])], tspan=ts)
    calculate_β_NR_sesam3(p,s), sol[s.β_NR,end]
end

@testset "no N recycling" begin
    #include("test/test_sesam3CN_sol.jl")
    include("test_sesam3CN_sol.jl")
end;

@testset "N recycled" begin
    p2 = copy(p0)
    p2[s.ρ_NBtvr] = 0.1  # 10% N recycled during turnover
    prob2 = ODEProblem(sp, u0, tspan, p2)
    sol = sol2 = solve(prob2, Rodas4());
    #include("test/test_sesam3CN_sol.jl")
    include("test_sesam3CN_sol.jl")
end


