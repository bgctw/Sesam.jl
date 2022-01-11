using Sesam
using ModelingToolkit, DifferentialEquations

@named s = sesam1()
@named pl = plant_const()

@parameters t
D = Differential(t)

sp = compose(ODESystem([
    s.i_L ~ pl.i_L,
    s.i_IN ~ pl.i_IN,
    s.β_Ni ~ pl.β_Ni,
    s.u_PlantNmax ~ pl.u_PlantNmax,
    s.k_PlantN ~ pl.k_PlantN,
  ], t; name=:sp), s, pl)
sp = structural_simplify(sp)
states(sp)
equations(sp)

p = Dict(
    s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
        #which respire and corresponding amount of N must be mineralized)
    s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
        # added to assimilable, i.e. uptake instead of R)
    s.a_E => 0.001*365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    s.m => 0.005*365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    s.τ => 1/60*365,  ##<< biomass turnover rate (12 days)    
    s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    s.k_R => 1/(20.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    s.k_mN => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    #
    # ----------- N -------------
    s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NE => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    #s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    s.l_N => 0.0,       
    s.ν_N =>  0.9,     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    pl.β_Ni0 => 25,
    pl.i_IN0 => 0,   ##<< input of mineral N,
    # for N uptake: take the defaults which take as much N as supplied by litter
    # pl.u_PlantNmax0 => Inf32, # only constrained by k_PlantN in min function
    # pl.k_PlantN0 => 10.57, #0.0289652*365     ##<< plant uptake rate first order of IN
)

u0 = Dict(
    s.B => 17,
    s.L => 100,
    s.R => 1100,
    #s.cumresp => 0.0,
    s.α_R => 0.1, # enzyme synthesis into L # TODO model by optimality
    s.I_N => 0.04, ##<< inorganic pool gN/m2 
)
u0[s.L_N] = u0[s.L]/p[pl.β_Ni0]
u0[s.R_N] = u0[s.R]/p[s.β_NB]

tspan = (0.0,100.0)    

#prob = ODEProblem(sp,[t for t in u0], tspan, [t for t in p])
#prob = ODEProblem(sp, remove_units(u0), tspan, remove_units(p))
prob = ODEProblem(sp, u0, tspan, p)
#prob = ODEProblem(sp,u0, tspan, p, jac=true)
sol = solve(prob)

using Plots
plot(sol, vars=[L,B])
plot(sol, vars=[I_N])
plot(sol)
plot(sol, vars=[β_NL])
# increase of C in system is the same as input
plot(sol, vars=[s.C_synBC, p[s.β_NB]*s.N_synBN, s.syn_B])
plot(sol, vars=[β_NL])


Dict(p)[i_L]
# total C stock equals input
#plot(sol, vars=[R + L + B + cumresp])
#plot!(t -> Dict(p)[i_L]*t + sum(getindex.(Ref(Dict(u0)),[B,L,R])))
# biomass C balance
# plot(sol, tspan=(0.0,2.0), vars=[u_C, r_B+syn_B+syn_Enz])
# plot(sol,  vars=[u_C, r_B+syn_B+syn_Enz])
#
# biomass N balance
 plot(sol, tspan=(0.0,2.0), vars=[u_NOM, r_B+syn_B+syn_Enz])
# plot(sol,  vars=[u_C, r_B+syn_B+syn_Enz])

DifferentialEquations v7.0.0
IfElse v0.1.1
ModelingToolkit v8.0.0
Unitful v1.10.1