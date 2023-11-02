# important to initialize R_L correctly for steady plant, otherwise
# it consumes (initialized with too high C/N ratio) inorganic N

@named s = sesam3CN()
#@named s = sesam3(use_seam_revenue=true)
@named pl = plant_face(; t1 = 20, t2 = 70)

@named sp = plant_sesam_system(s, pl)
equations(sp)

p = pC = Dict(s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
    #which respire and corresponding amount of N must be mineralized)
    s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
    # added to assimilable, i.e. uptake instead of R)
    s.a_E => 0.001 * 365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    s.m => 0.005 * 365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    s.τ => 1 / 60 * 365,  ##<< biomass turnover rate (12 days)    
    s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    s.k_R => 1 / (20.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    s.k_mN_L => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.k_mN_R => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    #pl.β_Ni0 => 25.0,
    pl.β_Ni0 => 30.0,
    pl.i_IN0 => 0.0,   ##<< input of mineral N,
    #
    # P from plant model parameters not used in CN-Sesam soil model
    pl.β_Pi0 => Inf, #25*20, ## leaf litter N:P ~20(massratio Kang10)
    pl.i_IP0 => Inf, #0.65,   ##<< input of mineral P, weathering: Table3 mixed sedimentary rocks 0.65g/m2/yr Hartmann14 10.1016/j.chemgeo.2013.10.025
    pl.s_EP0 => Inf, # 0.5, # plant 1/20 of typical total microbial enzyme synthesis flux    
    pl.u_PlantPmax0 => Inf,
    pl.k_PlantP0 => Inf)
pN = Dict(s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    #s.l_N => 0.0,       
    s.ν_N => 0.9     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    # for N uptake: take the defaults which take as much N as supplied by litter
    # pl.u_PlantNmax0 => Inf32, # only constrained by k_PlantN in min function
    # pl.k_PlantN0 => 10.57, #0.0289652*365     ##<< plant uptake rate first order of IN
)
p = merge(pC, pN)

u0 = u0C = Dict(s.B => 17,
    #s.B => 34,
    s.L => 100,
    #s.L => 110,
    s.R => 1100,
    #s.R => 3250,
    #s.cumresp => 0.0,
    s.α_R => 0.5)
u0C[s.α_L] = 1.0 - u0C[s.α_R]
u0N = Dict(s.I_N => 0.04, ##<< inorganic pool gN/m2 
    s.L_N => u0[s.L] / p[pl.β_Ni0],
    #s.R_N => u0[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
    s.R_N => u0[s.R] / 7.0)
u0 = merge(u0C, u0N)
#u0[s.R]/u0[s.R_N] # smaller p[s.β_NB]

tspan_sim = (-500, 120.0)
prob = ODEProblem(sp, u0, tspan_sim, p)
#sol = solve(prob, saveat=vcat([450],500:0.5:600))
sol = solve(prob)

using Plots
tspan = first(sol.t) .+ (0, 5)
tspan = first(sol.t) .+ (0, 500)
tspan = (last(sol.t) - 5, last(sol.t))
tspan = (0, last(sol.t))
plot(sol; tspan)
plot(sol, vars = [s.L]; tspan)
plot(sol, vars = [s.R]; tspan)
plot(sol, vars = [s.R / s.β_NR]; tspan)
plot(sol, vars = [s.β_NR])
plot(sol, vars = [s.R + s.L]; tspan)
#plot!(sol, vars=[s.R + s.L]; tspan)

plot(sol, vars = [s.u_PlantN, s.u_PlantNmax]; tspan)
plot(sol, vars = [s.I_N]; tspan)

plot(sol, vars = [s.i_L / s.β_Ni, s.u_PlantN, s.Φ_N, s.Φ_Nu, s.Φ_NB, s.Φ_Ntvr]; tspan)

plot(sol, vars = [s.lim_C, s.lim_N]; tspan)
plot(sol, vars = [s.α_L, s.α_R]; tspan)
plot!(sol, vars = [s.α_L, s.α_R]; tspan)
plot(sol, vars = [s.invest_Ln, s.invest_Rn]; tspan)
plot!(sol, vars = [s.return_Ln, s.return_Rn]; tspan, linestyle = :dash)

plot(sol, vars = [s.α_RC, s.α_RN]; tspan)

sol.t
i = 1
sol[s.α_RC][i]
sol[s.α_RT][i]
