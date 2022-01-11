@named s = sesam1()
@named pl = plant_const()

@parameters t
D = Differential(t)

# @named s = Sesam.sesam1C()
# sp = compose(ODESystem([
#     s.i_L ~ pl.i_L,
#     # colimitation
#     s.syn_B ~ s.C_synBC,
#   ], t; name=:sp), s, pl)

sp = compose(ODESystem([
    s.i_L ~ pl.i_L,
    s.i_IN ~ pl.i_IN,
    s.β_Ni ~ pl.β_Ni,
    s.u_PlantNmax ~ pl.u_PlantNmax,
    s.k_PlantN ~ pl.k_PlantN,
    # colimitation
    #s.syn_B ~ s.C_synBC,
    s.syn_B ~ min(s.C_synBC, s.β_NB*s.N_synBN), # TODO add P limitation
  ], t; name=:sp), s, pl)
sp = structural_simplify(sp)
states(sp)
equations(sp)

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
    s.k_mN => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    pl.β_Ni0 => 25,
    pl.i_IN0 => 0,   ##<< input of mineral N,
)
pN = Dict(
    s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NE => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
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
u0N = Dict(
    s.I_N => 0.04, ##<< inorganic pool gN/m2 
    s.L_N => u0[s.L]/p[pl.β_Ni0],
    s.R_N => u0[s.R]/p[s.β_NB],
    )
u0 = merge(u0C, u0N)    

tspan = (0.0,100.0)    

#prob = ODEProblem(sp,[t for t in u0], tspan, [t for t in p])
#prob = ODEProblem(sp, remove_units(u0), tspan, remove_units(p))
prob = ODEProblem(sp, u0, tspan, p)
#prob = ODEProblem(sp,u0, tspan, p, jac=true)
sol = solve(prob)

@testset "microbial dC balance" begin
    #plot(sol,  vars=[s.u_C, s.r_B+s.syn_B+s.syn_Enz])
    @test all(isapprox.(
        sol[s.u_C], 
        sol[s.r_B+s.syn_B+s.syn_Enz], 
        rtol = 1e-6))
end;

@testset "system dC balance" begin
    change = sol[s.dB + s.dL + s.dR]
    output = sol[s.r_tot] 
    input = sol[s.i_L]
    #plot(sol.t, [input, output, change, change .+ output])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;


@testset "microbial dN balance" begin
    uptake = p[s.ν_N]*sol[s.u_NOM]
    usage = sol[s.syn_Enz]/p[s.β_NE] + sol[s.syn_B]/p[s.β_NB] + sol[s.Φ_NB]
    #plot(sol.t, [uptake, usage])
    @test all(isapprox.(uptake, usage, rtol = 1e-6))
end;

@testset "system dN balance" begin
    change = sol[s.dB]/p[s.β_NB] + sol[s.dL_N + s.dR_N + s.dI_N]
    output = sol[s.u_PlantN + s.leach_N] 
    input = sol[s.i_L/s.β_Ni + s.i_IN]
    #sol[s.leach_N]
    #sol[s.i_IN]
    bo = 1:10
    bo = 1:length(input)
    #plot(sol.t[bo], [input[bo], output[bo], change[bo], (change .+ output)[bo]], label = ["input" "output" "change" "change+output"])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;

