@named s = sesam3()
@named pl = plant_const()

@named sp = plant_sesam_system(s,pl)
states(sp)

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
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    pl.β_Ni0 => 25,
    pl.i_IN0 => 0,   ##<< input of mineral N, 
    pl.β_Pi0 => 25*20, ## leaf litter N:P ~20(massratio Kang10)
    pl.i_IP0 => 0.65,   ##<< input of mineral P, weathering: Table3 mixed sedimentary rocks 0.65g/m2/yr Hartmann14 10.1016/j.chemgeo.2013.10.025
    pl.s_EP0 => 0.5, # plant 1/20 of typical total microbial enzyme synthesis flux
)
pN = Dict(
    s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    #s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    s.l_N => 0.0,       
    s.ν_N =>  0.9,     # microbial N use efficiency accounting for apparent 
)
pP = Dict(
    s.k_LP => pC[s.k_L], # TODO 1/x years, assume same rate as depolymerizing enzyme
    s.k_RP => pC[s.k_R], # TODO
    s.i_BP => pN[s.i_BN],       ##<< potential immobilization flux rate 
    s.β_PEnz => 50.0,     # TODO Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_PB => 40.0, # Sterner02: low P in microbial cell walls, more in genetic machinary and energy compounds
    s.l_P => 0.01,      # almost no leaching       
    s.ν_P =>  0.9,     # microbial P use efficiency accounting for apparent 
    # mineralization, pertains only to SOM depolimerization, biomineralization is 
    # all mineralized
    # s.k_mN_Pl => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
)
p = merge(pC, pN, pP)

u0 = u0C = Dict(
    s.B => 17,
    s.L => 100,
    s.R => 1100,
    #s.cumresp => 0.0,
    s.α_R => 0.1, 
    s.α_RP => 0.1, 
    s.α_LP => 0.1, 
)
u0C[s.α_L] = 1.0 - u0C[s.α_R] - u0C[s.α_LP] - u0C[s.α_RP]
u0N = Dict(
    s.I_N => 0.04, ##<< inorganic pool gN/m2 
    s.L_N => u0[s.L]/p[pl.β_Ni0],
    s.R_N => u0[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
    )
u0P = Dict(
    s.I_P => 0.04, ##<< TODO inorganic pool gN/m2 
    s.L_P => u0[s.L]/p[pl.β_Pi0],
    s.R_P => u0[s.R]/calculate_β_PR_sesam3(p,s) #p[s.β_NB],
    )
u0 = merge(u0C, u0N, u0P)    
#u0[s.R]/u0[s.R_N] # smaller p[s.β_NB]

tspan = (0.0,200.0)    

#prob = ODEProblem(sp,[t for t in u0], tspan, [t for t in p])
#prob = ODEProblem(sp, remove_units(u0), tspan, remove_units(p))
prob = ODEProblem(sp, u0, tspan, p)
#prob = ODEProblem(sp,u0, tspan, p, jac=true)
#sol = sol_sesam3 = solve(prob);
sol = sol_sesam3 = solve(prob, Tsit5());
#sol = sol_sesam3 = solve(prob, Tsit5(), callback=PositiveDomain(prob.u0));

i_plot = () -> begin
    #using Plots
    ts = tspan
    ts = (2,2.2)
    plot(sol)
    plot(sol, vars=[s.R])
    plot(sol, vars=[s.calculate_β_PR_sesam3])
    plot(sol, vars=[s.lim_C, s.lim_N, s.lim_P], tspan=ts)
    plot(sol, vars=[s.α_L, s.α_R, s.α_LP, s.α_RP])
    plot(sol, vars=[s.α_L, s.α_R, s.α_LP, s.α_RP], tspan=ts)
    plot(sol, vars=[s.revenue_L, s.revenue_R, s.revenue_LP, s.revenue_RP], tspan=ts)
    plot(sol, vars=[s.dec_RP, s.dec_RPPlant], tspan=ts)
    plot(sol, vars=[s.I_P])
    plot(sol, vars=[s.u_PlantP, s.u_immPPot])
    plot(sol, vars=[s.p_uPmic])

    plot(sol, vars=[p[s.a_E] * s.B])
    plot(sol, vars=[s.α_LP, s.α_RP])
    plot(sol, vars=[s.α_LP, s.α_RP], tspan=(68,92))
    plot(sol, vars=[s.α_LPT, s.α_RPT], tspan=(80,92))
end

@testset "non-negative pools" begin
    #st = first(states(sp))
    for st in states(sp)
        #@show st
        #@test all(sol[st] .>= 0.0) 
        @test all(sol[st] .>= -eps(Float64)*100) 
    end
    # minimum(sol[s.α_LP])
    # minimum(sol[s.α_RP])
end;

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
    usage = sol[s.syn_Enz]/p[s.β_NEnz] + sol[s.syn_B]/p[s.β_NB] + sol[s.Φ_NB]
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

@testset "balanced allocation sums to one" begin
    #plot(sol, vars=[s.α_LT,s.α_RT, s.α_LT+s.α_RT])
    @test all(isapprox.(sol[s.α_LT + s.α_RT + s.α_LPT + s.α_RPT], 1.0, rtol = 1e-8))
    @test all(isapprox.(sol[s.α_L + s.α_R + s.α_LP + s.α_RP], 1.0, rtol = 1e-8))
end;

@testset "microbial dP balance" begin
    uptake = p[s.ν_P]*sol[s.u_POM];
    usage = sol[s.syn_Enz]/p[s.β_PEnz] + sol[s.syn_B]/p[s.β_PB] + sol[s.Φ_PB];
    #plot(sol.t, [uptake, usage])
    @test all(isapprox.(uptake, usage, rtol = 1e-6))
end;

@testset "system dP balance" begin
    change = sol[s.dB]/p[s.β_PB] + sol[s.dL_P + s.dR_P + s.dI_P]
    output = sol[s.u_PlantP + s.leach_P] 
    input = sol[s.i_L/s.β_Pi + s.i_IP]
    #sol[s.leach_N]
    #sol[s.i_IN]
    bo = 1:10
    bo = 1:length(input)
    #plot(sol.t[bo], [input[bo], output[bo], change[bo], (change .+ output)[bo]], label = ["input" "output" "change" "change+output"])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;

