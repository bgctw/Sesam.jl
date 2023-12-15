using Sesam, Test
using Sesam: Sesam as CP
using ModelingToolkit, OrdinaryDiffEq
using ComponentArrays
using MTKHelpers
test_path = splitpath(pwd())[end] == "test" ? "." : "test"
@show test_path

include(joinpath(test_path,"testset_utils.jl"))

@named s = sesam3()
@named pl = plant_const_balanced()
@named plc = plant_const(; name = :pl)

@named sp = plant_sesam_system(s, pl)
states(sp)

sr = sesam3(; use_proportional_revenue = true, name = :s)
@named spr = plant_sesam_system(sr, pl)

i_inspect_equations = () -> begin
    tmp_s = Dict([e.lhs for e in equations(s)] .=> equations(s))
    tmp_pl = Dict([e.lhs for e in equations(pl)] .=> equations(pl))
    [k ∈ keys(tmp_s) for k in keys(tmp_pl)]
    tmp = merge(tmp_s, tmp_pl)
    tmp3 = sort(string.(states(s)))
    tmp2 = sort(replace.(string.(keys(tmp)), "Differential(t)(" => ""))
    tmp4 = [(s, eq) for (s, eq) in zip(tmp3[1:length(tmp3)], tmp2)]
    tmp4[(1 * 20) .+ (1:20)]
    tmp4[(4 * 20) .+ (1:20)]
    tmp4[100:end]
end

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
    #s.k_mN_R => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    s.ρ_CBtvr => 0.0,    # no carbon resorption
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    pl.β_Ni0 => 25,
    pl.i_IN0 => 0,   ##<< input of mineral N, 
    pl.β_Pi0 => 25 * 20, ## leaf litter N:P ~20(massratio Kang10)
    pl.i_IP0 => 0.65,   ##<< input of mineral P, weathering: Table3 mixed sedimentary rocks 0.65g/m2/yr Hartmann14 10.1016/j.chemgeo.2013.10.025
    pl.s_EP0 => 0.5)
pC[s.k_mN_R] = pC[s.k_mN_L]
pN = Dict(s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    s.ρ_NBtvr => 0.0,    # no nitrogen resorption
    #s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    s.l_N => 0.0,
    s.ν_N => 0.9)
pP = Dict(
    #s.k_LP => pC[s.k_L], # TODO 1/x years, assume same rate as depolymerizing enzyme
    #s.k_RP => pC[s.k_R], # TODO
    s.i_BP => pN[s.i_BN],       ##<< potential immobilization flux rate 
    s.β_PEnz => 50.0,     # TODO Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_PB => 40.0, # Sterner02: low P in microbial cell walls, more in genetic machinery and energy compounds
    s.ρ_PBtvr => 0.0,    # no phosphorus resorption
    s.l_P => 0.01,      # almost no leaching       
    s.ν_P => 0.9,     # microbial P use efficiency accounting for apparent 
    # mineralization, pertains only to SOM depolimerization, biomineralization is 
    # all mineralized
    s.β_Pm => 500 # at a c:P ratio of 500 depolimerization rate decreases to 1/2
    # s.k_mN_Pl => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
)
p = p0 = CP.get_updated_Penz_pars(merge(pC, pN, pP), s)

u0 = u0C = Dict(s.B => 17,
    s.L => 100,
    s.R => 1100,
    #s.cumresp => 0.0,
    s.α_R => 0.1,
    s.α_P => 0.1)
#u0C[s.α_L] = 1.0 - u0C[s.α_R] - u0C[s.α_LP] - u0C[s.α_RP]
u0N = Dict(s.I_N => 0.04, ##<< inorganic pool gN/m2 
    s.L_N => u0[s.L] / p[pl.β_Ni0],
    s.R_N => u0[s.R] / calculate_β_NR_sesam3(p, s))
u0P = Dict(s.I_P => 0.04, ##<< TODO inorganic pool gN/m2 
    s.L_P => u0[s.L] / p[pl.β_Pi0],
    s.R_P => u0[s.R] / calculate_β_PR_sesam3(p, s))
u00 = u0 = merge(u0C, u0N, u0P)
#u0[s.R]/u0[s.R_N] # smaller p[s.β_NB]

tspan = (0.0, 200.0)

#prob = ODEProblem(sp,[t for t in u0], tspan, [t for t in p])
#prob = ODEProblem(sp, remove_units(u0), tspan, remove_units(p))
prob = ODEProblem(sp, u0, tspan, p)
#prob = ODEProblem(sp,u0, tspan, p, jac=true)
#sol = sol_sesam3 = solve(prob, Tsit5());
sol = sol_sesam3 = solve(prob, Rodas4());
#sol = sol_sesam3 = solve(prob, Tsit5(), callback=PositiveDomain(prob.u0));

probr = ODEProblem(spr, u0, tspan, p)
solr = solve(probr, Rodas4());

i_plot = () -> begin
    #import StatsPlots
    #using Plots
    ts = tspan
    ts = extrema(sol.t)
    ts = (0, 0.06)
    ts = (2, 30)
    Plots.plot(sol)
    Plots.plot(sol, idxs = [s.R])
    Plots.plot(sol, idxs = [s.calculate_β_PR_sesam3])
    Plots.plot(sol, idxs = [s.lim_C, s.lim_N, s.lim_P], tspan = ts)
    Plots.plot(sol, idxs = [s.α_L, s.α_R, s.α_P], tspan = ts)
    Plots.plot(sol, idxs = [s.d_L, s.d_R, s.d_P], tspan = ts)
    Plots.plot(sol, idxs = [s.du_L, s.du_R, s.du_P, s.mdu], tspan = ts)
    Plots.plot(sol, idxs = [s.dα_R, s.dα_P], tspan = ts)

    Plots.plot!(solr, idxs = [s.α_L, s.α_R, s.α_P], tspan = ts)
    Plots.plot(sol,
        idxs = [s.revenue_L, s.revenue_R, s.revenue_LP, s.revenue_RP],
        tspan = ts)
    Plots.plot(sol, idxs = [sr.return_L, sr.return_R, sr.return_P], tspan = ts)
    Plots.plot(sol, idxs = [sr.dec_RP_P, sr.dec_PPlant], tspan = ts)
    Plots.plot(sol, idxs = [s.I_P])
    Plots.plot(sol, idxs = [s.u_PlantP, s.u_immPPot])
    Plots.plot(sol, idxs = [s.p_uPmic])

    Plots.plot(sol, idxs = [p[s.a_E] * s.B])
    Plots.plot(sol, idxs = [s.B])
    Plots.plot(sol, idxs = [s.α_LP, s.α_RP])
    Plots.plot(sol, idxs = [s.α_LP, s.α_RP], tspan = (68, 92))
    Plots.plot(sol, idxs = [s.α_LPT, s.α_RPT], tspan = (80, 92))

    Plots.plot(sol, idxs = [s.β_PR], tspan = ts)
    Plots.plot(sol,
        idxs = [
            s.tvr_B / p[s.β_PB],
            s.tvr_B / p[s.β_PBtvr],
            s.resorp_P,
            s.tvr_B / p[s.β_PBtvr] + s.resorp_P,
        ],
        tspan = ts)

    Plots.plot(sol,
        idxs = [s.β_PBtvr, p[s.β_PB] * (1 - p[s.ρ_CBtvr]) / (1 - p[s.ρ_PBtvr])],
        tspan = ts)
    calculate_β_PR_sesam3(p, s), sol[s.β_PR, end]
end

@testset "calculate_β_NR_sesam3 symbols" begin
    β_NR = calculate_β_NR_sesam3(p, s)
    p_sym = Dict(Symbol(k) => v for (k, v) in p)
    β_NR2 = calculate_β_NR_sesam3(p_sym)
    @test β_NR2 == β_NR
end;

@testset "calculate_β_PR_sesam3 symbols" begin
    β_PR = calculate_β_PR_sesam3(p, s)
    p_sym = Dict(Symbol(k) => v for (k, v) in p)
    β_PR2 = calculate_β_PR_sesam3(p_sym)
    @test β_PR2 == β_PR
end;

@testset "calculate_propR" begin
    p_sym = ComponentVector(pl₊β_Ni0 = 60.0, s₊β_NB = 16, s₊β_NEnz = 3.1,
        s₊a_E = 0.365, s₊κ_E = 0.8, s₊ϵ_tvr = 0.45, s₊τ = 6.1)
    pR = calculate_propR(p_sym, 25.0)
    #@test pR ≈ 0.51 atol=0.005 # without accounting for enzyme contribution
    @test pR≈0.44 atol=0.005
end;

@testset "compute_mean_du3" begin
    # all three
    @test CP.compute_mean_du3(0.7, 1 / 3, 0.4, 1 / 3, 0.4, 1 / 3) == 0.5
    # only du1 and du2
    # mean is slightly less than 0.5 to account for decreasing third
    @test 0.4 < CP.compute_mean_du3(0.6, 0.45, 0.4, 0.45, 0.01, 0.1) < 0.5
    # only du1 and du3
    @test 0.4 < CP.compute_mean_du3(0.6, 0.45, 0.01, 0.1, 0.4, 0.45) < 0.5
    # only du3 and du3
    @test 0.4 < CP.compute_mean_du3(0.01, 0.1, 0.6, 0.45, 0.4, 0.45) < 0.5
    # only du1
    # slightly less than 0.7 to account for decreasing others
    @test 0.5 < CP.compute_mean_du3(0.7, 0.9, 0.1, 0.1, 0.1, 0.1) < 0.7
    # only du2
    @test 0.5 < CP.compute_mean_du3(0.1, 0.1, 0.7, 0.9, 0.1, 0.1) < 0.7
    # only du3
    @test 0.5 < CP.compute_mean_du3(0.1, 0.1, 0.1, 0.1, 0.7, 0.9) < 0.7
end;

@testset "no P recycled" begin
    #include("test/test_sesam3_sol.jl")
    include(joinpath(test_path,"test_sesam3_sol.jl"))
end;

@testset "P recycled" begin
    p2 = copy(p)
    p2[s.ρ_PBtvr] = 0.1  # 10% N recycled during turnover
    prob2 = ODEProblem(sp, u0, tspan, p2)
    sol = sol2 = solve(prob2, Rodas4())
    include(joinpath(test_path,"test_sesam3_sol.jl"))
end;

@testset "dalpha on P gradient" begin
    # example from sesam_LRP_deriv.Rmd
    dL = 0.7
    dR = 0.5
    dP = 1.0
    B0 = 1.0
    s_EP = 0.0 #0.2
    kwargs = (;)
    sLRP = CP.sesam_const_dLRP(dL, dR, dP; name = :s, kwargs...)
    @named spLRP = plant_sesam_system(sLRP, pl)
    sLRP_r = CP.sesam_const_dLRP_relative(dL, dR, dP; name = :s, kwargs...)
    @named spLRP_r = plant_sesam_system(sLRP_r, pl)
    states(spLRP_r)
    u0LRP = copy(u0)
    u0LRP[s.B] = B0
    pLRP = copy(p0)
    pLRP[pl.s_EP0] = s_EP
    pLRP[s.a_E] = 0.1
    pLRP[s.k_mN_L] = pLRP[s.k_mN_R] = pLRP[s.k_mN_P] = pLRP[s.a_E] * B0 / 2
    probLRP = ODEProblem(spLRP, u0LRP, (0, 5), pLRP)
    probLRP_r = ODEProblem(spLRP_r, u0LRP, (0, 5), pLRP)
    sol = solLRP_r = solve(probLRP, Rodas4());
    sol = solLRP = solve(probLRP_r, Rodas4());
    #push!(LOAD_PATH, "/User/homes/twutz/julia/scimltools/")
    #using ComponentArrays
    #using MTKHelpers
    popt = ComponentVector(state = (s₊B = 1,), par = (s₊d_L0 = dL, s₊d_R0 = dR, s₊d_P0 = dP))
    pset = ODEProblemParSetter(spLRP, popt)
    pset_r = ODEProblemParSetter(spLRP_r, popt)
    #pset_αr = ODEProblemParSetter(spLRP_r, ComponentVector(state=ComponentVector(), par=(s₊α_R = 1 / 3, s₊α_P = 1 / 3)))

    calc_alpha3_proptoderiv = (dL, dR, dP, B0 = 1, s_EP = 0; use_proportional_revenue = false, kwargs...) -> begin
        popt = ComponentVector(state=(s₊B = B0,), par=(s₊d_L0 = dL, s₊d_R0 = dR, s₊d_P0 = dP))
        prob_u = use_proportional_revenue ?
                 remake(probLRP_r, popt, pset_r) :
                 remake(probLRP, popt, pset)
        get_paropt_labeled(pset, prob_u)
        get_paropt_labeled(pset, probLRP)
        sol = solve(prob_u, Rodas4(), saveat = [prob_u.tspan[2]])
        xE = map(v -> sol[v][end], [s.α_L, s.α_R, s.α_P])
        xE
    end
    # calc_alpha3_proptoderiv(dL, dR, 0)
    # regression to values of sesam_LRP_deriv.Rmd
    @test isapprox(calc_alpha3_proptoderiv(dL, dR, 0)[1], 0.5839202, atol = 1e-5)
    @test isapprox(calc_alpha3_proptoderiv(dL, dR, 0)[3], 0.0, atol = 1e-5)
    @test isapprox(calc_alpha3_proptoderiv(dL, dR, 0.3)[3], 0.15, atol = 0.01)
    @test isapprox(calc_alpha3_proptoderiv(dL, dR, 1)[3], 0.48, atol = 0.01)
    #
    @test isapprox(calc_alpha3_proptoderiv(dL, dR, 0; use_proportional_revenue = true)[1],
        0.56,
        atol = 0.01)
    @test isapprox(calc_alpha3_proptoderiv(dL, dR, 0.3; use_proportional_revenue = true)[3],
        0.23,
        atol = 0.01)

    tmpf = () -> begin
        #using Plots
        #using LaTeXStrings
        dPs = range(0, stop = 1, length = 51)[2:end]
        tmp = map(dPi -> calc_alpha3_proptoderiv(dL, dR, dPi), dPs)
        plot(dPs, getindex.(tmp, 1), label = "α_L")
        plot!(dPs, getindex.(tmp, 2), label = "α_R")
        plot!(dPs, getindex.(tmp, 3), label = "α_P")
        #L is a macro and must be resolved also during precompilation
        # but LaTeXStrings should not be part of Project.toml
        #xlabel!("potential biominaralization "*L"d_P \ (g \, m^{-2} \, yr^{-1})")
        xlabel!("potential biominaralization " * "d_P (g/m2/yr)")
        tmp2 = map(dPi -> calc_alpha3_proptoderiv(dL,
                dR,
                dPi;
                use_proportional_revenue = true),
            dPs)
        plot!(dPs, getindex.(tmp2, 1), label = "α_L rel", linestyle = :dot)
        plot!(dPs, getindex.(tmp2, 2), label = "α_R rel", linestyle = :dot)
        plot!(dPs, getindex.(tmp2, 3), label = "α_P rel", linestyle = :dot)
        tmp3 = map(dPi -> calc_alpha3_proptoderiv(dL, dR, dPi, 5), dPs)
        plot!(dPs, getindex.(tmp3, 1), label = "α_L highB", linestyle = :dash)
        plot!(dPs, getindex.(tmp3, 2), label = "α_R highB", linestyle = :dash)
        plot!(dPs, getindex.(tmp3, 3), label = "α_P highB", linestyle = :dash)
        #
        tmp4 = map(dPi -> calc_alpha3_proptoderiv(dL,
                dR,
                dPi,
                5;
                use_proportional_revenue = true),
            dPs)
        plot!(dPs, getindex.(tmp4, 1), label = "α_L rel,highB", linestyle = :dashdot)
        plot!(dPs, getindex.(tmp4, 2), label = "α_R rel,highB", linestyle = :dashdot)
        plot!(dPs, getindex.(tmp4, 3), label = "α_P rel,highB", linestyle = :dashdot)
        tmp5 = map(dPi -> calc_alpha3_proptoderiv(dL,
                dR,
                dPi,
                0.1;
                use_proportional_revenue = true),
            dPs)
        plot!(dPs, getindex.(tmp5, 5), label = "α_L rel,lowB", linestyle = :dashdot)
        plot!(dPs, getindex.(tmp5, 2), label = "α_R rel,lowB", linestyle = :dashdot)
        plot!(dPs, getindex.(tmp5, 3), label = "α_P rel,lowB", linestyle = :dashdot)
    end
end;
