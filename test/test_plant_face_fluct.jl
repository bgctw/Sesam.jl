@parameters t
D = Differential(t)
@variables x(t) # integrator

@named plf = plant_face_fluct()

function rep_plf(pl;name, simplify=true)
    @parameters t
    D = Differential(t)
        sts = @variables x(t)
    # integrate of i_L to check annual litter inputs
    tmp = compose(ODESystem([
        D(x) ~ pl.i_L
        ],t,sts,[];name), pl)
    simplify ? structural_simplify(tmp) : tmp
 end
 @named plf_rep = rep_plf(plf)
# @named plf_rep = embed_system(plf)
 equations(plf_rep)
 
p_plf = Dict(
    plf.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    plf.β_Ni0 => 25,
    plf.i_IN0 => 0,   ##<< input of mineral N,
    plf.t1 => 0.0,
    plf.t2 => 3.0,
)
u_plf = Dict(
    plf.Lagr => 0.2,
    x => 0.0,
)
#tspan = (-500.0, 120.0)
tspan = (-5.0, 5.0)
prob = ODEProblem(plf_rep, u_plf, tspan, p_plf)
#sol = solve(prob, saveat = 0:1:tspan[2])
#solver = CompositeAlgorithm((Tsit5(),Vern7()), (integrator) -> (Int(integrator.dt<0.001) + 1))
solver = Vern7() # implicit method
# solver = AutoVern9(Rodas4())
# solver = Vern9()
#using  deSolveDiffEq; solver = deSolveDiffEq.lsoda()

# solver = Tsit5() # takes forever and produces wrong solution
sol = solve(prob, solver, tspan=tspan,
    reltol = 1e-5,  # simulate litterfall with higher accuracy
    #progress = true, progress_steps = 1,
    #saveat = 0.005 # uncomment for proper plotting
);
# sol = sol_euler = solve(prob, Euler(), tspan=tspan,
#     dt = 0.005,
#     progress = true, progress_steps = 1,
#     #saveat = 0.005
# );

i_plot = () -> begin
    # using Plots
    ts = tspan
    Plots.plot(sol, tspan=ts, vars=[plf.i_L, plf.i_L_annual], xlab="time (yr)", ylab="litter input (g/m2/yr)")
    Plots.plot(sol, tspan=ts, vars=[plf.β_Ni], xlab="time (yr)", ylab="CN litter input (g/g)")
    Plots.plot(sol, tspan=ts, vars=[plf.inc_period], xlab="time (yr)", ylab="inc_period")
    series_sol(sol, )

    ts = (-2,5)
    ts = (-0.1,+0.1)
    ts = (0.7,1)
    ts = (-0.3,0)
    ts = (0,10) .+ tspan[1]
    ts = tspan
    Plots.plot(sol, tspan=ts, vars=[plf.i_Lagr, plf.dec_Lagr, plf.d_Lagr])

    x = (9.5/12:1/(12*10):11.5/12) .+ 1
    Plots.scatter(x, first.(sol.(x)))
    plot(sol, tspan=extrema(x), vars=[plf.i_Lagr, plf.dec_Lagr])

    plot(sol, vars=[plf.i_L / plf.β_Ni])
end

@testset "integration of litterfall over one year" begin
    isapprox(sol(0, idxs=x) - sol(-1, idxs=x), 400, atol=0.05)
    isapprox(sol(5, idxs=x) - sol(4, idxs=x), 480, atol=0.05)
end;

i_plot_integrated_annual = () -> begin
    Plots.plot(sol, vars=[x])
    ts = range(first(tspan)+1, last(tspan), length=200)
    xann = sol.(ts, idxs=x) .- sol.(ts .- 1.0, idxs=x)
    Plots.plot(ts, xann)
end


# dlit = shifloNormal(10/12, 11/12)
# #plot(dlit, xlim=(0,1))
# plot(dlit)
# ts = 0:1/24:3
# t = 10.5*1/12Sophienstraße 2, 07743 Jena 
# ts = (9.5:0.05:11.5)/12
# plot(ts*12, iL.(ts,1))

# # check that integral matches annual litterfall
# using QuadGK
# lann = 3.2
# f1(x) = iL(x,lann)
# @test isapprox(first(quadgk(f1, 0, 1, rtol=1e-3)), lann, rtol=1e-2)

i_inspect_dlit = () -> begin
    autumn_start=8.5/12; autumn_end=11.5/12
    d_lit_agr = LocationScale(autumn_start, autumn_end - autumn_start, 
        #     LogitNormal(0,sqrt(2)), 
        fit_mode_flat(LogitNormal, 0.3; peakedness = 3)
        #    Uniform()
    )
    xs = 0.6:0.005:1
    Plots.plot(xs .* 12 .+ 1, pdf.(d_lit_agr, xs))
end

