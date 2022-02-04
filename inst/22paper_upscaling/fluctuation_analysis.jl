# inspect effect of averaging litter input compare to fluctuating litter input

using Sesam
using ModelingToolkit, DifferentialEquations
using DataFrames, Tables
using Distributions
using Chain
using MTKHelpers

# name both systems s, so that can use same keys in parameters
ss = sesam3(;name=:s)
s = se = seam3(;name=:s)
@named pl = plant_face_fluct()

p = pC = Dict(
    s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
        #which respire and corresponding amount of N must be mineralized)
    s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
        # added to assimilable, i.e. uptake instead of R)
    s.a_E => 0.001*365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    s.m => 0.005*365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    s.τ => 1/60*365,  ##<< biomass turnover rate (12 days)    
    s.k_L => 1.0,       ##<< 1/(x years)   
    #s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    s.k_R => 1/(40.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    ss.k_mN => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    s.k_N => 60, 
    s.k_m => 0.05, 
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    #pl.β_Ni0 => 25.0,
    pl.β_Ni0 => 30.0,
    #pl.i_IN0 => 0.0,   ##<< input of mineral N,
    pl.i_IN0 => 0.7,   ##<< 7kg/ha/yr
    pl.k_Lagr => 12/2, # above ground litter turnover of 2 month
    pl.k_PlantN0 => 100.0, 
    #pl.k_PlantN0 => 2, # try with a lower rate - cannot resupply
)
pN = Dict(
    s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    #s.l_N => 0.0,       
    s.ν_N =>  0.9,     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    # for N uptake: take the defaults which take as much N as supplied by litter
    # pl.u_PlantNmax0 => Inf32, # only constrained by k_PlantN in min function
    # pl.k_PlantN0 => 10.57, #0.0289652*365     ##<< plant uptake rate first order of IN
)
p = merge(pC, pN)

get_litter_input_fake_system = () -> begin
    #pls = structural_simplify(pl)
    #@unpack Lagr, i_L0, β_Ni0, i_IN0 = pls
    #ps = Dict(i_L0 => p[pl.i_L0], β_Ni0  => p[pl.β_Ni0], i_IN0 => p[pl.i_IN0])
    #@named plse = embed_system(pl)
    prob_pl = ODEProblem(plse, [pl.Lagr => p[pl.i_L0]/2 / p[pl.k_Lagr]], (-500,200), p)
    sol_pl = solve(prob_pl, Vern7())
    plf = plant_face_fluct_fake(;name=pl.name,sys=pl, sol=sol_pl);
end
plf = get_litter_input_fake_system()

#sim_u0steady = () -> begin
    u0C = Dict(
        s.B => 17,
        s.L => 100,
        s.R => 1100,
        #s.cumresp => 0.0,
        s.α_R => 0.1, # enzyme synthesis into L # TODO model by optimality
        pl.Lagr => p[pl.i_L0]/2 / p[pl.k_Lagr]
    )
    u0C[s.α_L] = 1.0 - u0C[s.α_R]
    u0C2 = Dict(
        s.E_L => u0C[s.α_L] * p[s.a_E] * u0C[s.B] / (p[s.k_N]),
        s.E_R => u0C[s.α_R] * p[s.a_E] * u0C[s.B] / (p[s.k_N]),
    )
    u0N = Dict(
        s.I_N => 0.04, ##<< inorganic pool gN/m2 
        s.L_N => u0C[s.L]/p[pl.β_Ni0],
        s.R_N => u0C[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
        )
    u0p = Dict(
        pl.i_L => 0.0, pl.β_Ni => 0.0, pl.u_PlantNmax => 0.0, pl.k_PlantN => 0.0, pl.i_IN => 0.0, #pl.s_one => 1.0,
    )
    u00 = merge(u0C, u0C2, u0N, u0p)    
    #u00[s.R]/u00[s.R_N] # smaller p[s.β_NB]
    tspan = (0.0,500.0)
    @named ssp0 = plant_sesam_system(se,plf)
    prob0 = ODEProblem(ssp0, u00, tspan, p)
    # pss_tmp = ProblemParSetter(ssp0, (pl.share_autumn, pl.Lagr))
    # prob0s = update_statepar(pss_tmp, (0.0, 0.0), prob0)
    #prob0s = deepcopy(prob0)
    sol = sol_sesam3s0 = solve(prob0);
    k = first(states(ssp0))
    u0s = Dict(k => sol[k,end] for k in states(ssp0))
    #print(u0s) # outputs code that can b pasted for u0
end

u0 = sim_u0steady()
#u0 = merge(u0, Dict(s.I_N => u0[s.I_N]*20))

#---------- simulate sesam3 with varying litter inputs
using CairoMakie, AlgebraOfGraphics
set_aog_theme!()  # to get the consistent colors

import IterTools
using DataFrames, Chain
variants = @chain IterTools.product((:seam, :sesam),(:seasonal, :annual)) begin
    DataFrame()
    rename!([:enzyme,:litter])
    leftjoin(DataFrame(enzyme=[:seam, :sesam], linestyle=[:dash, :solid]), on=:enzyme)
    transform(:enzyme => (x -> Makie.current_default_theme().palette.color[][1:length(x)]) => :color)
    transform(AsTable(Cols(:enzyme,:litter)) => ByRow(t -> "$(t.enzyme)_$(t.litter)") => :label)
end
variants[!,:sol] = Vector{Any}(fill(missing,nrow(variants)))
typeof(variants.sol)

tspan = (-80.0,50.0)
solver = Vern7() # implicit method
#solver = Vern9() # implicit method
#solver = Rodas5()
@named sep = plant_sesam_system(se,pl)
probe = ODEProblem(sep, u0, tspan, p)
@named ssp = plant_sesam_system(ss,pl)
probs = ODEProblem(ssp, u0, tspan, p)

# first both seam3 and sesam3 with fluctuating litter input
sol = variants[findfirst(variants.label .== "seam_seasonal"),:sol] = sol_seam3f = 
    solve(probe, solver; tspan, abstol=1e-8);
sol = variants[findfirst(variants.label .== "sesam_seasonal"),:sol] = sol_sesam3f = 
    solve(probs, solver; tspan, abstol=1e-8);

plot(sol_seam3f, vars=[s.I_N])
plot(sol_sesam3f, vars=[s.I_N])

tend = maximum(sol_sesam3f.t)
ts = tend .+ (-0.5, +0.5)
tse = tend .+ (-0.5, 0.0)
plot(sol_sesam3f, tspan=(-75,maximum(sol_sesam3f.t)), vars=[s.I_N])
plot!(sol_seam3f, tspan=(-75,maximum(sol_sesam3f.t)+1), vars=[s.I_N])

plot(sol_sesam3f, tspan=ts, vars=[s.I_N])
plot!(sol_sesam3f, tspan=tse, vars=[s.I_N], xlim=ts)
plot!(sol_sesam3f, tspan=ts, vars=[s.i_L], xlim=ts)
plot!(sol_seam3f, tspan=ts, vars=[s.i_L], xlim=ts)

plot(sol_seam3f, tspan=ts, vars=[s.α_R], xlim=ts)
plot!(sol_sesam3f, tspan=tse, vars=[s.α_R], xlim=ts)

plot(sol_seam3f, tspan=ts, vars=[s.E_R], xlim=ts)
plot!(sol_sesam3f, tspan=tse, vars=[s.E_R], xlim=ts)

plot(sol_seam3f, tspan=ts, vars=[s.C_synBC, s.C_synBN], xlim=ts)
plot!(sol_sesam3f, tspan=tse, vars=[s.C_synBC, s.C_synBN], xlim=ts)

plot(sol_seam3f, tspan=ts, vars=[s.w_C, s.w_N], xlim=ts)
plot!(sol_sesam3f, tspan=tse, vars=[s.w_C, s.w_N], xlim=ts)



plot(sol_sesam3f, tspan=(-75,maximum(sol_sesam3f.t)), vars=[s.R])
plot(sol_sesam3f, tspan=(-75,maximum(sol_sesam3f.t)), vars=[s.u_PlantN, s.u_PlantNmax])

plot(sol_sesam3s, tspan=(-75,maximum(sol_sesam3f.t)), vars=[s.u_PlantN, s.u_PlantNmax])
plot(sol_seam3s, tspan=(-75,maximum(sol_sesam3f.t)+1), vars=[s.u_PlantN, s.u_PlantNmax])
plot(sol_sesam3s, tspan=(-10,0), vars=[s.I_N])
plot(sol_sesam3s, tspan=(-10,0), vars=[s.leach_N])
plot(sol_sesam3s, tspan=(-10,0), vars=[s.R/s.R_N])
plot(sol_sesam3s, vars=[s.R])



# next parameterize plant model for non-fluctuating litter input
# this is achieved by setting the share of autumn input to zero
pse_s = ProblemParSetter(sep, (pl.share_autumn, pl.Lagr))
probe_s = update_statepar(pse_s, (0.0, 0.0), probe);
pss_s = ProblemParSetter(ssp, (pl.share_autumn, pl.Lagr))
probs_s = update_statepar(pss_s, (0.0, 0.0), probs);

sol = variants[findfirst(variants.label .== "seam_annual"),:sol] = sol_seam3s = 
    solve(probe_s, solver; tspan, abstol=1e-8);
sol = variants[findfirst(variants.label .== "sesam_annual"),:sol] = sol_sesam3s = 
    solve(probs_s, solver; tspan, abstol=1e-8);




function plot_vars(vars, tspan=tspan; kwargs...)
    pl = plot(sol_seam3f; tspan=ts, vars=vars, label="seam", linestyle=:dash, kwargs...)
    plot!(pl, sol_seam3s; tspan=ts, vars=vars, label="seam annual", linestyle=:dash, kwargs...)
    plot!(pl, sol_sesam3s; tspan=ts, vars=vars, label="sesam annual", kwargs...)
    plot!(pl, sol_sesam3f; tspan=ts, vars=vars, label="sesam", xlab="time (yr)", kwargs...)    
end


cm2inch(x) = x/2.54
cm2inch.((8.3,7))
inch2px(x, dpi=Plots.DPI) = x * dpi
dpi_pub = 300
#plot!(dpi=dpi_pub, size = inch2px.(cm2inch.((8.3,7)), dpi_pub))
figpath = "inst/22paper_upscaling"

i_plot_Plots = () -> begin
    #using Plots
    pla = (
        palette = :Dark2_6,
        dpi = dpi_pub,
        size = inch2px.(cm2inch.((8.3,7)), dpi_pub),
    )
    ts = (-5.0,min(5.0,maximum(sol_sesam3f.t)))
    plot(sol_seam3f; tspan=ts, vars=[pl.i_L, pl.i_L_annual], xlab="time (yr)", 
        ylab="litter input (g/m2/yr)", pla...)
    savefig(joinpath(figpath,"fluct_litterinput.pdf"))

    plot_vars([s.E_R]; ylab="enzyme pool E_R (g/m2)", tspan=(-2,0), legend=:topleft, pla...)
    savefig(joinpath(figpath,"fluct_E_R.pdf"))

    plot_vars([s.leach_N]; ylab="N leaching (g/m2/yr)", tspan=ts, legend=:topleft, pla...)
    savefig(joinpath(figpath,"fluct_Nleach.pdf"))
end

sol = sol_sesam3f
ts = (-5.0,min(5.0,maximum(sol.t)))
tsi = (1:length(sol.t))[first(ts) .<= sol.t .<= last(ts)]

# vars=[pl.i_L, pl.i_L_annual]
using CairoMakie, AlgebraOfGraphics
set_aog_theme!()
legend_theme = Theme(
    Legend=(
        bgcolor=(:white,0.8),
        framevisible=true,
        rowgap=0.2,
        framewidth=0.2,
        framecolor=(:lightgrey,0.8),
        padding=(5.0,5.0,2.0,2.0),
        patchsize=(20.0,10.0),
        margin=(2,2,2,2),
))
update_theme!(legend_theme)
update_theme!(
    figure_padding = 1,
    Lines = (linewidth = 0.8,),    
    Series = (linewidth=0.8,),
    )
axis_theme = Theme(
    Axis=(
        #xlabelsize=20, 
        xgridvisible=true,
        ygridvisible=true,
        #xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1.0, ytickalign=1.0, yticksize=3, xticksize=3,
        xticklabelpad=0,
        xlabelpadding=0,
        #ylabelpadding=0,
        #leftspinevisible = false,
        rightspinevisible = true,
        # bottomspinevisible = false,
        topspinevisible = true,        
        topspinecolor = :darkgray,
        rightspinecolor = :darkgray,
    )    
)    
update_theme!(axis_theme)


#Makie.current_default_theme().palette.color[]
#Makie.current_default_theme().Lines
#Makie.current_default_theme().Legend
#Makie.current_default_theme().Axis

function plotm_vars!(ax, vars, tspan=tspan; variants = variants, legend_position=:rr, kwargs...)
    #row = first(Tables.namedtupleiterator(select(variants, :sol, :label, :linestyle, :color)));
    for row in Tables.namedtupleiterator(select(variants, :sol, :label, :linestyle, :color))
        @show typeof(row)
        series_sol!(ax, row.sol, vars; tspan, labels=[row.label], 
        linestyle=row.linestyle,
        solid_color=row.color,
        )
    end
    axislegend(ax, unique=true, position = legend_position)
    display(fig)
    fig
end
#save("tmp.pdf", fig, pt_per_unit = 1)

fig, ax = pdf_figure2(xlabel = "Time (yr)", ylabel="enzyme pool E_R (g/m2)");
plotm_vars!(ax, [s.E_R], (-2,0); variants = variants[[1,2],:], legend_position=:lt)
savefig(joinpath(figpath,"fluct_E_R.pdf"))

fig, ax = pdf_figure2(xlabel = "Time (yr)", ylabel="litter input (g/m2/yr)");
ts = (-5.0,min(5.0,maximum(sol_seam3f.t)))
series_sol!(ax, sol_seam3f, [pl.i_L, pl.i_L_annual], tspan=ts, labels=["seasonal","annual"], linewidth=0.8)
axislegend(ax, unique=true, valign = :top, halign=:left, margin=(2,2,2,2))
display(fig)
save(joinpath(figpath,"fluct_litterinput.pdf"), fig, pt_per_unit = 1)

fig, ax = pdf_figure2(xlabel = "Time (yr)", ylabel="N leaching (g/m2/yr)");
plotm_vars!(ax, [s.leach_N], (-20,5); variants = variants[[1,3,4],:], legend_position=:lt)
display(fig)
save(joinpath(figpath,"fluct_Nleach.pdf"), fig, pt_per_unit = 1)

fig, ax = pdf_figure2(xlabel = "Time (yr)", ylabel="R (g/m2)");
plotm_vars!(ax, [s.R], (first(tspan),0); variants = variants[[1,2,3,4],:], legend_position=:lt)
display(fig)

fig, ax = pdf_figure2(xlabel = "Time (yr)", ylabel="N leaching (g/m2/yr)");
plotm_vars!(ax, [s.leach_N], (-20,5); variants = variants[[1],:], legend_position=:lt)
display(fig)
series_sol!(ax, )

sol = sol_seam3f
plot(sol, tspan=(-80,-0), vars=[s.R])
i_plot = () -> begin
    # using Plots
    plot(sol, vars=[s.R])
    ts = tspan
    ts = (-5.0,min(5.0,maximum(sol_sesam3f.t)))
    plot(sol_seam3f; tspan=ts, vars=[pl.i_L, pl.i_L_annual], xlab="time (yr)", 
        ylab="litter input (g/m2/yr)", pla...)
    
    #plot(sol, tspan=ts, vars=[pl.Lagr])
        

    plot_vars([s.R]; ylab="residue pool R (g/m2)", tspan=ts, pla...)
    #plot_vars([s.dR]; ylab="change residue pool dR (g/m2)", tspan=ts)
    plot_vars([s.L]; ylab="litter pool L (g/m2)", tspan=ts, pla...)
    plot_vars([s.L+s.R]; ylab="SOM stocks (L+R) (g/m2)", tspan=ts, pla...)
    plot_vars([s.r_tot]; ylab="respiration (g/m2/yr)", tspan=ts, pla...)
    plot_vars([s.leach_N]; ylab="N leaching (g/m2/yr)", tspan=ts, legend=:bottomright, pla...)
    plot_vars([s.α_R]; ylab="enzyme allocation to E_R (g/g)", tspan=ts, pla...)
    plot_vars([s.E_R]; ylab="residue enzyme pool E_R (g/m2)", tspan=ts, pla...)

    plot(sol_sesam3f, vars = [pl.i_L, pl.i_L_annual]; ylab="litter input (g/m2/yr)", 
        tspan=ts, xlab="time (yr)")
    
    ts = (-2,0)
    # almost no difference between steady enzyme and explicit enzyme
    plot_vars([s.E_R]; ylab="residue enzyme pool E_R (g/m2)", tspan = ts, ylim=(0.065,0.075))

    ts = (-2,5)
    ts = (-1,1)
    ts = (0,10) .+ tspan[1]
    plot(sol, tspan=ts, vars=[pl.i_L, pl.i_L_annual])
    plot(sol, tspan=ts, vars=[pl.i_Lagr, pl.dec_Lagr, pl.d_Lagr])

    x = (9.5/12:1/(12*10):11.5/12) .+ 1
    scatter(x, first.(sol.(x)))
    plot(sol, tspan=extrema(x), vars=[pl.i_Lagr, pl.dec_Lagr])

    plot(sol, vars=[pl.i_L / pl.β_Ni])
end

using CairoMakie
lines(sol.t, sol[s.L])
scatter(sol_seam3f.t[1:20], sol_seam3f[s.L,1:20])

