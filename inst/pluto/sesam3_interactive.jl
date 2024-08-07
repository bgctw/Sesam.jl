### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes")].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 0232943f-74f8-40ec-96a3-495d5942e897
if occursin("Sesam", Base.current_project())
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.develop("MTKHelpers")
    Pkg.instantiate()
    using PlutoLinks: @revise
    using Sesam, MTKHelpers
    # currently using @revise enters an endless loop of reloading
    #@revise using Sesam
    #@revise using MTKHelpers
else
    using Sesam, MTKHelpers
end

# ╔═╡ d6f85fa6-97cc-40e9-841a-053563713993
begin
    using Plots, PlutoUI
    import PlutoUI: combine
    gr()
end;

# ╔═╡ 2c1e57e8-2069-4268-ab11-9b514054c9f6
using Suppressor

# ╔═╡ 01be1fe6-39f5-4ecb-9d9e-3c5bcf39c936
using ModelingToolkit, OrdinaryDiffEq, DataFrames, Tables, Distributions, Chain

# ╔═╡ f1a94d8f-927e-40f3-9af6-667890a6e733
using DistributionFits

# ╔═╡ 032c42f2-5103-11eb-0dce-e7ec59924648
html"<button onclick='present()'>present</button>"

# ╔═╡ 827fbc3a-512d-11eb-209e-cd74ddc17bae
begin
    struct TwoColumn{A, B}
        left::A
        right::B
    end
    function Base.show(io, mime::MIME"text/html", tc::TwoColumn)
        write(io,
            """
            <div style="display: flex;">
            	<div style="flex: 50%;">
            """)
        show(io, mime, tc.left)
        write(io,
            """
            	</div>
            	<div style="flex: 50%;">
            """)
        show(io, mime, tc.right)
        write(io,
            """
           		</div>
           	</div>
           """)
    end
end

# ╔═╡ 6e60edbd-82eb-4283-91ec-8fcb33a64ff7
md"## Testing local image"

# ╔═╡ a9ff79e4-f0dd-4c23-8755-379ab84475bf
pwd()

# ╔═╡ a392d095-19f6-495f-8c59-a93eace2b9a6
@suppress LocalResource("fig/logo.png", :width => 200)

# ╔═╡ 6382a73e-5102-11eb-1cfb-f192df63435a
md"""
# Sesam3 interactive exploration
**Thomas Wutzler 2022**
"""

# ╔═╡ 26a15193-a85f-443c-b007-2d14915e69f7
md"## Setup prior parameter distributions"

# ╔═╡ 078f11d4-9630-4e5c-8048-ef0c0232294f
md"Define the System components"

# ╔═╡ 36b322a5-ef4a-4697-90a8-720d8c224219
begin
    tspinup = 2000.0
    tface = 100.0
    #@named s = sesam3(use_seam_revenue=true)
    @named s = sesam3CN(use_seam_revenue = false)
    @named pl = plant_face(t1 = 0.0, t2 = tface)
    @named sp = plant_sesam_system(s, pl)
end;

# ╔═╡ 351ebf17-102c-4497-aab0-1e0a3f55f799
md"Define parameters and their distributions"

# ╔═╡ 7c778bf7-de69-4c14-99b6-b76cdf06596c
begin
    p = pC = Dict(s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
        #which respire and corresponding amount of N must be mineralized)
        s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
        # added to assimilable, i.e. uptake instead of R)
        s.a_E => 0.001 * 365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
        s.m => 0.005 * 365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
        s.τ => 1 / 60 * 365,  ##<< biomass turnover rate (12 days)    
        s.k_L => 1.0,       ##<< 1/(x years)   
        #s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
        s.k_R => 1 / (40.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
        s.k_mN_L => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
        s.k_mN_R => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
        s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
        #i_L => t -> 1 - exp(-t),  # litter input
        pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
        #pl.β_Ni0 => 25.0,
        pl.β_Ni0 => 30.0,
        #pl.i_IN0 => 0.0,   ##<< input of mineral N,
        pl.i_IN0 => 0.7,   ##<< 7kg/ha/yr
        #pl.k_Lagr => 12/2, # above ground litter turnover of 2 month
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
end;

# ╔═╡ 8247cdc0-8b23-4b1f-a4ee-6e6ac09a2232
begin
    cols = (:par, :dType, :mode, :upper)
    paramsModeUpperRows = [
        (s.β_NB, LogNormal, 8.0, 16.0),
        (s.β_NEnz, LogNormal, 3.0, 3.5),
        (s.k_mN_L, LogNormal, 60 * 0.05, 120 * 2.0),
        (s.k_mN_R, LogNormal, 60 * 0.05, 120 * 2.0),
        (s.κ_E, LogNormal, 0.7, 0.9),
        (s.k_R, LogNormal, 1 / 100, 1 / 10),
        (s.k_L, LogNormal, 12 / 18, 12 / 4), # tvr time between 18 months and 4 months
        (s.a_E, LogNormal, 0.001 * 365, 0.005 * 365),
        (s.m, LogNormal, 0.005 * 365, 0.02 * 365),
        (s.τ, LogNormal, 1.0 / 60 * 365, 1.0 / 5 * 365),
        (s.ϵ, LogitNormal, 0.5, 0.7),
        (s.ϵ_tvr, LogitNormal, 0.3, 0.9),
        (s.i_BN, LogNormal, 0.4, 4),
        (s.l_N, LogNormal, 0.96, 5 * 0.96),
        (s.ν_N, LogNormal, 0.9, 0.99),
        #(:kIPlant, LogNormal, 10.57 , 20)
    ]
    df_dist = rename!(DataFrame(columntable(paramsModeUpperRows)), collect(cols))
end;

# ╔═╡ ad1b457f-f98e-4556-83f2-477e2b857a07
#df_dist

# ╔═╡ 383cfe5b-922e-4151-9375-6093523d7f03
begin
    f1v = (dType, mode, upper) -> fit(dType, @qp_m(mode), @qp_uu(upper))
    transform!(df_dist, Cols(:dType, :mode, :upper) => ByRow(f1v) => :dist)
end;

# ╔═╡ 5c6d02b6-0aed-4ad3-843f-7e46367171cf
begin
    transform!(df_dist, :dist => ByRow(d -> quantile(d, 0.025)) => :lower)
    transform!(df_dist, :par => ByRow(pname -> p[pname]) => :default)
    select!(df_dist, Not(:dist), :dist)
end

# ╔═╡ d45a7ad5-4803-4846-ad2e-a9d4fd294413
md"""
## Modifying parameters

Select the parameters to explore with sliders below. (Below the plot)
"""

# ╔═╡ 1d6e347e-1a1b-4f2a-bc48-524c8d36f771
md"## Plotting the solution"

# ╔═╡ 7f1a89ce-1b84-45f0-af85-87b1050b18ba
begin
    # tstart_scrubb = @bind tstart NumberField(tstart_coarse.+(-20:20),default=tstart_coarse)
    # tend_scrubb = @bind tend NumberField(tend_coarse.+(-20:20),default=tend_coarse)	
    tstart_scrubb = @bind tstart NumberField((-tspinup):(tface + 20), default = -2)
    tend_scrubb = @bind tend NumberField((-tspinup):(tface + 20), default = +10)
end;

# ╔═╡ 9f251cbd-777e-4904-814d-6cef789fc560
@bind b_settimes Button("Set times programmatically")

# ╔═╡ 4b08433e-b8d8-4342-9153-2e4dad20356e
md"""
tstart: $(tstart_scrubb), tend: $(tend_scrubb)
"""

# ╔═╡ fdcefaa1-7ebf-49a1-ade5-00302ddeb6c9
@bind vars_disp MultiCheckBox(vcat(unknowns(sp), getproperty.(observed(sp), :lhs)))

# ╔═╡ ba4b67ab-b3fe-44f8-aeaf-925f76b005b5
md"""
Try modifying 
- β_NB between 11.6 and 11.8
- s.ϵ_tvr from 0.5 to 0.509
"""

# ╔═╡ 10214a12-e0e9-40aa-8a65-e8a9661125f6
@bind pars_mod MultiCheckBox(df_dist.par)

# ╔═╡ 9fe4cc1e-59d7-43c1-b10c-764756f5b4e8
function parameter_input(pars_mod)
    return PlutoUI.combine() do Child
        inputs = [begin
            row = subset(df_dist, :par => ByRow(x -> symbol(x) == symbol(name)))[1, :]
            el = Scrubbable(range(row.lower, row.upper, length = 100);
                default = row.default,
                format = ".3")
            # el = NumberField(range(first(row.lower), first(row.upper), length = 30);default=first(row.mode))			
            md""" $(name): $(
            	Child(name, el)
            )"""
        end
                  for name in pars_mod]
        md"""
        $(inputs) 
        """
    end
end

# ╔═╡ 854ab5c5-a572-41e8-90ee-bf8a63d7adf9
@bind popt parameter_input(string.(pars_mod))

# ╔═╡ 66dae3d2-9ba5-491b-a3cd-77fcfbb0ab78
ps = ODEProblemParSetter(sp, propertynames(popt))

# ╔═╡ a8b21bd2-3914-46b2-8b11-f9279f88e732
begin
    tspan_sim = (-tspinup, tface) # simulate 500 yrs spinup, increase at yr 20
    prob0 = prob = ODEProblem(sp, u0, tspan_sim, p) #2ms
    if length(popt) != 0
        prob = remake(prob0, ps, values(popt))
        p_upd = label_par(ps, prob.p)
    end
    sol0t = sol = solve(prob, Tsit5())
end;

# ╔═╡ 2ddb2869-4917-4973-8b3f-51dce7eba302
length(vars_disp) != 0 && plot(sol0t, tspan = (tstart, tend), vars = vars_disp)

# ╔═╡ b565d593-1369-4efe-8704-90f1611526c4
propertynames(popt)

# ╔═╡ 61e14646-568d-49a9-8d5e-afe9555e29ff
popt

# ╔═╡ 0b578017-ed4f-4dbf-80c8-7c5492799d17
plot(sol0t,
    tspan = (tstart, tend),
    vars = [s.dec_R p_upd[symbol(s.ϵ_tvr)] * s.tvr_B +
                    (1 - p_upd[symbol(s.κ_E)]) * s.tvr_Enz],
    legend = :bottomright);
vline!([0], linestyle = :dash, color = :lightgray, label = nothing);

# ╔═╡ 2448b012-864d-4797-b768-33174f9fd08f
plot(sol0t, tspan = (tstart, tend), vars = [s.α_L, s.α_R], legend = :bottomright);
vline!([0], linestyle = :dash, color = :lightgray, label = nothing);

# ╔═╡ 773273a8-70fc-4897-8174-8b06aff47013
plot(sol0t, tspan = (tstart, tend), vars = [s.lim_C, s.lim_N]);
vline!([0], linestyle = :dash, color = :lightgray, label = nothing);

# ╔═╡ 0cbe8944-5ebd-4804-863d-e567c703b89f
plot(sol0t, tspan = (tstart, tend), vars = [s.I_N]);
vline!([0], linestyle = :dash, color = :lightgray, label = nothing);

# ╔═╡ 1b5413de-3d32-4086-b80c-11c10ee87dac
plot(sol0t, tspan = (tstart, tend), vars = [s.Φ_NB, s.Φ_N - s.i_L / s.β_Ni]);
hline!([0], linestyle = :dash, color = :lightgray, label = nothing);
vline!([0], linestyle = :dash, color = :lightgray, label = nothing);

# ╔═╡ 255018b9-28c0-442c-8fa9-9698e7d7b513
md"## Thanks"

# ╔═╡ Cell order:
# ╠═0232943f-74f8-40ec-96a3-495d5942e897
# ╟─d6f85fa6-97cc-40e9-841a-053563713993
# ╟─032c42f2-5103-11eb-0dce-e7ec59924648
# ╟─827fbc3a-512d-11eb-209e-cd74ddc17bae
# ╠═6e60edbd-82eb-4283-91ec-8fcb33a64ff7
# ╠═a9ff79e4-f0dd-4c23-8755-379ab84475bf
# ╠═2c1e57e8-2069-4268-ab11-9b514054c9f6
# ╠═a392d095-19f6-495f-8c59-a93eace2b9a6
# ╟─6382a73e-5102-11eb-1cfb-f192df63435a
# ╟─26a15193-a85f-443c-b007-2d14915e69f7
# ╠═01be1fe6-39f5-4ecb-9d9e-3c5bcf39c936
# ╠═f1a94d8f-927e-40f3-9af6-667890a6e733
# ╟─078f11d4-9630-4e5c-8048-ef0c0232294f
# ╠═36b322a5-ef4a-4697-90a8-720d8c224219
# ╟─351ebf17-102c-4497-aab0-1e0a3f55f799
# ╟─7c778bf7-de69-4c14-99b6-b76cdf06596c
# ╠═8247cdc0-8b23-4b1f-a4ee-6e6ac09a2232
# ╠═ad1b457f-f98e-4556-83f2-477e2b857a07
# ╠═383cfe5b-922e-4151-9375-6093523d7f03
# ╟─5c6d02b6-0aed-4ad3-843f-7e46367171cf
# ╟─d45a7ad5-4803-4846-ad2e-a9d4fd294413
# ╟─a8b21bd2-3914-46b2-8b11-f9279f88e732
# ╟─66dae3d2-9ba5-491b-a3cd-77fcfbb0ab78
# ╠═b565d593-1369-4efe-8704-90f1611526c4
# ╟─1d6e347e-1a1b-4f2a-bc48-524c8d36f771
# ╟─7f1a89ce-1b84-45f0-af85-87b1050b18ba
# ╠═9f251cbd-777e-4904-814d-6cef789fc560
# ╠═4b08433e-b8d8-4342-9153-2e4dad20356e
# ╠═fdcefaa1-7ebf-49a1-ade5-00302ddeb6c9
# ╠═2ddb2869-4917-4973-8b3f-51dce7eba302
# ╟─ba4b67ab-b3fe-44f8-aeaf-925f76b005b5
# ╠═10214a12-e0e9-40aa-8a65-e8a9661125f6
# ╟─854ab5c5-a572-41e8-90ee-bf8a63d7adf9
# ╟─61e14646-568d-49a9-8d5e-afe9555e29ff
# ╟─9fe4cc1e-59d7-43c1-b10c-764756f5b4e8
# ╠═0b578017-ed4f-4dbf-80c8-7c5492799d17
# ╠═2448b012-864d-4797-b768-33174f9fd08f
# ╠═773273a8-70fc-4897-8174-8b06aff47013
# ╠═0cbe8944-5ebd-4804-863d-e567c703b89f
# ╠═1b5413de-3d32-4086-b80c-11c10ee87dac
# ╟─255018b9-28c0-442c-8fa9-9698e7d7b513
