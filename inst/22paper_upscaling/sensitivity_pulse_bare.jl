# inspect sensitivity of increase in SOM stocks after 50yrs of increased
# C input

#using Pkg; Pkg.activate("inst/22paper_upscaling")
using Sesam
#push!(LOAD_PATH, expanduser("~/julia/scimltools/")) # access local package repo
using ModelingToolkit, DifferentialEquations
using DataFrames, Tables
using Distributions
using Chain
using MTKHelpers
import ComponentArrays as CA
using OrderedCollections
import SubglobalSensitivityAnalysis as SSA

@named s = sesam3_revMM(use_seam_revenue = false)
@named pl = plant_const()
@named sp = plant_sesam_system(s, pl)

p = pC = Dict(s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
    #which respire and corresponding amount of N must be mineralized)
    s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
    # added to assimilable, i.e. uptake instead of R)
    s.a_E => 0.001 * 365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    s.m => 1.825, #0.005*365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    s.τ => 1.0 / 60 * 365,  ##<< biomass turnover rate (12 days)    
    s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    s.k_R => 1.0 / (20.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    s.k_mN_L => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.k_mN_R => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #
    pl.i_L0 => 0.0)
pN = Dict(s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    s.l_N => 0.0, #0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    #s.l_N => 0.0,       
    s.ν_N => 0.9,     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    # for N uptake: take the defaults which take as much N as supplied by litter
    pl.u_PlantNmax0 => 0.0,
    pl.k_PlantN0 => 0.0,
    pl.i_IN0 => 0.0,
    pl.i_IP0 => 0.0,
    #pl.β_Ni0, # need to be set dependent on scenario
    pl.β_Pi0 => NaN)
p = merge(pC, pN)

u0C = Dict(s.B => 0.3,
    #s.B => 34.0,
    s.L => 100.0,
    #s.L => 110.0,
    #s.R => 1200.0,
    s.R => 0.1,  # to avoid initial C limitation, start from near zeor R pool
    #s.R => 3250.0,
    #s.cumresp => 0.0,
    s.α_R => 0.5)
u0C[s.α_L] = 1.0 - u0C[s.α_R]

#----------------- compute range for given CN
convert_symbol(numdict::Dict{Num}) = Dict(symbol(k) => val for (k, val) in numdict)
function convert_symbol(numdict::OrderedDict{Num})
    OrderedDict(symbol(k) => val for (k, val) in numdict)
end
convert_symbol(p)

function get_u0(cn0, u0C, p_labeled)
    u0N = OrderedDict(
        #s.I_N => 1.0, ##<< inorganic pool gN/m2 
        s.I_N => 0.001, ##<< inorganic pool gN/m2 
        s.L_N => u0C[s.L] / cn0,
        #s.R_N => u0[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
        s.R_N => u0C[s.R] / p_labeled[:s₊β_NB],
        pl.β_Ni0 => 30.0)
    u0 = OrderedDict(merge(u0C, u0N))
end

#------- run a single simulation and extract range of CUE ----
cn0 = 50.0
u0_cn0 = get_u0(cn0, u0C, convert_symbol(p))
#tend = 4.0 # must be extended for lower k_L 
tend = 800.0 # will end earlier by callback if biomass gets very low
tspan_sim = (0, tend) # simulate 500 yrs spinup, increase at yr 20
prob0 = ODEProblem(sp, Dict(u0_cn0), tspan_sim, p)

# in index to update initial states and params specify distributions already
# in big parts redundant to sensitivity_face
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
using Random;
Random.seed!(0);
# :par should hold symbols instead of Nums
df_dist = SSA.fit_distributions(paramsModeUpperRows)
transform!(df_dist, :par => identity => :par_num)
transform!(df_dist, :par => ByRow(symbol) => :par)

names_omit = [:s₊l_N]
names_opt = setdiff(df_dist.par, names_omit)
#upd! = SystemParUpdater(collect(names_opt), sp)
pset = ProblemParSetter(sp, CA.Axis(names_opt))
popt0 = get_paropt_labeled(pset, prob0)
prob = update_statepar(pset, popt0, prob0)

u00 = get_u0(50, u0C, convert_symbol(p))
pset_state = ProblemParSetter(sp, CA.Axis(symbol.(keys(u00))))

cb_nonnegative_biomass = PositiveDomain(copy(prob0.u0))

# α_R can also become negative
# cb_nonnegative_biomass = ManifoldProjection() do resid,u,p,t
#     ul = label_state(pset, u)
#     residl = label_state(pset, resid)
#     residl .= zero(residl)
#     residl[:s₊B] = max(zero(eltype(residl)), -ul[:s₊B]) 
#     parent(residl) # return plan unerlying vector
# end

cb_terminate_lowbiomass = DiscreteCallback((u, t, integrator) -> label_state(pset, u)[:s₊B] <
                                                                 1e-8,
    integrator -> terminate!(integrator))

function simExp(popt; prob0 = prob0, cn0 = cn0)
    #saveat=[0.0,tend] # save at start of increase and tface years later
    #2ms
    #prob0 = ODEProblem(sp, u0, tspan_sim, p; jac=true) # actually slower: 4ms
    prob = update_statepar(pset, popt, prob0)
    popt_labeled = label_par(pset, prob.p)
    # also update u0 for given p 
    u0o_dict = get_u0(cn0, u0C, popt_labeled)
    prob = update_statepar(pset_state, collect(values(u0o_dict)), prob)
    #solve(prob, Rodas5(), callback = cb_nonnegative_biomass);
    solve(prob,
        Rodas5(),
        abstol = 1e-10,
        callback = CallbackSet(cb_nonnegative_biomass, cb_terminate_lowbiomass))
end

function compute_range(x)
    ex = extrema(x)
    ex[2] - ex[1]
end

function sim_range(popt)
    sol = simExp(popt)
    range_CUE = compute_range(sol[s.syn_B / s.u_C])
    (; range_CUE, sol)
end

(r0, sol) = sim_range(get_paropt(pset, prob0));

tmp_fig4_pulse50_90 = () -> begin
    # tmp = map([50,70,90]) do cn0
    #     u0_cn0 = get_u0(cn0, u0C, convert_symbol(p))    
    #     cn0, u0_cn0
    # end
    # df = rename!(DataFrame(tmp), [:cn0, :u0])
    df = DataFrame(cn0 = [50, 70, 90])

    #tmp = simExp(df.u0[1])
    prob4 = remake(prob0, tspan = (0, 4))
    transform!(df, :cn0 => ByRow(cn0 -> simExp(popt0; prob0 = prob4, cn0 = cn0)) => :sol)
    tmpf = () -> begin
        sol0t = df.sol[1]
        plot(sol0t, idxs = [s.L + s.R])
        plot(sol0t, idxs = [s.syn_B])
        plot!(df.sol[2], idxs = [s.syn_B])
        plt = plot(df.sol[1], idxs = [s.syn_B / s.u_C])
        for i in 2:nrow(df)
            print(i)
            plot!(plt, df.sol[i], idxs = [s.syn_B / s.u_C])
        end
        display(plt)
        plot(sol0t, idxs = [s.lim_C, s.lim_N])
        plot(sol, idxs = [s.lim_C, s.lim_N])
        plot(sol, idxs = [s.syn_B / s.u_C])
    end
    sol0 = solve(prob0; saveat)
end

tmp_experiment_simperiod = () -> begin
    tend = 4.0
    probe = remake(prob0, tspan = (0, tend))
    sol = simExp(popt0; prob0 = probe)
    plot(sol, idxs = [s.B])
    plot(sol, idxs = [s.L])
    plot(sol, idxs = [s.R])
    plot(sol, idxs = [s.L_N])
    plot(sol, idxs = [s.lim_C, s.lim_N])
    plot(sol, idxs = [s.lim_C, s.lim_N], tspan = (sol.t[end] - 2, sol.t[end]))
    plot(sol, idxs = [s.I_N])
    x = label_state(pset, sol[63])
end

#--------- sensitivity analysis ------------------
# https://bgctw.github.io/SubglobalSensitivityAnalysis.jl/dev/reload_design/
if !@isdefined N
    N = 100
    N = 5000
    #N = 5001 # replicate
    #N = 10_000
    #N = 25_000
end

δ_cp = 0.2
δ_cp = 0.1

fsens = (popt) -> begin
    local r, sol_p = sim_range(popt)
    sol_p.retcode ∈ (:Success, :Terminated) ? (; r) : missing
end

df_dist_opt = subset(df_dist, :par => ByRow(x -> x ∈ names_opt))
SSA.set_reference_parameters!(df_dist_opt, Dict(propertynames(popt0) .=> values(popt0)))

scen_str = "N$(N)_deltacp$(Integer(δ_cp*100))"
fname_sens = joinpath("tmp", "sens_pulse_bare_$(scen_str).rds")
fname_output = joinpath("tmp", "sens_pulse_bare_$(scen_str).feather")
fname_plot = joinpath("tmp", "sens_pulse_bare_$(scen_str).pdf")

estim_file = SSA.SobolTouati(; rest = SSA.RSobolEstimator("sens_touati", fname_sens))

SSA.calculate_parbounds!(df_dist_opt)
X1 = SSA.get_uniform_cp_sample(df_dist_opt, N);
X2 = SSA.get_uniform_cp_sample(df_dist_opt, N);
cp_design = SSA.generate_design_matrix(estim_file, X1, X2);
SSA.reload_design_matrix(estim_file);

q_design = SSA.transform_cp_design_to_quantiles(df_dist_opt, cp_design);
fsens(first(eachrow(q_design)))
res = map(r -> fsens(r), eachrow(q_design));
import Feather
Feather.write(fname_output, res)

y = [tup[:r] for tup in res];
df_sobol = SSA.estimate_sobol_indices(estim_file, y, df_dist_opt.par)
df_sobol[!, :target] .= :r;

#--------- plotting
using CategoricalArrays
using AlgebraOfGraphics

# plot the sensitivity
R"library(ggplot2)"
R"library(dplyr)"
dfp = transform(df_sobol,
    # :index => ByRow(string) => :index,
    # :par => ByRow(string) => :par
    :index => (x -> categorical(string.(x))) => :index,
    :par => (x -> categorical(string.(x))) => :par,
    :target => (x -> categorical(string.(x))) => :target)
recode!(dfp.target,
    "csom0" => "SOM",
    "csomD" => "ΔSOM",
    "csomDrel" => "ΔSOM/SOM",
    "r" => "ΔCUE")
levels!(dfp.target, ["ΔCUE"])    # index
levels!(dfp.index, ["total", "first_order"]);
# sort by decreasing effect withing index/target
sort!(dfp, [:index, :target, :value])
# relevel par so that index is preserved in plot
levels!(dfp.par,
    subset(dfp,
        AsTable([:target, :index]) => ByRow(x -> values(x) == ("ΔCUE", "total"))).par);

R"""
# https://stackoverflow.com/questions/12768176/unicode-characters-in-ggplot2-pdf-output
library(ggplot2)
library(dplyr)
cm2inch=1/2.54
fname=$(fname_plot)
#cairo_pdf(paste0("inst/22paper_upscaling/saltsens_",$(string(target)),".pdf"), width=8.3*cm2inch, 
cairo_pdf(fname, width=8.3*cm2inch, height=7*cm2inch, family="DejaVu Sans")
#png(paste0("tmp.png"))
p1 = $(dfp) %>%  
 #filter(index == "First Order") %>%  
 #filter(index == "Total") %>%  
 ggplot(aes(value, par, shape=index)) + geom_point() + facet_grid(.~target) +
 #xlab(paste0("Sobol sensitivity index of ",$(string(target)))) + 
 xlab(paste0("Sobol sensitivity indices")) + 
 scale_shape_discrete(solid = F) +
 scale_x_continuous(breaks=seq(0,1,by=0.5), limits=c(NA,1)) +
 theme_bw() +
  theme(legend.title = element_blank()) +
 #theme(axis.text.x = element_text(angle = 35, vjust = 0.8)) +
  #theme(legend.position = "bottom") +
  theme(legend.position = "none") +
  guides(linetype = guide_legend(nrow = 0.7, byrow = TRUE)) +
# theme(legend.position = c(0.05,0.05), legend.justification = c(0,0)) +
 theme(axis.title.y = element_blank()) 
print(p1)
dev.off()
"""
