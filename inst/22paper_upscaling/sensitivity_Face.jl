# inspect sensitivity of increase in SOM stocks after 50yrs of increased
# C input

using Sesam
#push!(LOAD_PATH, expanduser("~/julia/scimltools/")) # access local package repo
using ModelingToolkit, DifferentialEquations
using DataFrames, Tables
using Distributions
using Chain
using MTKHelpers
import ComponentArrays as CA
import SubglobalSensitivityAnalysis as SSA

tspinup = 500.0; tface=100.0
#@named s = sesam3(use_seam_revenue=true)
@named s = sesam3_revMM(use_seam_revenue=false)
@named pl = plant_face(t1=0.0,t2=tface)
@named sp = plant_sesam_system(s,pl)
#equations(sp)

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
    s.l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    #s.l_N => 0.0,       
    s.ν_N =>  0.9,     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    # for N uptake: take the defaults which take as much N as supplied by litter
    # pl.u_PlantNmax0 => Inf32, # only constrained by k_PlantN in min function
    # pl.k_PlantN0 => 10.57, #0.0289652*365     ##<< plant uptake rate first order of IN
)
p = merge(pC, pN)

u0 = u0C = Dict(
    s.B => 17,
    #s.B => 34,
    s.L => 100,
    #s.L => 110,
    s.R => 1100,
    #s.R => 3250,
    #s.cumresp => 0.0,
    s.α_R => 0.5, 
)
u0C[s.α_L] = 1.0 - u0C[s.α_R]
u0N = Dict(
    s.I_N => 0.04, ##<< inorganic pool gN/m2 
    s.L_N => u0[s.L]/p[pl.β_Ni0],
    #s.R_N => u0[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
    s.R_N => u0[s.R]/7.0
    )
u0 = merge(u0C, u0N)    


tspan_sim = (-tspinup,tface) # simulate 500 yrs spinup, increase at yr 20
saveat=[0.0,tface] # save at start of increase and tface years later
prob0 = ODEProblem(sp, u0, tspan_sim, p) #2ms
#prob0 = ODEProblem(sp, u0, tspan_sim, p; jac=true) # actually slower: 4ms
sol0t = sol = solve(prob0)
#plot(sol0t, vars=[s.L + s.R], tspan=(-20,tface))
sol0 = solve(prob0; saveat)

# using BenchmarkTools
# @btime sol0 = solve(prob0; saveat)

# sensitiviy measure is the change in SOM stocks after 50 yrs increased input
# compared to equilibrium
function som_change(sol) 
    !(sol.retcode ∈ (:Success, :Terminated)) && return(
        (;csom0=missing, csomF=missing, csomD=missing, 
        I_N0 = missing, I_NF = missing, B0 = missing))
    csom = sol[s.L + s.R]
    I_N = sol[s.I_N]
    i0 = findfirst(>=(0.0),sol.t)
    iF = findlast(<=(tface),sol.t)
    csom0 = csom[i0]
    csomF = csom[iF]
    csomD = csomF - csom0
    (;csom0, csomF, csomD, I_N0 = I_N[i0], I_NF = I_N[iF], B0 = sol[s.B][i0])
end
som_change(sol0)


# "Get the symbolic represenation of the Num omitting the namespaces."
# parsymbol(num) = Symbol(first((match(r"₊(.+)$", string(num))).captures))
# parsymbol(s.k_L)
# parsyssymbol(num) = Symbol(replace(string(num), "₊" => "."))
# parsyssymbol(s.k_L)

cols = (:par, :dType, :mode, :upper)
parmsModeUpperRows = [
     (s.β_NB, LogNormal, 8. , 16.),
     (s.β_NEnz, LogNormal, 3.0 , 3.5),
     (s.k_mN_L, LogNormal, 60*0.05 , 120*2.),
     (s.k_mN_R, LogNormal, 60*0.05 , 120*2.),
     (s.κ_E, LogNormal, 0.7 , 0.9),
     (s.k_R, LogNormal, 1/100, 1/10),
     (s.k_L, LogNormal, 12/18, 12/4), # tvr time between 18 months and 4 months
     (s.a_E, LogNormal, 0.001*365 , 0.005*365),
     (s.m, LogNormal, 0.005*365 , 0.02*365),
     (s.τ, LogNormal, 1.0/60*365 , 1.0/5*365),
     (s.ϵ, LogitNormal, 0.5 , 0.7),
     (s.ϵ_tvr, LogitNormal, 0.3 , 0.9),
     (s.i_BN, LogNormal, 0.4 , 4),
     (s.l_N, LogNormal, 0.96 , 5*0.96),
     (s.ν_N, LogNormal, 0.9 , 0.99),
     #(:kIPlant, LogNormal, 10.57 , 20)
]
df_dist = SSA.fit_distributions(parmsModeUpperRows)
transform!(df_dist, :par => identity => :par_num)
transform!(df_dist, :par => ByRow(symbol) => :par)

names_opt_all = df_dist.par 
names_omit = []
names_opt = setdiff(names_opt_all, names_omit)
#upd! = SystemParUpdater(collect(names_opt), sp)
pset = ProblemParSetter(sp, CA.Axis(symbol.(names_opt)))
popt0 = get_paropt_labeled(pset, prob0)

function sim_face(popt)
    prob = update_statepar(pset, popt, prob0)
    solve(prob; saveat)
end
sol_p = sim_face(popt0)
som_change(sol_p) # SOM change different with updated parameters


#---- repeat with different δ_cp 
δ_cp = 0.2
δ_cp = 0.1

df_dist_opt = subset(df_dist, :par => ByRow(x -> x ∈ names_opt))
SSA.set_reference_parameters!(df_dist_opt, Dict(propertynames(popt0) .=> values(popt0)))

i_outputLatex = () -> begin
    #produce latex strings for prior distributions in Table A1 in Wutzler22
    #using CategoricalArrays, Chain
    tmp = @chain df_dist_opt begin
        # prescribe sort order of :par to match paper
        # for this create ordered Categorical vector and set levels
        transform(:par => ByRow(string ∘ parsymbol) => :pars) 
        transform(:pars => (x -> CategoricalArray(x,ordered=true)) => :pars)
        # transform(:pars => (x -> levels!(x,["β_NB","β_NEnz","k_R","k_L","κ_E","a_E","k_mN","τ","m","ϵ","ϵ_tvr","ν_N","i_BN","l_N"])) => :pars)
        transform(:pars => (x -> levels!(x,["β_NB","β_NEnz","k_R","k_L","κ_E","a_E","k_mN_R","k_mN_L","τ","m","ϵ","ϵ_tvr","ν_N","i_BN","l_N"])) => :pars)
        sort(:pars)
        # getproperty(:pars)
        # join("\",\"")
        transform(:dist => ByRow(d -> (Base.typename(typeof(d)).name, d.μ,d.σ,quantile(d,0.025))) => [:typename, :μ, :σ, :lower]) 
        transform([:typename,:μ,:σ] => ByRow((n,m,s) ->  "$n($(round(m,digits=2)),$(round(s,sigdigits=2)))") => :typestr)
        select(:pars,:typestr,:lower, :mode,:upper)
    end
    # format latex table
    tmp2 = @chain tmp begin
        transform([:lower, :mode,:upper] => ByRow((x...) -> round.(x, sigdigits=2)) => [:lower, :mode,:upper])
        transform([:typestr,:lower, :mode,:upper] => ByRow((x...) -> join(x," & ")) => :latex)
        getproperty(:latex)      
    end
    print(join(tmp2,"\n"))
end


#------------ get the design matrix from sensitivity package in R
using RCall
if !@isdefined N 
    N = 100
    #N = 5000
    N = 5000 # replicate
    #N = 10_000
    #N = 25_000
    #N = 50_000
    # N = 100_000
end

scen_str = "N$(N)_deltacp$(Integer(δ_cp*100))"
fname_sens = joinpath("tmp","sens_face_$(scen_str).rds")
fname_output = joinpath("tmp","sens_face_$(scen_str).feather")
fname_plot = joinpath("tmp","sens_face_$(scen_str).pdf")

estim_file = SSA.SobolTouati(;rest=SSA.RSobolEstimator("sens_touati", fname_sens))

SSA.calculate_parbounds!(df_dist_opt)
X1 = SSA.get_uniform_cp_sample(df_dist_opt, N);
X2 = SSA.get_uniform_cp_sample(df_dist_opt, N);
cp_design = SSA.generate_design_matrix(estim_file, X1, X2);
#SSA.reload_design_matrix(estim_file);
q_design = SSA.transform_cp_design_to_quantiles(df_dist_opt, cp_design);

fsens = (popt) -> begin
    local sol_p = sim_face(popt)
    som_change(sol_p) # SOM change different with updated parameters
end
fsens(first(eachrow(q_design)))

res = map(r -> fsens(r), eachrow(q_design));
describe(DataFrame(res))
import Feather
Feather.write(fname_output, res)

df_sobol = vcat(map(propertynames(res[1])) do target
    y = [tup[target] for tup in res]
    df_sobol =  SSA.estimate_sobol_indices(estim_file, y, df_dist_opt.par)
    transform!(df_sobol, [] => ByRow(() -> target) => :target)
end...)

#-------------- plot the results
using CategoricalArrays

# plot the sensitivity
using RCall
R"library(ggplot2)"
R"library(dplyr)"
dfp = transform(df_sobol, 
    # :index => ByRow(string) => :index,
    # :par => ByRow(string) => :par
    :index => (x -> categorical(string.(x))) => :index,
    :par => (x -> categorical(string.(x))) => :par,
    :target => (x -> categorical(string.(x))) => :target,
    )
recode!(dfp.target, "csom0" => "SOM", "csomD" => "ΔSOM", "csomDrel" => "ΔSOM/SOM", "r" => "ΔCUE")    
subset!(dfp, :target => ByRow(x -> x ∈ ["SOM","ΔSOM"]))    
levels!(dfp.target, ["SOM","ΔSOM"])    
levels!(dfp.index, ["total","first_order"]);    
# sort by decreasing effect withing index/target
sort!(dfp, [:index, :target, :value])   
# relevel par so that index is preserved in plot
levels!(dfp.par, subset(dfp, AsTable([:target,:index]) => ByRow(
     x -> values(x) == ("ΔSOM","total"))).par);

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

