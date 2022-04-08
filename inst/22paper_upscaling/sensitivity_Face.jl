# inspect sensitivity of increase in SOM stocks after 50yrs of increased
# C input

using Sesam
using ModelingToolkit, DifferentialEquations
using DataFrames, Tables
using Distributions
using Chain
using MTKHelpers

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
df_dist = rename!(DataFrame(columntable(parmsModeUpperRows)), collect(cols))


# names_opt = [
#     s.τ,
#     s.k_L,
#     s.k_R,
# ]
names_opt_all = df_dist.par 
names_opt = copy(names_opt_all)


ps = ProblemParSetter(sp, names_opt)
popt = get_paropt_labeled(ps, prob0)
prob = update_statepar(ps, popt, prob0)
# upd! = SystemParUpdater(collect(names_opt), sp)
# prob = deepcopy(prob0)
# popt0 = getindex.(Ref(p), names_opt) 
# popt = getindex.(Ref(p), names_opt) .* 1.2
# upd!(prob.p, popt)
sol_p = solve(prob);
#plot!(sol_p, vars=[s.L + s.R], tspan=(-20,tface))
som_change(sol_p) # SOM change different with updated parameters


#------ set up prior distributions of parameters and intervals of sens.analysis
using Random; Random.seed!(0);
using DistributionFits
# df_dist: moved up to access :par

#d1 = fit(LogNormal, @qp_m(8.), @qp_u(16.))
f1v = (dType, mode, upper) -> fit(dType, @qp_m(mode), @qp_uu(upper))
transform!(df_dist, Cols(:dType,:mode,:upper) => ByRow(f1v) => :dist)


"""
compute the values at quantiles ±δ_cp around x
with δ_cp difference in the cumulated probability.

A wider distribution prior distribution will result in a wider intervals.
"""
function calculate_parbounds(dist, x; δ_cp = 0.1 )
    cp_par = cdf(dist, x)
    cp_lower = max(0.005, cp_par - δ_cp)
    cp_upper = min(0.995, cp_par + δ_cp)
    qs = quantile.(dist, (cp_lower, cp_upper))
    (cp_par, cp_lower, cp_upper, qs...)
end
irow = 1
dist = df_dist[irow,:dist]
par = s.β_NB
x = p[par]
calculate_parbounds(dist, x)

#---- repeat with different δ_cp 
δ_cp = 0.1
δ_cp = 0.2
calculate_parbounds(dist, x; δ_cp)

f2v = (par, dist) -> begin
    ref = p[par]
    cp_par, cp_sens_lower, cp_sens_upper, sens_lower, sens_upper = calculate_parbounds(dist,ref; δ_cp)
    (;ref, sens_lower, sens_upper, cp_par, cp_sens_lower, cp_sens_upper)
end
transform!(df_dist, [:par,:dist,] => ByRow(f2v) => AsTable)
select(df_dist, :par, :mode, :ref, :sens_lower, :sens_upper, :upper)

i_outputLatex = () -> begin
    #produce latex strings for prior distributions in Table A1 in Wutzler22
    #using CategoricalArrays, Chain
    tmp = @chain df_dist begin
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
    N = 5001 # replicate
    #N = 10_000
    #N = 25_000
    #N = 50_000
    # N = 100_000
end

# ranges of cumulative probabilities given to R


#-------- second pass with fewer vars but larger N for variance
#names_omit = [s.k_mN_L, s.k_mN_L, s.l_N, s.ϵ_tvr, s.ν_N]
names_omit = []
names_opt = setdiff(names_opt_all, names_omit)
#upd! = SystemParUpdater(collect(names_opt), sp)
ps = ProblemParSetter(sp, names_opt)
popt0 = get_paropt_labeled(ps, prob0)
prob = update_statepar(ps, popt0, prob0)

#N = 10 # for testing setting y
#N = 50_000 * 8 # 50_000 took about 1/2 hour one night but sobolowen takes more

# transform :par
# df_dist_opt = @chain df_dist begin
#     select(:par) 
#     transform(:par => ByRow(strip_namespace ∘ symbol) => :par)
#     subset(:par => ByRow(x -> x in strip_namespace.(symbol.(names_opt))))
# end

# ranges of cumulative probabilities given to R
df_cfopt = @chain df_dist begin
    select(:par, :cp_sens_lower, :cp_sens_upper, :dist) 
    transform(:par => ByRow(strip_namespace ∘ symbol) => :par) 
    subset(:par => ByRow(x -> x in strip_namespace.(symbol.(names_opt))))
end

# here use sobolowen for better ci estimates    
cp_design = rcopy(R"""
library(sensitivity)
# for sobolowen X1,X2,X3 need to be data.frames, and need to convert
# design matrix (now also a data.frame) to array
set.seed(0815)
N = $(N)
δ_cp = $(δ_cp)
dfr <- $(select(df_cfopt, :par, :cp_sens_lower, :cp_sens_upper))
get_sample <- function(){
  setNames(data.frame(sapply(1:nrow(dfr), function(i){
    runif($(N), min = dfr$cp_sens_lower[i], max = dfr$cp_sens_upper[i])
  })), dfr$par)
}
#plot(density(get_sample()$k_L))
#lines(density(get_sample()$k_R))
#sensObject <- sobolSalt(NULL,get_sample(), get_sample(), nboot=100) 
sensObject <- soboltouati(NULL,get_sample(), get_sample(), nboot=100) # will be used down
# sobolowen returned fluctuating results on repeated sample matrices
#sensObject <- sobolowen(NULL,get_sample(), get_sample(), get_sample(), nboot=100) 
saveRDS(sensObject, paste0("sensObject_",N,"_",δ_cp*100,".rds"))
str(sensObject$X)
data.matrix(sensObject$X)
""");

R"""
δ_cp = $(δ_cp)
library(sensitivity)
sensObject = readRDS(paste0("sensObject_",N,"_",δ_cp*100,".rds"))
str(sensObject$X)
"""
cp_design
# cp_design=rcopy(R"""
# sensObject$X
# """)
size(cp_design) # nsample x npar
size(cp_design,1)*2/1000/3600 # hours one takes about 2ms

# transform cumulative probabilities back to quantiles
q_design = similar(cp_design)
for i_par in 1:size(cp_design,2)
    q_design[:,i_par] .= quantile.(df_cfopt.dist[i_par], cp_design[:,i_par])
end
#q_design

#---------- compute the response to changed parameters
using Feather
#nsamp = 10
nsamp = size(q_design,1)
#i_samp = 1
som_change0 = som_change(sol0)
som_changes = Array{Union{Float64,Missing}}(missing,length(som_change0), nsamp)
println("computing $nsamp samples.")
for i_samp in 1:nsamp
    popt = q_design[i_samp,:]
    #upd!(prob0, popt)
    prob = update_statepar(ps, popt, prob0)
    sol_p = solve(prob; saveat);
    #plot!(sol_p, vars=[s.L + s.R], tspan=(-20,tface))
    if sol_p.retcode == :Success 
        som_changes[:,i_samp] .= values(som_change(sol_p))
    end
    mod(i_samp,10_000) == 0 && print("$i_samp, ")
end
df_som_changes = DataFrame(transpose(som_changes), collect(keys(som_change0)))
fname = "df_som_changes_$(N)_$(Integer(δ_cp*100)).feather"
Feather.write(fname, df_som_changes)

describe(df_som_changes.csomD)
minimum(df_som_changes.csomD)

#-------- tell the results to sensitivity in R 
using RCall, Feather, DataFrames
df_som_changes = Feather.read(fname)
transform!(df_som_changes, [:csom0, :csomD] => ByRow((csom0,csomD) -> csomD/csom0) => :csomDrel)
names(df_som_changes)
target = :csomD  # sensitivity to what
target = :csom0  # sensitivity to what
#target = :csomDrel  # sensitivity to what
y = df_som_changes[!,target]
R"""
N = $(N)
library(sensitivity)
sensObject = readRDS(paste0("sensObject_",N,"_",δ_cp*100,".rds"))
#sensObject = readRDS(paste0("sensObject_",N,".rds"))
tell(sensObject, $(df_som_changes[!,target]))
"""

#--------------  extracting df_S and plotting 
df_S,df_T = rcopy(R"""
l <- list(sensObject$S, sensObject$T)
lapply(l, function(o){
    colnames(o) <- gsub(" ","", colnames(o)); o
})
""")

# rcopy(R"head(sensObject$y)")
# df_som_changes[1:6,target]
df_T[!,:par] .= df_cfopt.par
df_T[!,:order] .= :T
df_S[!,:par] .= df_cfopt.par
df_S[!,:order] .= :S
df_S[!,:originalT] .= df_T.original
sort!(df_T, :original, rev=true)
sort!(df_S, :originalT, rev=true) # sort same as df_T
par_ordered = df_T.par
using Plots
Plots.scatter(string.(df_T.par), df_T.original, label="Total", yerror=(df_T.original .- df_T.
min_c_i_,df_T.max_c_i_ .- df_T.original))
# must be the same sort order as df_T, otherwise x-lables change
Plots.scatter!(string.(df_S.par), df_S.original, label="First order", yerror=(df_S.original .- df_S.
min_c_i_, df_S.max_c_i_ .- df_S.original))

# df_ST = vcat(df_S, df_T)
# sort!(df_ST, :original, rev=true)
df_ST = sort!(vcat(select(df_S, Not(:originalT)), df_T), :original, rev=true)
df_ST[!,:target] .= target
# # both S1 and T1 with the same order
# df_S1 = subset(df_ST, :order => ByRow(==(:S)))
# df_T1 = subset(df_ST, :order => ByRow(==(:T)))

# scatter(string.(df_ST.par), df_ST.original, label="Total", yerror=(df_T.min_c_i_ .- df_T.original,df_T.max_c_i_ .- df_T.original))
# scatter!(string.(df_S.par), df_S.original, label="First order", yerror=(df_S.min_c_i_ .- df_S.original, df_S.max_c_i_ .- df_S.original))



# store the effects for this target
# feather does not support Symbols nor missings
using JLD2

fnameEffects = "df_ST_$(target)_$(N)_$(Integer(δ_cp*100)).jld2"
#fnameEffects = "df_ST_$(target)_$N.jld2"
jldsave(fnameEffects; df_ST)

#tmp = load_object(fnameEffects)

using CategoricalArrays

rcopy(R"""
sensObject$nboot
dim(sensObject$X)
""")

# plot the sensitivity
R"library(ggplot2)"
R"library(dplyr)"
df_STs = transform(df_ST, 
    # :order => ByRow(string) => :order,
    # :par => ByRow(string) => :par
    :order => (x -> categorical(string.(x))) => :order,
    :par => (x -> categorical(string.(x))) => :par,
    :target => (x -> categorical(string.(x))) => :target,
    )
levels!(df_STs.order, ["T","S"])    
recode!(df_STs.order, "T" => "Total", "S" => "First Order")    
levels!(df_STs.par, reverse(string.(par_ordered)))    # ordered decreasing
#levels!(df_STs.target, ["csomD", "csom0", "csomDrel"])    # order
levels!(df_STs.target, ["csom0", "csomD", "csomDrel"])    # order
recode!(df_STs.target, "csom0" => "SOM", "csomD" => "ΔSOM", "csomDrel" => "ΔSOM/SOM")    # also changes order
R"""
# https://stackoverflow.com/questions/12768176/unicode-characters-in-ggplot2-pdf-output
cm2inch=1/2.54
fname=paste0("inst/22paper_upscaling/saltsens_",$(string(target)),"_",$(Integer(δ_cp*100)),".pdf")
#cairo_pdf(paste0("inst/22paper_upscaling/saltsens_",$(string(target)),".pdf"), width=8.3*cm2inch, 
cairo_pdf(fname, width=8.3*cm2inch, height=7*cm2inch, family="DejaVu Sans")
#png(paste0("tmp.png"))
p1 = $(df_STs) %>%  
 #filter(order == "First Order") %>%  
 #filter(order == "Total") %>%  
 ggplot(aes(original, par, shape=order)) + geom_point() + facet_grid(.~target) +
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

#----------- plot all sensitivities together
#df_ST = map([:csom0, :csomD, :csomDrel]) do target
#do not plot the csomDrel - its not very different from csomD
df_ST = map([:csom0, :csomD]) do target
    fnameEffects = "df_ST_$(target)_$(N)_$(Integer(δ_cp*100)).jld2"
    #fnameEffects = "df_ST_$(target)_$N.jld2"
    df_STt = load_object(fnameEffects)
end
df_ST = vcat(df_ST...)
target = :combined
# invoke ggplot  above

# order parameters by total effects of target csomD
using Chain
par_ordered = @chain df_ST begin
    subset(:target => ByRow(==(:csomD)))
    subset(:order => ByRow(==(:T)))
    sort(:original, rev=true)
    getproperty(:par)
end;

# 3D plot of 3 most important parameters
# is = sample(1:N, 1000)
# qs = DataFrame(q_design[is,:], collect(parsymbol.(names_opt)))
# y = df_som_changes.csomD[is]
# qs.δSOM = y


# using CairoMakie
# const CM = CairoMakie
# CM.scatter(qs.β_NB, df_som_changes.com0, color=qs.β_NEnz)
is = sample(1:N, 2000)
df_q = DataFrame(q_design, collect(strip_namespace.(paroptsyms(ps))))
df_qs = DataFrame(q_design[is,:], collect(strip_namespace.(paroptsyms(ps))))
df_s = df_som_changes[is,:]
scatter(df_qs.k_R, df_s.csom0, zcolor=df_qs.ϵ_tvr, xlab="k_R (1/yr)", ylab="SOM stocks (g/m2)", label=nothing, colorbar_title = " \nϵ_tvr (g/g)", right_margin = 6Plots.mm)


scatter(df_qs.β_NB, df_s.csomD, zcolor=df_qs.ϵ_tvr, xlab="β_NB (g/g)", ylab="Δ SOM stocks (g/m2)", label=nothing, colorbar_title = " \nϵ_tvr (g/g)", right_margin = 6Plots.mm)



scatter(df_qs.β_NB, df_s.B0)
scatter(df_qs.k_L, df_s.csom0)
 #
@chain df_ST begin
    subset(:target => ByRow(==(:csomD)))
    subset(:order => ByRow(==(:T)))
end

i_inspect_Ndeposition = () -> begin
    minb, maxb = extrema(df_q.β_NB)
    x = range(minb, stop=maxb, length=12)
    ps_β_NB = ProblemParSetter(sp, (s.β_NB,))
    prob_ndep = prob0 # without N deposition
    ndep = 0.7 # 7 kg/ha/yr Hueso11 7 to 12 kg/ha/yr 
    prob_ndep = update_statepar(ProblemParSetter(sp, (pl.i_IN0,)),(ndep,), prob0)
    xi = x[3]
    solps = map(x) do xi
        probp = update_statepar(ps_β_NB, (xi,), prob_ndep)
        solp = solve(probp);
    end
    #scs0 = DataFrame(som_change.(solps)); scs0[!,:ndep] .= 0.0
    scs = DataFrame(som_change.(solps)); scs0[!,:ndep] .= ndep

    label_par(ps,prob_ndep.p)

    plot(x, scs.csom0)
    plot!(x, scs0.csom0, label = "ndep=0")
    # not a significant change of dynamics, only I levels are different

    plot(x, scs.csomD)
    plot!(x, scs0.csomD, label = "ndep=0")
    # slightly shifted but same dynamics
end

i_inspect_high_csom0 = () -> begin
    # what parameters yield the high initial stocks with low CN_B?
    using LinearAlgebra
    popt1 = get_paropt_labeled(ps, update_statepar(psx, (x[1],), prob_ndep))
    i_close = argmin(mapslices(pi -> norm(pi .- popt1), q_design; dims=2)[:,1])
    popt = q_design[i_close,:]
    df_tmp = DataFrame(transpose(hcat(popt0s, popt)), parsymbol.(names_opt))
    df_som_changes[i_close,:] # same low csom0

    # compare original solution with update β_NB to low value to (popt2)
    # parameter with low β_NB (popt1) that yields a high csom0
    i_low_β_BN = findall(df_q.β_NB .<= 10)
    df_som_changes[i_low_β_BN,:csom0]
    im = i_low_β_BN[argmax(df_som_changes[i_low_β_BN,:csom0])]
    df_som_changes[im,:]
    popt1 = df_q[im,:]
    #df_qs_low = subset(df_qs, :β_NB => ByRow(<=(10)))
    probp = update_statepar(ps, collect(popt1), prob0)
    solp = solve(probp)
    som_change(solp)

    ps_β_NB = ProblemParSetter(sp, (s.β_NB,))
    probp2 = update_statepar(ps_β_NB, (9.94,), prob0)
    popt2 = get_paropt_labeled(ps, probp2)
    label_par(ps,probp2.p)
    solp2 = solve(probp2)
    som_change(solp2)

    # which parameters are very different?
    rename!(DataFrame(vcat(Tuple(popt1), Tuple(popt2))), names(df_qs))
    # ------ modify k_R
    probp3 = update_statepar( ProblemParSetter(sp, (s.k_R,)), (popt1.k_R,), probp2)
    solp3 = solve(probp3)
    som_change(solp3)

    plot(solp)
    tspanp = (-450,0)
    # I_N similar
    plot(solp, tspan=tspanp, vars=[s.I_N])
    plot!(solp2, tspan=tspanp, vars=[s.I_N])
    # lower biomass with high csom0
    plot(solp, tspan=tspanp, vars=[s.B], label="high coms0")
    plot!(solp2, tspan=tspanp, vars=[s.B])
    plot!(solp2, tspan=tspanp, vars=[s.B])
    plot!(solp3, tspan=tspanp, vars=[s.B])
end

i_vary_single_parameters = () -> begin
    par_v = :k_R
    par_v = :β_NB
    par_v = :τ
    par_v = :κ_E
    minb, maxb = extrema(df_q[!,par_v])
    x = range(minb, stop=maxb, length=12)
    ps_v = ProblemParSetter(sp, (getproperty(s, par_v),))
    xi = x[3]
    solps = map(x) do xi
        probp = update_statepar(ps_v, (xi,), prob0)
        solp = solve(probp);
    end
    scs = DataFrame(som_change.(solps))
    scs[!,par_v] = x
    # using StatsPlots #@df
    # @df scs plot(:k_R, [:csom0])
    plot(x, scs.csom0, xlab=string(par_v), ylab="csom0", label=nothing)
end



