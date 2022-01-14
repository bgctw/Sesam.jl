# inspect sensitivity of increase in SOM stocks after 50yrs of increased
# C input

using Sesam
using ModelingToolkit, DifferentialEquations
using DataFrames, Tables

tspinup = 500.0; tface=100.0
#@named s = sesam3(use_seam_revenue=true)
@named s = sesam3(use_seam_revenue=false)
@named pl = plant_face(t1=0.0,t2=tface)
@named sp = plant_sesam_system(s,pl)

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
    s.k_mN => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #i_L => t -> 1 - exp(-t),  # litter input
    pl.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    #pl.β_Ni0 => 25.0,
    pl.β_Ni0 => 30.0,
    pl.i_IN0 => 0.0,   ##<< input of mineral N,
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
saveat=[0.0,tface] # save at start of increase and 50yrs later
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
    csom0 = csom[findlast(==(0.0),sol.t)]
    csomF = csom[findfirst(==(tface),sol.t)]
    csomD = csomF - csom0
    (;csom0, csomF, csomD)
end
som_change(sol0)

"Get the symbolic represenation of the Num omitting the namespaces."
parsymbol(num) = Symbol(first((match(r"₊(.+)$", string(num))).captures))
parsymbol(s.k_L)
parsyssymbol(num) = Symbol(replace(string(num), "₊" => "."))
parsyssymbol(s.k_L)

cols = (:par, :dType, :mode, :upper)
parmsModeUpperRows = [
     (s.β_NB, LogNormal, 8. , 16.),
     (s.β_NEnz, LogNormal, 3.0 , 3.5),
     (s.k_mN, LogNormal, 60*0.05 , 120*2.),
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


names_opt = [
    s.τ,
    s.k_L,
    s.k_R,
]
names_opt_all = df_dist.par 
names_opt = copy(names_opt_all)


upd! = SystemParUpdater(collect(names_opt), sp)
prob = deepcopy(prob0)
popt = getindex.(Ref(p), names_opt) .* 1.2
upd!(prob.p, popt)
sol_p = solve(prob) 
#plot!(sol_p, vars=[s.L + s.R], tspan=(-20,tface))
som_change(sol_p) # SOM change different with updated parameters


#------ set up prior distributions of parameters and intervals of sens.analysis
using Random, Distributions; Random.seed!(0);
using DistributionFits
# df_dist: moved up to access :par

#d1 = fit(LogNormal, @qp_m(8.), @qp_u(16.))
f1 = (dType, mode, upper) -> fit(dType, @qp_m(mode), @qp_uu(upper))
transform!(df_dist, [:dType,:mode,:upper] => ByRow(f1) => :dist)


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

f2 = (par, dist) -> begin
    ref = p[par]
    cp_par, cp_sens_lower, cp_sens_upper, sens_lower, sens_upper = calculate_parbounds(dist,ref)
    (;ref, sens_lower, sens_upper, cp_par, cp_sens_lower, cp_sens_upper)
end
transform!(df_dist, [:par,:dist,] => ByRow(f2) => AsTable)
select(df_dist, :par, :mode, :ref, :sens_lower, :sens_upper, :upper)

#------------ get the design matrix from sensitivity package in R
using RCall
using Pipe
if !@isdefined N 
    N = 100
    N = 2000
    N = 10_000
    N = 50_000
    # N = 100_000
end

# ranges of cumulative probabilities given to R
df_cfopt = @pipe df_dist |> 
    select(_,:par, :cp_sens_lower, :cp_sens_upper, :dist) |> 
    transform(_, :par => ByRow(parsymbol) => :par) |>
    subset(_, :par => ByRow(x -> x in parsymbol.(names_opt)))

cp_design = rcopy(R"""
library(sensitivity)
set.seed(0815)
dfr <- $(select(df_cfopt, :par, :cp_sens_lower, :cp_sens_upper))
get_sample <- function(){
  sapply(1:nrow(dfr), function(i){
    runif($(N), min = dfr$cp_sens_lower[i], max = dfr$cp_sens_upper[i])
  })
}
#plot(density(get_sample()$k_L))
#lines(density(get_sample()$k_R))
#sensObject <- sobolSalt(NULL,get_sample(), get_sample(), nboot=100) 
sensObject <- soboltouati(NULL,get_sample(), get_sample(), nboot=100) # will be used down
sensObject$X
""");
size(cp_design) # nsample x npar

# transform cumulative probabilities back to quantiles
q_design = similar(cp_design)
for i_par in 1:size(cp_design,2)
    q_design[:,i_par] .= quantile.(df_cfopt.dist[i_par], cp_design[:,i_par])
end
#q_design

#---------- compute the response to changed parameters
nsamp = 10
nsamp = size(q_design,1)
#i_samp = 1
som_change0 = som_change(sol0)
som_changes = Array{Union{Float64,Missing}}(missing,length(som_change0), nsamp)
println("computing $nsamp samples.")
for i_samp in 1:nsamp
    popt = q_design[i_samp,:]
    upd!(prob.p, popt)
    sol_p = solve(prob; saveat) 
    #plot!(sol_p, vars=[s.L + s.R], tspan=(-20,tface))
    if sol_p.retcode == :Success 
        som_changes[:,i_samp] .= values(som_change(sol_p))
    end
    mod(i_samp,1000) == 0 && print("$i_samp, ")
end
df_som_changes = DataFrame(transpose(som_changes), collect(keys(som_change0)))
describe(df_som_changes.csomD)

#-------- tell the results to sensitivity in R 
@rput df_som_changes;
R"""
tell(sensObject, df_som_changes$csomD)
"""

R"""
rownames(sensObject$S) <- rownames(sensObject$T) <- dfr$par
#print(sensObject)
pdf(paste0("inst/22paper_upscaling/sobol_N=", $N, ".pdf"), width=9, height=5)
if (inherits(sensObject,"sobolSalt")) {
    plot(sensObject, choice = 1)
    plot(sensObject, choice = 2)
} else {
    plot(sensObject)
}
dev.off()
sensObject
"""


R"""
y <- rnorm(88000)
tell(sensObject,y)
x$T
"""

R"""
str(sensObject)
"""

R"""
sensObject$nboot
"""

df_S,df_T = rcopy(R"""
l <- list(sensObject$S, sensObject$T)
lapply(l, function(o){
    colnames(o) <- gsub(" ","", colnames(o)); o
})
""")

df_S[!,:par] = df_cfopt.par
sort!(df_S, :original, rev=true)
df_T[!,:par] = df_cfopt.par
sort!(df_T, :original, rev=true)
scatter(string.(df_T.par), df_T.original, label="Total", yerror=(df_T.min_c_i_,df_T.max_c_i_))
scatter!(string.(df_S.par), df_S.original, label="First order", yerror=(df_S.min_c_i_,df_S.max_c_i_))




#-------- second pass with fewer vars but larger N for variance
names_omit = [s.k_mN, s.l_N, s.ϵ_tvr, s.ν_N]
names_opt = setdiff(names_opt_all, names_omit)
upd! = SystemParUpdater(collect(names_opt), sp)


N = 10 # for testing setting y
N = 50_000 * 8 # 50_000 took about 1/2 hour one night but sobolowen takes more

# ranges of cumulative probabilities given to R
df_cfopt = @pipe df_dist |> 
    select(_,:par, :cp_sens_lower, :cp_sens_upper, :dist) |> 
    transform(_, :par => ByRow(parsymbol) => :par) |>
    subset(_, :par => ByRow(x -> x in parsymbol.(names_opt)))

# here use sobolowen for better ci estimates    
cp_design = rcopy(R"""
library(sensitivity)
# for sobolowen X1,X2,X3 need to be data.frames, and need to convert
# design matrix (now also a data.frame) to array
set.seed(0815)
dfr <- $(select(df_cfopt, :par, :cp_sens_lower, :cp_sens_upper))
get_sample <- function(){
  data.frame(sapply(1:nrow(dfr), function(i){
    runif($(N), min = dfr$cp_sens_lower[i], max = dfr$cp_sens_upper[i])
  }))
}
#plot(density(get_sample()$k_L))
#lines(density(get_sample()$k_R))
#sensObject <- sobolSalt(NULL,get_sample(), get_sample(), nboot=100) 
#sensObject <- soboltouati(NULL,get_sample(), get_sample(), nboot=100) # will be used down
sensObject <- sobolowen(NULL,get_sample(), get_sample(), get_sample(), nboot=500) 
data.matrix(sensObject$X)
""");
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
nsamp = 10
nsamp = size(q_design,1)
#i_samp = 1
som_change0 = som_change(sol0)
som_changes = Array{Union{Float64,Missing}}(missing,length(som_change0), nsamp)
println("computing $nsamp samples.")
for i_samp in 1:nsamp
    popt = q_design[i_samp,:]
    upd!(prob.p, popt)
    sol_p = solve(prob; saveat) 
    #plot!(sol_p, vars=[s.L + s.R], tspan=(-20,tface))
    if sol_p.retcode == :Success 
        som_changes[:,i_samp] .= values(som_change(sol_p))
    end
    mod(i_samp,10_000) == 0 && print("$i_samp, ")
end
df_som_changes = DataFrame(transpose(som_changes), collect(keys(som_change0)))
describe(df_som_changes.csomD)

#-------- tell the results to sensitivity in R 
@rput df_som_changes;
R"""
tell(sensObject, df_som_changes$csomD)
"""

# extracting df_S and plotting as above
