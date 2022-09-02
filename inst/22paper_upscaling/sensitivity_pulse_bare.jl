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

@named s = sesam3_revMM(use_seam_revenue=false)
@named pl = plant_const()
@named sp = plant_sesam_system(s,pl)

p = pC = Dict(
    s.ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
        #which respire and corresponding amount of N must be mineralized)
    s.κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
        # added to assimilable, i.e. uptake instead of R)
    s.a_E => 0.001*365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    s.m => 1.825, #0.005*365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    s.τ => 1/60*365,  ##<< biomass turnover rate (12 days)    
    s.k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    s.k_R => 1/(20.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    s.k_mN_L => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    s.k_mN_R => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
    # /yr enzyme turnover 60 times a year
    s.ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    #
    pl.i_L0 => 0,
)
pN = Dict(
    s.i_BN => 0.4, ##<< potential immobilization flux rate 
    s.β_NEnz => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    s.β_NB => 11.0,
    s.l_N => 0, #0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    #s.l_N => 0.0,       
    s.ν_N =>  0.9,     # microbial N use efficiency accounting for apparent 
    ## minceralization of N during uptake in heterogeneous soils
    # for N uptake: take the defaults which take as much N as supplied by litter
    pl.u_PlantNmax0 => 0.0, 
    pl.k_PlantN0 => 0.0, 
    pl.i_IN0 => 0.0,
    pl.i_IP0 => 0.0,
    #pl.β_Ni0, # need to be set dependend on scenario
    pl.β_Pi0 => NaN,
)
p = merge(pC, pN)

u0 = u0C = Dict(
    s.B => 0.3,
    #s.B => 34,
    s.L => 100,
    #s.L => 110,
    s.R => 1200,
    #s.R => 3250,
    #s.cumresp => 0.0,
    s.α_R => 0.5, 
)
u0C[s.α_L] = 1.0 - u0C[s.α_R]

function get_u0(cn0, u0C, p)
    u0N = Dict(
        s.I_N => 1, ##<< inorganic pool gN/m2 
        s.L_N => u0C[s.L]/cn0,
        #s.R_N => u0[s.R]/calculate_β_NR_sesam3(p,s) #p[s.β_NB],
        s.R_N => u0C[s.R]/p[s.β_NB],
        pl.β_Ni0 => 30.0 #cn0,
        )
    u0 = merge(u0C, u0N)  
end  
u0_50 = get_u0(50, u0C, p)

function simExp(u0_cn0)
    tend = 4.0
    tspan_sim = (0,tend) # simulate 500 yrs spinup, increase at yr 20
    saveat=[0.0,tend] # save at start of increase and tface years later
    prob0 = ODEProblem(sp, u0_cn0, tspan_sim, p) #2ms
    #prob0 = ODEProblem(sp, u0, tspan_sim, p; jac=true) # actually slower: 4ms
    sol0t = sol = solve(prob0, );
end
sol50 = simExp(u0_50);

function tmp_fig4_pulse50_90()
    tmp = map([50,70,90]) do cn0
        u0_cn0 = get_u0(cn0, u0C, p)    
        cn0, u0_cn0
    end
    df = rename!(DataFrame(tmp), [:cn0, :u0])

    #tmp = simExp(df.u0[1])
    transform!(df, :u0 => ByRow(simExp) => :sol)
    function tmpf()
        sol0t = df.sol[1];
        plot(sol0t, idxs=[s.L + s.R])
        plot(sol0t, idxs=[s.syn_B])
        plot!(df.sol[2], idxs=[s.syn_B])
        pl = plot(df.sol[1], idxs=[s.syn_B/s.u_C])
        for i = 2:nrow(df) 
            print(i)
            plot!(pl, df.sol[i], idxs=[s.syn_B/s.u_C])
        end
        display(pl)
        plot(sol0t, idxs=[s.lim_C, s.lim_N])
    end
    sol0 = solve(prob0; saveat)
end

CUEmin = minimum(sol50[s.syn_B/s.u_C])

