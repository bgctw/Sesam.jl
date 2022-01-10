using ModelingToolkit, DifferentialEquations, IfElse

@parameters t ϵ_tvr κ_E a_E m τ k_L k_R k_mN ϵ
@parameters α_R0
@parameters i_L0
@variables B(t) L(t) R(t) α_R(t) i_L(t) cumresp(t)
@variables syn_Enz(t) r_M(t) tvr_B(t) dec_LPot(t) dec_L(t) dec_RPot(t) dec_R(t) u_C(t)
@variables C_synBCt(t) C_synBC(t) r_G(t) r_tvr(t) syn_B(t)
@variables C_tot(t) r_B(t) r_GEnz(t) r_O(t)
@variables α_L(t)

D = Differential(t)

eqs = [
    D(B) ~ syn_B - tvr_B,
    D(L) ~ -dec_L + i_L,
    D(R) ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
    D(α_R) ~ 0.0,
    syn_Enz ~ a_E*B,
    r_M ~ m*B,
    tvr_B ~ τ*B,
    dec_LPot ~ k_L * L,
    dec_L ~ dec_LPot*(α_L * syn_Enz)/(k_mN + α_L*syn_Enz),
    dec_RPot ~ k_R * R,
    dec_R ~ dec_RPot*(α_R * syn_Enz)/(k_mN + α_R*syn_Enz),
    u_C ~ dec_L + dec_R + κ_E*syn_Enz,
    C_synBCt ~ u_C - syn_Enz/ϵ - r_M,
    C_synBC ~ IfElse.ifelse(C_synBCt > 0.0, ϵ*C_synBCt, C_synBCt), 
    syn_B ~ C_synBC, # TODO add more limitations
    α_L ~ (1-α_R),
    # auxilaries
    D(cumresp) ~ r_B + r_tvr,
    r_tvr ~ (1-ϵ_tvr)*tvr_B,
    r_B ~ r_GEnz + r_G + r_M + r_O,
    r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
    r_G ~ IfElse.ifelse(syn_B > 0, (1-ϵ)/ϵ * syn_B, 0.0), 
    r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
    i_L ~ i_L0 # TODO specify as Pin and compose with different system
    ]

p = [
    ϵ_tvr => 0.45,   # carbon use efficiency of microbial tvr (part by predators 
        #which respire and corresponding amount of N must be mineralized)
    κ_E => 0.8,     ##<< amount of recycling enzyme turnover by biomass (
        # added to assimilable, i.e. uptake instead of R)
    a_E => 0.001*365,   ##<< C biomass allocated to enzymes 1/day /microbial biomass 
    m => 0.005*365,    ##<< maintenance respiration rate   1/day /microbial biomass,    
    τ => 1/60*365,  ##<< biomass turnover rate (12 days)    
    k_L => 5.0,       ##<< 1/(x years)   # formerly 1 year
    k_R => 1/(20.0),        ##<< 1/(x years) # to demonstrate changes on short time scale
    k_mN => 0.05 * 60, # enzyme half-saturation constant, in magnitude of enzymes * 
        # /yr enzyme turnover 60 times a year
    ϵ => 0.5,      ##<< carbon use efficiency for growth respiration
    i_L0 => 400.0,         # g/m2 input per year (half NPP)
    #i_L => t -> 1 - exp(-t),  # litter input
]    

@named sesam1 = ODESystem(eqs)    
sesam1 = structural_simplify(sesam1)
states(sesam1)
equations(sesam1)

u0 = [
    B => 17,
    L => 100,
    R => 1100,
    cumresp => 0.0,
    α_R => 0.3, # enzyme synthesis into L # TODO model by optimality
    ]

tspan = (0.0,100.0)    

prob = ODEProblem(sesam1,u0, tspan, p)
#prob = ODEProblem(sesam1,u0, tspan, p, jac=true)
sol = solve(prob, Tsit5())

using Plots
plot(sol, vars=[L,B])
# increase of C in system is the same as input

Dict(p)[i_L]
# total C stock equals input
#plot(sol, vars=[R + L + B + cumresp])
#plot!(t -> Dict(p)[i_L]*t + sum(getindex.(Ref(Dict(u0)),[B,L,R])))
# uptake equals usage of C
# plot(sol, tspan=(0.0,2.0), vars=[u_C, r_B+syn_B+syn_Enz])









