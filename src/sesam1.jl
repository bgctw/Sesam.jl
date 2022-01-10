using ModelingToolkit, DifferentialEquations, IfElse

@parameters t 
D = Differential(t)

@parameters ϵ_tvr κ_E a_E m τ k_L k_R k_mN ϵ
@parameters α_R0
@parameters i_L0
@variables B(t) L(t) R(t) α_R(t) i_L(t) cumresp(t)
@variables syn_Enz(t) r_M(t) tvr_B(t) dec_LPot(t) dec_L(t) dec_RPot(t) dec_R(t) u_C(t)
@variables C_synBCt(t) C_synBC(t) r_G(t) r_tvr(t) syn_B(t)
@variables C_tot(t) r_B(t) r_GEnz(t) r_O(t)
@variables α_L(t)

@variables L_N(t) R_N(t) I_N(t)
@variables i_IN(t) Φ_N(t) Φ_Nu(t) Φ_NB(t) u_PlantN(t) u_NOM(t) u_immNPot(t) 
@variables u_N(t) N_synBN(t) M_ImbN(t) β_NL(t) β_NR(t) β_Ni(t)
@parameters β_NE β_NB β_Ni0 u_PlantNmax k_PlantN l_N ν_N i_BN I_IN0


p = Dict(
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
    #
    # ----------- N -------------
    β_NE => 3.1,     # Sterner02: Protein (Fig. 2.2.), high N investment (low P) need 
    β_Ni0 => 25,
    β_NB => 11.0,
    u_PlantNmax => Inf32, # only constrained by k_PlantN in min function
    k_PlantN => 10.57, #0.0289652*365     ##<< plant uptake rate first order of IN
    l_N => 0.96,       #0.00262647*365     ##<< leaching rate of mineralN lN IN
    ν_N =>  0.9,     # microbial P use efficiency accounting for apparent 
    ## minceralization of P during uptake in heterogeneous soils
    i_BN => 0.4, ##<< potential immobilization flux rate 
    I_IN0 => 0,   ##<< input of mineral N,
)
# constrain to experiment
p[k_PlantN] = 1000.0
p[u_PlantNmax] = p[i_L0]/p[β_Ni0]
p[l_N] = 0.0

eqs = [
    D(B) ~ syn_B - tvr_B,
    D(L) ~ -dec_L + i_L,
    D(R) ~ -dec_R + ϵ_tvr*tvr_B + (1-κ_E)*syn_Enz,
    D(α_R) ~ 0.0,
    # C cycle
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
    # N cycle
    β_NL ~ L/L_N,
    β_NR ~ R/R_N,
    D(L_N) ~ -dec_L/β_NL + i_L/β_Ni,
    D(R_N) ~ -dec_R/β_NR + ϵ_tvr*tvr_B/β_NB + (1-κ_E)*syn_Enz/β_NE,
    u_PlantN ~ min(u_PlantNmax, k_PlantN*I_N), # TODO specify as Pin 
    D(I_N) ~ i_IN - u_PlantN - l_N*I_N + Φ_N,
    Φ_N ~ Φ_Nu + Φ_NB + r_tvr/β_NB,
    Φ_Nu ~ (1-ν_N) * u_NOM,
    u_N ~ ν_N * u_NOM + u_immNPot,
    u_NOM ~ dec_L/β_NL + dec_R/β_NR + κ_E*syn_Enz/β_NE,
    u_immNPot ~ i_BN * I_N,
    N_synBN ~ u_N - syn_Enz/β_NE,
    M_ImbN ~ u_N - (syn_B/β_NB + syn_Enz/β_NE),
    Φ_NB ~ M_ImbN - u_immNPot,
    # co-limitation
    #syn_B ~ C_synBC,
    syn_B ~ min(C_synBC, β_NB*N_synBN), # TODO add P limitation
    α_L ~ (1-α_R),
    # excell fluxes
    r_O ~ u_C - (syn_Enz/ϵ +  syn_B + r_G + r_M), 
    # auxilaries
    #D(cumresp) ~ r_B + r_tvr,
    r_tvr ~ (1-ϵ_tvr)*tvr_B,
    r_B ~ r_GEnz + r_G + r_M + r_O,
    r_GEnz ~ (1-ϵ)/ϵ * syn_Enz,
    r_G ~ IfElse.ifelse(syn_B > 0, (1-ϵ)/ϵ * syn_B, 0.0), 
    i_L ~ i_L0, # TODO specify as Pin and compose with different system
    i_IN ~ I_IN0, # TODO specify as Pin and compose with different system
    β_Ni ~ β_Ni0,
    ]


@named sesam1 = ODESystem(eqs)    
sesam1 = structural_simplify(sesam1)
states(sesam1)
equations(sesam1)

u0 = Dict(
    B => 17,
    L => 100,
    R => 1100,
    cumresp => 0.0,
    α_R => 0.1, # enzyme synthesis into L # TODO model by optimality
    I_N => 0.04, ##<< inorganic pool gN/m2 
)
u0[L_N] = u0[L]/p[β_Ni0]
u0[R_N] = u0[R]/p[β_NB]

# u0[β_NL] = p[β_Ni0]
# u0[β_NR] = p[β_NB]











