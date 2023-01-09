using Sesam
import Sesam as CP
using ModelingToolkit, DifferentialEquations

# resembling parameters and states in Sesam.R test_modSesam3P.R

@named s = sesam3()
@named pl = plant_const_balanced()
@named plc = plant_const(;name=:pl)

@named sp = plant_sesam_system(s,pl)
states(sp)

sr = sesam3(;use_proportional_revenue=true, name=:s)
@named spr = plant_sesam_system(sr,pl)

sFS = CP.sesam_fixed_substrates(s)
@named spFS = plant_sesam_system(sFS,plc)
sFS_r = CP.sesam_fixed_substrates(sr) 
@named spFS_r = plant_sesam_system(sFS_r,plc)
states(spFS_r)

@parameters β_NR0 β_PR0
parms = Dict(
    s.β_NB => 7.16
    #, cnBW => 10    ##<< C/N ratio of cell walls (that go to R, )
    ,s.β_NEnz => 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
    #,cnE => 7.16
    ,β_NR0 => 4.5     ##<< between micr and enzyme signal
    ,pl.β_Ni0 => 30      ##<< N poor substrate
    #,kN => 0.05     ##<< (per day) enzyme turnover
    #,kN => 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
    ,s.k_mN_L => 0.3*0.01*365  ##<< /yr enzyme turnover 1% turning over each day
    ,s.κ_E => 0.8      ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R)
    #,kR => 0.2      ##<< substrate decomposition rate N-rich (here simulate large N stock)
    #,kL => 1        ##<< substrate decomposition rate N-poor
    #,kR => 5e-3     ##<< substrate decomposition rate N-rich (here simulate large N stock)
    #,kL => 10e-3    ##<< substrate decomposition rate N-poor
    #,kL => 5e-2     ##<< substrate decomposition rate N-poor
    #,aE => 0.05     ##<< C-uptake allocated to enzymes
    #,kR => 1/(5*365)      ##<< 5 years
    #,kL => 1/(0.5*365)    ##<< 1/2 years
    ,s.k_R => 1/(50)    ##<< 1/(x years)
    ,s.k_L => 1/(1)     ##<< 1/(x years)
    #,aE => 0.003*365 ##<< C biomass allocated to enzymes gC/day /microbial biomass
    ,s.a_E => 0.001*365 ##<< C biomass allocated to enzymes gC/day /microbial biomass
    #,km => 0.3       ##<< enzyme half-saturation constant
    #,km => 0.03     ##<< enzyme half-saturation constant
    #,km => 14       ##<< enzyme half-saturation constant
    ,s.m => 0.02*365   ##<< maintenance respiration rate   gC/day /microbial biomass
    ,s.τ => 1/60*365 ##<< biomass turnover rate (12 days)
    ,s.ϵ => 0.5      ##<< carbon use efficiency
    ,s.ϵ_tvr => 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
    #,iR => 0         ##<< input modelled explicitely
    ,pl.i_L0 => 300       ##<< g/m2 input per year (half NPP)
    #,plantNUp => 300/70*1/4  ##<< plant N uptake balancing N inputs
    ,plc.u_PlantNmax0 => 0   ##<< plant N uptake balancing N inputs
    #,useFixedAlloc => FALSE    ##<< set to true to use fixed enzyme allocation (alpha => 0.5)
    ,pl.k_PlantN0 => 10.57 #0.0289652*365         ##<< plant uptake iP IN
    ,s.i_BN => 0.4 #0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iBN IN
    ,pl.i_IN0 => 0         ##<< input of mineral N
    ,s.l_N => 0.96 #0.00262647*365       ##<< leaching rate of mineralN lN IN
    , s.ν_N => 0.9       ##<< microbial N use efficiency
    #, isEnzymeMassFlux => FALSE  ##<< steady state B solution neglects enyzme mass fluxes
    #, isEnzymeMassFlux => TRUE  ##<< steady state B solution accounts for enyzme mass fluxes
    , s.ν_P => 0.3      ##<< microbial uptake of depolymerized P, (1-nuP) is mineralized
    , s.β_PEnz => 50
    , s.β_PB => 40
    #, cpBW => 50
    , β_PR0 => 40
    #, pl.β_Pi0 => 40*3
    , pl.β_Pi0 => 40*6  # to make real P limitation possible
    , s.i_BP => 0.4 #0.38 * 10.57 # start with same as N
    , pl.s_EP0 => 0.3*0.01*365 /20  ##<< 1/10 of kmn: /yr enzyme turnover 1% turning over each day
    ,s.β_Pm => 500 # at a c:P ratio of cpm biomineralization rate decreases to 1/2
)
parms[s.k_mN_R] = parms[s.k_mN_L]  
merge!(parms, Dict(
    #  kmR => kmL => km
     # eps1 => eps2 => eps
      #cnER => cnEL => cnE
      #kNR => kNL => kN
      s.l_P => parms[s.l_N],       # leaching rate of inorganic P equals that of N
      s.ν_P => parms[s.ν_N],     # mineralization of P during decomposiition equals that of N
      pl.k_PlantP0 => parms[pl.k_PlantN0],  # plant uptake rate of P equals that of N
      pl.i_IP0 => 10.0,
      pl.i_IN0 => 10.0,
      plc.u_PlantNmax0 => parms[pl.i_L0] / parms[pl.β_Ni0],	# same litter input as plant uptake
      pl.k_PlantN0 => 0.0, 
      plc.u_PlantNmax0 => 0.0,			# no plant uptake
    ))
parms = CP.get_updated_Penz_pars(parms, s)    

x0C = ComponentVector(R=7000, L=200)
x0 = x0Orig = Dict( 
    s.B => 20                     ##<< microbial biomass
    , s.R => x0C.R                 ##<< N rich substrate
    , s.R_N => x0C.R/parms[β_NR0]    ##<< N rich substrate N pool
    , s.R_P => x0C.R/parms[β_PR0]
    , s.L => x0C.L                  ##<< N poor substrate
    , s.L_N => x0C.L/parms[pl.β_Ni0]     ##<< N poor substrate N pool
    , s.L_P => x0C.L/parms[pl.β_Pi0]     ##<< N poor substrate P pool
    , s.I_N =>  1                   ##<< inorganic N pool
    , s.I_P =>  1                  ##<< inorganic P pool
    , s.α_L => 0.4              ##<< initial community composition
    , s.α_R => 0.5              ##<< initial community composition
    , s.α_P => 0.1              ##<< initial community composition
)
#
x0C_Nlim = ComponentVector(R=1000, L=200)
x0Nlim = Dict( 
    s.B => 20                     ##<< microbial biomass
  , s.R => x0C_Nlim.R                 ##<< N rich substrate
  , s.R_N => x0C_Nlim.R/parms[β_NR0]    ##<< N rich substrate N pool
  , s.R_P => x0C_Nlim.R/parms[β_PR0]
  , s.L => x0C_Nlim.L                  ##<< N poor substrate
  , s.L_N => x0C_Nlim.L/parms[pl.β_Ni0]     ##<< N poor substrate N pool
  , s.L_P => x0C_Nlim.L/parms[pl.β_Pi0]      ##<< N poor substrate P pool
  , s.I_N =>  1e-12                   ##<< inorganic N pool
  , s.I_P =>  10                  ##<< inorganic P pool
  , s.α_L => 0.4              ##<< initial community composition
  , s.α_R => 0.5              ##<< initial community composition
  , s.α_P => 0.1              ##<< initial community composition
)


@testset "regression to R version, fixed substrates" begin
    tspan = (0,2100)
    probFS = ODEProblem(spFS, x0, tspan, parms)
    probFS_r = ODEProblem(spFS_r, x0, tspan, parms)
    sol = solLRP_r = solve(probFS_r, Rodas4());
    sol = solLRP = solve(probFS, Rodas4());
    # see Sesam R project test_modSesam3P.R
    @test isapprox(sol[s.B][end],11.1535271, atol=1e-5)
    @test isapprox(sol[s.α_R][end],0.4316128, atol=1e-5)
end;


@testset "regression to R version, fixed substrates Nlim" begin
    tspan = (0,8100)
    probFS_r = ODEProblem(spFS_r, x0Nlim, tspan, parms)
    sol = solLRP_r = solve(probFS_r, Rodas4());
    probFS = ODEProblem(spFS, x0Nlim, tspan, parms)
    sol = solLRP = solve(probFS, Rodas4());
    # see Sesam R project test_modSesam3P.R
    @test isapprox(sol[s.B][end],5.4536778 , atol=1e-5)
    @test isapprox(sol[s.α_R][end],0.3791681, atol=1e-5)
end;

@testset "regression to R version, fixed substrates CPlim" begin
    x0Plim = copy(x0Nlim)
    x0Plim[s.I_N] = 10.0
    x0Plim[s.I_P] = 0.002
    tspan = (0,8100)
    probFS_r = ODEProblem(spFS_r, x0Plim, tspan, parms)
    sol = solLRP_r = solve(probFS_r, Rodas4());
    probFS = ODEProblem(spFS, x0Plim, tspan, parms)
    sol = solLRP = solve(probFS, Rodas4());
    # see Sesam R project test_modSesam3P.R
    @test isapprox(sol[s.B][end],0.01057638  , atol=1e-5)
    @test isapprox(sol[s.α_P][end],0.68570386, atol=1e-5)
end;

@testset "regression to R version, variable substrates" begin
    tspan = (0,8000)
    prob = ODEProblem(sp, x0, tspan, parms)
    prob_r = ODEProblem(spr, x0, tspan, parms)
    sol = solLRP_r = solve(prob_r, Rodas4());
    sol = solLRP = solve(prob, Rodas4());
    # see Sesam R project test_modSesam3P.R
    @test isapprox(sol[s.B][end],16.6604961, atol=1e-5)
    @test isapprox(sol[s.α_R][end],0.2085514, atol=1e-5)
end;

@testset "regression to R version, variable substrates Nlim" begin
    tspan = (0,8000)
    parmsNlim = merge(parms, Dict(
      pl.i_IN0 => 1e-5,
      pl.β_Ni0 => 40.0, 
    ))
    prob = ODEProblem(sp, x0Nlim, tspan, parmsNlim)
    prob_r = ODEProblem(spr, x0Nlim, tspan, parmsNlim)
    sol = solLRP_r = solve(prob_r, Rodas4());
    sol = solLRP = solve(prob, Rodas4());
    # see Sesam R project test_modSesam3P.R
    @test isapprox(sol[s.B][end],15.6249494  , atol=1e-5)
    @test isapprox(sol[s.α_R][end],0.4151781, atol=1e-5)
end;

@testset "regression to R version, variable substrates CPlim" begin
    tspan = (0,8000)
    x0Plim = copy(x0)
    x0Plim[s.I_P] = 0.002
    parmsPlim = merge(parms, Dict(
      pl.i_IP0 => 1e-5,
    ))
    prob = ODEProblem(sp, x0Plim, tspan, parmsPlim)
    prob_r = ODEProblem(spr, x0Plim, tspan, parmsPlim)
    sol = solLRP_r = solve(prob_r, Rodas4());
    sol = solLRP = solve(prob, Rodas4());
    # see Sesam R project test_modSesam3P.R
    # TODO not equal
    @test isapprox(sol[s.B][end],5.4722157    , atol=1e-5)
    @test isapprox(sol[s.α_P][end],0.4034128, atol=1e-5)
end;


tmpf = () -> begin
    # using Plots
    #plot(sol, idxs=[s.L]) # constant
    [sol[s.dα_R][1], sol[s.dα_P][1]]
    [sol[s.C_synBC][1],  sol[s.C_synBN][1], sol[s.C_synBP][1]]
    sol[s.r_M][1]
    sol[s.syn_Enz][1]
    sol[s.resorp_N][1]
    [sol[s.syn_B][1], sol[s.tvr_B][1]]
    sol[s.dB][1]
    sol[s.ν_TP][1]
    [sol[s.lim_C][1], sol[s.lim_N][1], sol[s.lim_P][1]]
    sol[s.dec_RPPot][1]
    [sol[s.ν_TN][1], sol[s.ν_TP][1]]
    sol[s.τsyn][1]
    sol[s.mdu][1]
    sol[s.d_R][1]
    [sol[s.dec_LPot][1],sol[s.dec_RPot][1],sol[s.dec_LPPot][1],sol[s.dec_RPPot][1]]
    [sol[s.du_L][1],sol[s.du_R][1],sol[s.du_P][1]]
    [sol[s.ω_L][1],sol[s.ω_R][1],sol[s.ω_P][1]]
    [sol[s.dI_N][1], sol[s.dI_P][1]]
    [sol[s.p_uNmic][1], sol[s.p_uPmic][1]]
    [sol[s.u_immNPot][1], sol[s.u_immPPot][1]]
    [sol[s.u_PlantN][1], sol[s.u_PlantP][1]]
    [sol[s.u_PlantNmax][1], sol[s.u_PlantPmax][1]]
    getindex.(getindex.(Ref(sol), [s.i_IP, s.u_PlantP, s.leach_P, s.dec_LP_P, s.dec_RP_P, s.Φ_P]),1)
    getindex.(getindex.(Ref(sol), [s.dec_LPPot, s.dec_RPPot, s.lim_enz_P, s.s_EP]),1)

    sol[s.mdu][1]
    sol[s.mdu][1]
    sol[s.mdu][1]
    parms[s.k_mN_P]
    parms[s.β_PB]
    parms[s.ν_N]
    ts = (0,0.02)
    ts = (0,2)
    ts = (0,0.07)
    ts = (0,min(20.0, sol.t[end]))
    ts = (0,min(200.0, sol.t[end]))
    ts = (0,min(2000.0, sol.t[end]))
    plot(sol, idxs=[s.B], tspan=ts)
    plot(sol, idxs=[s.α_L], tspan=ts)
    plot(sol, idxs=[s.α_R], tspan=ts)
    plot(sol, idxs=[s.dα_R], tspan=ts)
    plot(sol, idxs=[s.α_L, s.α_R, s.α_P], tspan=ts)
    plot(sol, idxs=[s.lim_C, s.lim_N, s.lim_P], tspan=ts)
    plot(sol, idxs=[s.lim_LP, s.lim_RP], tspan=ts)
    plot(sol, idxs=[s.L], tspan=ts)
    plot(sol, idxs=[s.R], tspan=ts)
    plot(sol, idxs=[s.R_P], tspan=ts)
    plot(sol, idxs=[s.I_P], tspan=ts)
    plot(sol, idxs=[s.β_PR], tspan=ts)
    plot(sol, idxs=[s.tvr_B, s.syn_Enz, s.tvr_Enz], tspan=ts)
    #plot(sol, idxs=[s.β_PBtvr], tspan=ts)
    #plot(sol, idxs=[(parms[s.ϵ_tvr]*s.tvr_B + (1-parms[s.κ_E])*s.tvr_Enz)/(parms[s.ϵ_tvr]*s.tvr_B/s.β_PBtvr + (1-parms[s.κ_E])*s.tvr_Enz/parms[s.β_PEnz])], tspan=ts) # input to R is of C:P 40, but dec_RP decreases C:P of R pool

    # test julia solution in R
    popt = ComponentVector(s₊B = 1)
    pset = ProblemParSetter(sp, popt)
    #pset_r = ProblemParSetter(spLRP_r, popt)
    label_state(pset,sol[end])

end
