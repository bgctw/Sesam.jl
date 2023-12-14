using Sesam
using ModelingToolkit, OrdinaryDiffEq

# see test_sesam for constructing and solving the system

using Plots
plot(sol, vars = [s.L, s.B])
plot(sol, vars = [s.I_N])
plot(sol)
# syn_B is the minimum of C_synBC and C_synBN
plot(sol, vars = [s.C_synBC, s.C_synBN, s.syn_B])
plot(sol, vars = [s.β_NR])

plot(sol, vars = [s.α_LT, s.α_RT, s.α_LT + s.α_RT])
plot(sol, vars = [s.revenue_L, s.revenue_R])
plot(sol, vars = [s.lim_C, s.lim_N])

# --------- why is R accumulating despite N limitation?
sol[s.u_PlantNmax]

# B is C limited
plot(sol, vars = [s.Φ_N, s.Φ_Nu, s.Φ_NB, s.Φ_Ntvr])
plot(sol, vars = [s.Φ_NB])
plot(sol, vars = [s.u_PlantN, s.u_PlantNmax])
# entire flux by tvr and u is taken up by plant
# only after R stocks increase towards a steady-state, enough replenishment
plot(sol, vars = [s.dec_L / s.β_NL, s.dec_R / s.β_NR])
# 
plot(sol, vars = [s.dec_R, s.dec_RPot]) # negative dec_R?
plot(sol, vars = [s.α_R * s.syn_Enz]) # negative dec_R?
#plot(sol, vars=[s.α_R, s.α_RT, s.α_L, s.α_LT]) # negative dec_R?
plot(sol, vars = [s.α_R, s.α_L]) # other pools than R near steady-state, slowly accumulating, with it shifting towards its decomposition

Dict(p)[i_L]
# total C stock equals input
#plot(sol, vars=[R + L + B + cumresp])
#plot!(t -> Dict(p)[i_L]*t + sum(getindex.(Ref(Dict(u0)),[B,L,R])))
# biomass C balance
# plot(sol, tspan=(0.0,2.0), vars=[u_C, r_B+syn_B+syn_Enz])
# plot(sol,  vars=[s.u_C, s.r_B+s.syn_B+syn_Enz])
#
# biomass N balance
plot(sol, tspan = (0.0, 2.0), vars = [u_NOM, r_B + syn_B + syn_Enz])
# plot(sol,  vars=[u_C, r_B+syn_B+syn_Enz])
