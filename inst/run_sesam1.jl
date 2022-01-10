tspan = (0.0,100.0)    

#prob = ODEProblem(sesam1,[t for t in u0], tspan, [t for t in p])
prob = ODEProblem(sesam1, u0, tspan, p)
#prob = ODEProblem(sesam1,u0, tspan, p, jac=true)
sol = solve(prob)

using Plots
plot(sol, vars=[L,B])
plot(sol, vars=[I_N])
plot(sol)
plot(sol, vars=[β_NL])
# increase of C in system is the same as input
plot(sol, vars=[C_synBC, p[β_NB]*N_synBN, syn_B])
plot(sol, vars=[β_NL])

Dict(p)[i_L]
# total C stock equals input
#plot(sol, vars=[R + L + B + cumresp])
#plot!(t -> Dict(p)[i_L]*t + sum(getindex.(Ref(Dict(u0)),[B,L,R])))
# biomass C balance
# plot(sol, tspan=(0.0,2.0), vars=[u_C, r_B+syn_B+syn_Enz])
# plot(sol,  vars=[u_C, r_B+syn_B+syn_Enz])
#
# biomass N balance
 plot(sol, tspan=(0.0,2.0), vars=[u_NOM, r_B+syn_B+syn_Enz])
# plot(sol,  vars=[u_C, r_B+syn_B+syn_Enz])

