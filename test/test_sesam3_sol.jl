# included from test_sesam3, sp, sol, p, u0 need to be defined

@testset "non-negative pools" begin
    #st = first(states(sp))
    for st in states(sp)
        #@show st
        #@test all(sol[st] .>= 0.0) 
        @test all(sol[st] .>= -eps(Float64)*100) 
    end
    # minimum(sol[s.α_LP])
    # minimum(sol[s.α_RP])
end;

@testset "microbial dC balance" begin
    #Plots.plot(sol,  vars=[s.u_C, s.r_B+s.syn_B+s.syn_Enz])
    @test all(isapprox.(
        sol[s.u_C], 
        sol[s.r_B+s.syn_B+s.syn_Enz], 
        rtol = 1e-6))
end;

@testset "system dC balance" begin
    change = sol[s.dB + s.dL + s.dR]
    output = sol[s.r_tot] 
    input = sol[s.i_L]
    #Plots.plot(sol.t, [input, output, change, change .+ output])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;


@testset "microbial dN balance" begin
    uptake = p[s.ν_N]*sol[s.u_NOM] + sol[s.resorp_N]
    usage = sol[s.syn_Enz]/p[s.β_NEnz] + sol[s.syn_B]/p[s.β_NB] + sol[s.Φ_NB]
    #Plots.plot(sol.t, [uptake, usage])
    @test all(isapprox.(uptake, usage, rtol = 1e-6))
end;

@testset "system dN balance" begin
    change = sol[s.dB]/p[s.β_NB] + sol[s.dL_N + s.dR_N + s.dI_N]
    output = sol[s.u_PlantN + s.leach_N] 
    input = sol[s.i_L/s.β_Ni + s.i_IN]
    #sol[s.leach_N]
    #sol[s.i_IN]
    bo = 1:10
    bo = 1:length(input)
    #Plots.plot(sol.t[bo], [input[bo], output[bo], change[bo], (change .+ output)[bo]], label = ["input" "output" "change" "change+output"])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;

@testset "balanced allocation sums to one" begin
    #Plots.plot(sol, vars=[s.α_LT,s.α_RT, s.α_LT+s.α_RT])
    #@test all(isapprox.(sol[s.α_LT + s.α_RT + s.α_PT], 1.0, rtol = 1e-8))
    @test all(isapprox.(sol[s.α_L + s.α_R + s.α_P], 1.0, rtol = 1e-8))
end;

@testset "microbial dP balance" begin
    uptake = p[s.ν_P]*sol[s.u_POM] + sol[s.resorp_P];
    usage = sol[s.syn_Enz]/p[s.β_PEnz] + sol[s.syn_B]/p[s.β_PB] + sol[s.Φ_PB];
    #Plots.plot(sol.t, [uptake, usage])
    @test all(isapprox.(uptake, usage, rtol = 1e-6))
    @test all(sol[s.resorp_P] .>= 0.0)
end;

@testset "system dP balance" begin
    change = sol[s.dB]/p[s.β_PB] + sol[s.dL_P + s.dR_P + s.dI_P]
    output = sol[s.u_PlantP + s.leach_P] 
    input = sol[s.i_L/s.β_Pi + s.i_IP]
    #sol[s.leach_N]
    #sol[s.i_IN]
    bo = 1:10
    bo = 1:length(input)
    bo = length(input) .- (10:-1:0)
    #Plots.plot(sol.t[bo], [input[bo], output[bo], change[bo], (change .+ output)[bo]], label = ["input" "output" "change" "change+output"])
    #Plots.plot(sol.t[bo], [input[bo], (change .+ output)[bo]], label = ["input" "change+output"])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;

