
@testset "non-negative pools" begin
    #st = first(unknowns(sp))
    for st in unknowns(sp)
        @test all(sol[st] .>= 0.0)
    end
end;

@testset "microbial dC balance" begin
    #plot(sol,  vars=[s.u_C, s.r_B+s.syn_B+s.syn_Enz])
    @test all(isapprox.(sol[s.u_C],
        sol[s.r_B + s.syn_B + s.syn_Enz],
        rtol = 1e-6))
end;

@testset "system dC balance" begin
    change = sol[s.dB + s.dL + s.dR]
    output = sol[s.r_tot]
    input = sol[s.i_L]
    #plot(sol.t, [input, output, change, change .+ output])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;

@testset "microbial dN balance" begin
    uptake = p[s.ν_N] * sol[s.u_NOM] + sol[s.resorp_N]
    usage = sol[s.syn_Enz] / p[s.β_NEnz] + sol[s.syn_B] / p[s.β_NB] + sol[s.Φ_NB]
    #plot(sol.t, [uptake, usage])
    @test all(isapprox.(uptake, usage, rtol = 1e-6))
end;

@testset "system dN balance" begin
    change = sol[s.dB] / p[s.β_NB] + sol[s.dL_N + s.dR_N + s.dI_N]
    output = sol[s.u_PlantN + s.leach_N]
    input = sol[s.i_L / s.β_Ni + s.i_IN]
    #sol[s.leach_N]
    #sol[s.i_IN]
    bo = 1:10
    bo = 1:length(input)
    #plot(sol.t[bo], [input[bo], output[bo], change[bo], (change .+ output)[bo]], label = ["input" "output" "change" "change+output"])
    @test all(isapprox.(input, change .+ output, rtol = 1e-6))
end;

@testset "balanced allocation sums to one" begin
    #plot(sol, vars=[s.α_LT,s.α_RT, s.α_LT+s.α_RT])
    #@test all(isapprox.(sol[s.α_LT + s.α_RT], 1.0, rtol = 1e-8))
    @test all(isapprox.(sol[s.α_L + s.α_R], 1.0, rtol = 1e-8))
end;
