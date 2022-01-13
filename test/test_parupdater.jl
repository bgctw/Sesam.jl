function simple_system(;name=:sys)
    @parameters t, p1, p2, p3
    D = Differential(t)
    @variables tmp(t)
    sys = ODESystem([D(tmp)~p1 + p2 + p3]; name)
    u0 = Dict(tmp => 1.0)
    p0 = Dict(p1 => 1.1, p2 => 1.2, p3 => 1.3)
    prob = ODEProblem(sys, u0, (0.0,1.0), p0)
    (sys, prob, u0, p0)
end

@testset "SystemParUpdater" begin
    sys, prob, u0, p0 = simple_system()
    parameters(sys)
    d_up = Dict(p3 => 0.3, p1 => 0.1)
    upd! = SystemParUpdater(d_up,sys)
    pvec = copy(prob.p)
    upd!(pvec, d_up)
    @test pvec[prob.p .== 1.1] == [0.1]
    @test pvec[prob.p .== 1.2] == [1.2]
    @test pvec[prob.p .== 1.3] == [0.3]
    #prob2 = remake(prob, p = pvec)
end;

