@parameters t
D = Differential(t)

@named plf = plant_face()
equations(plf)

function rep_plf(pl;name, simplify=true)
   sts = @variables x(t)
   tmp = compose(ODESystem([
       D(x) ~ pl.i_L
       ],t,sts,[];name), pl)
   simplify ? structural_simplify(tmp) : tmp
end
@named plf_rep = rep_plf(plf)

p_plf = Dict(
    plf.i_L0 => 400.0,         # g/m2 input per year (half NPP)
    plf.β_Ni0 => 25,
    plf.i_IN0 => 0,   ##<< input of mineral N,
)
tspan = (0, 120)
prob = ODEProblem(plf_rep, [0.0], tspan, p_plf)
sol = solve(prob, saveat = 0:1:tspan[2])

@testset "increased input" begin
    is_face_period = 20 .<= sol.t .< 20+50
    @test all(sol[plf.i_L][.!is_face_period] .== 400.0)
    @test all(sol[plf.β_Ni][.!is_face_period] .== 25.0)
    @test all(sol[plf.i_L][is_face_period] .== 1.2*400.0)
    @test all(sol[plf.β_Ni][is_face_period] .== 1.2*25.0)
end;

# using Plots
# plot(sol)
# plot(sol, vars=[plf.i_L])
# plot(sol, vars=[plf.β_Ni])