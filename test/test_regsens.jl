using Sesam

cols = (:par, :dType, :mode, :upper)
paramsModeUpperRows = [
    (:a, LogNormal, 0.001 * 365, 0.005 * 365),
    (:b, LogitNormal, 0.3, 0.9),
]
df_dist = rename!(DataFrame(columntable(paramsModeUpperRows)), collect(cols))

@testset "fitDistr" begin
    #st = first(unknowns(sp))
    for st in unknowns(sp)
        @test all(sol[st] .>= 0.0)
    end
end;
