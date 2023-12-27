using Sesam
using Test
using Aqua

@testset "Sesam.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(
            Sesam;
            #unbound_args = false, # does not recognize NamedTuple{K, NTuple{N,E}}
            stale_deps = (ignore = [:Requires],),
            ambiguities = false, # many ambiguities in Symbolic and Stats packages
        )
    end;
    @testset "ambiguities package" begin
        Aqua.test_ambiguities(Sesam;)
    end;
end
