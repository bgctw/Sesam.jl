using Sesam
using Test
using JET: JET

@testset "JET" begin
    @static if VERSION ≥ v"1.9.2"
        JET.test_package(Sesam; target_modules = (@__MODULE__,))
    end
end;
# JET.report_package(Sesam) # to debug the errors
# JET.report_package(Sesam; target_modules=(@__MODULE__,)) # to debug the errors
