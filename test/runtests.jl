using Sesam
using ModelingToolkit, DifferentialEquations
using Test
using DistributionFits
using MTKHelpers
using Sesam: Sesam as CP
using ComponentArrays


@testset "sesam3_protect" begin
    #include("test/test_sesam3_protect.jl");        
    include("test_sesam3_protect.jl")        
end;

@testset "sesam3" begin
    #include("test/test_sesam3.jl");        
    include("test_sesam3.jl")        
end;

@testset "sesam3 regression to R" begin
    #include("test/test_sesam3_regressionR.jl")        
    include("test_sesam3_regressionR.jl")        
end;


@testset "sesam3_revMM" begin
    #include("test/test_sesam3_revMM.jl")        
    include("test_sesam3_revMM.jl")        
end;

@testset "sesam3_CN" begin
    #include("test/test_sesam3CN.jl");        
    include("test_sesam3CN.jl")        
end;

@testset "seam3" begin
    #include("test/test_seam3.jl");        
    include("test_seam3.jl")        
end;

@testset "plant_face" begin
    include("test_plant_face.jl")        
end;

@testset "plant_face_fluct" begin
    include("test_plant_face_fluct.jl")        
end;

@testset "util" begin
    include("test_util.jl")        
end;


# @testset "parupdater" begin
#     include("test_parupdater.jl")        
# end;
