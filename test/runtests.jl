using Sesam
using ModelingToolkit, DifferentialEquations
using Test
using DistributionFits
using MTKHelpers

@testset "sesam3" begin
    include("test_sesam3.jl")        
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

@testset "seam3" begin
    include("test_seam3.jl")        
end;


# @testset "parupdater" begin
#     include("test_parupdater.jl")        
# end;
