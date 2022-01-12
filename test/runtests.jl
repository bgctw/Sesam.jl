using Sesam
using ModelingToolkit, DifferentialEquations
using Test

@testset "sesam3" begin
    include("test_sesam1.jl")        
end

@testset "plant_face" begin
    include("test_plant_face.jl")        
end
