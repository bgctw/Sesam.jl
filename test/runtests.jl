using Test, SafeTestsets
const GROUP = get(ENV, "GROUP", "All") # defined in in CI.yml
@show GROUP

@time begin
    if GROUP == "All" || GROUP == "Basic"
        #@safetestset "Tests" include("test/test_sesam3_protect.jl")
        @time @safetestset "sesam3_protect" include("test_sesam3_protect.jl")
        #@safetestset "Tests" include("test/test_sesam3.jl")
        @time @safetestset "sesam3" include("test_sesam3.jl")
        #@safetestset "Tests" include("test/test_sesam3_regressionR.jl")
        @time @safetestset "sesam3 regression to R" include("test_sesam3_regressionR.jl")
        #@safetestset "Tests" include("test/test_sesam3_revMM.jl")
        @time @safetestset "sesam3_revMM" include("test_sesam3_revMM.jl")
        #@safetestset "Tests" include("test/test_sesam3CN.jl")
        @time @safetestset "sesam3_CN" include("test_sesam3CN.jl")
        #@safetestset "Tests" include("test/test_seam3.jl")
        @time @safetestset "seam3" include("test_seam3.jl")
        #@safetestset "Tests" include("test/test_plant_face.jl")
        @time @safetestset "plant_face" include("test_plant_face.jl")
        #@safetestset "Tests" include("test/test_plant_face_fluct.jl")
        @time @safetestset "plant_face_fluct" include("test_plant_face_fluct.jl")
    end
end

