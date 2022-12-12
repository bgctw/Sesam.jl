module Sesam

using ModelingToolkit, DifferentialEquations, IfElse
#using Unitful
using DistributionFits: fit_mode_flat
using Distributions
using MTKHelpers
using StaticArrays

export plant_const, plant_const_balanced, 
    plant_face, plant_face_fluct, plant_face_fluct_fake,
    sesam3, seam3, sesam3CN, plant_sesam_system,
    sesam3_protect, sesam3_revMM
export calculate_β_NR_sesam3, calculate_β_PR_sesam3
#export smoothstep
#export SystemParUpdater

#include("units.jl")
include("util.jl")
include("sesam3_revMM.jl")
include("sesam3.jl")
include("sesam3P.jl")
include("seam3.jl")
include("sesam3_protect.jl")
include("plant_const.jl")
include("plant_face.jl")
include("plant_face_fluct.jl")
include("compose.jl")
include("override.jl")

end
