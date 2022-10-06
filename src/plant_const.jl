"""
System to be coupled with Sesam3 that provides constant inputs, i.e. drivers.
Need to specify parameters
- `i_L0`: litter input rate
- `β_Ni0`, `β_Pi0`: C:N and C:P ratio of litter input
- `i_IN0`, `i_IP0`: external inputs of inorganic N and P (N deposition, weathering)
optional:
- `u_PlantNmax0`, `u_PlantPmax0`: maximum plant uptake of inorganic N and P
- `k_PlantN0`, `k_PlantP0`: plant uptake rate per inorganic pool of N and P 
- `u_Pup`: plant inorganic P uplift (to potentially modify `u_PlantPmax0`)

Defaults of plant uptake are set so that uptake rate is very high (1e3) so that
the maximum is taken up, i.e. wins competition with microbes.
The maximum plan uptake is set to match the input by litter.
"""
function plant_const(;name)
    # @parameters t [unit = uT]
    # D = Differential(t)
    # @parameters i_L0 [unit = uQ/uT] i_IN0 [unit = uQ/uT] β_Ni0
    # @parameters u_PlantNmax0 [unit = uQ/uT] k_PlantN0 [unit = uT^-1]
    # @variables (begin
    #     i_L(t),[unit = uQ/uT], β_Ni(t), i_IN(t),[unit = uQ/uT],
    #     u_PlantNmax(t),[unit = uQ/uT], k_PlantN(t),[unit = uT^-1]
    # end) 
    @parameters t 
    D = Differential(t)
    ps = @parameters (i_L0,  i_IN0,   
        u_PlantNmax0, u_Pup, 
        k_PlantN0, β_Ni0, i_IP0, β_Pi0, u_PlantPmax0, k_PlantP0, s_EP0)
    states = @variables (begin
        i_L(t), 
        β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t),
        β_Pi(t), i_IP(t),
        u_PlantPmax(t), k_PlantP(t),
        s_EP(t)
    end) 
    eqs = [
        i_L ~ i_L0, 
        β_Ni ~ β_Ni0,
        i_IN ~ i_IN0, 
        u_PlantNmax ~ u_PlantNmax0,
        k_PlantN ~ k_PlantN0,
        β_Pi ~ β_Pi0,
        i_IP ~ i_IP0, 
        u_PlantPmax ~ u_PlantPmax0,
        k_PlantP ~ k_PlantP0,
        s_EP ~ s_EP0, 
    ]
    defaults=Pair{Num, Any}[
        # rate not constrained, to make u_PlantNmax0 effective
        k_PlantN0 => 1.0e3, 
        # take as much N as provided with litter
        u_PlantNmax0 => i_L0/β_Ni0, 
        # rate not constrained
        k_PlantP0 => 1.0e3, 
        # take as much P as provided with litter
        u_PlantPmax0 => i_L0/β_Pi0, 
        # defaults in order to work with the CN-only version
        β_Pi0 => 50.0, i_IP0 => 0.0, s_EP0 => 0.0,
        # uplift of inorganic P from deeper soil
        u_Pup => 0.0,
    ]
    ODESystem(eqs, t, states, ps; name, defaults)    
end


