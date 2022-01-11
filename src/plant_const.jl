function plant_const(;name)
    @parameters t 
    D = Differential(t)
    @parameters i_L0 i_IN0 β_Ni0 u_PlantNmax0 k_PlantN0
    @variables i_L(t) i_IN(t) β_Ni(t) u_PlantNmax(t) k_PlantN(t)
    eqs = [
        i_L ~ i_L0, 
        β_Ni ~ β_Ni0,
        i_IN ~ i_IN0, 
        u_PlantNmax ~ u_PlantNmax0,
        k_PlantN ~ k_PlantN0,
    ]
    defaults=Pair{Num, Any}[
        # rate not constrained
        k_PlantN0 => 1000, 
        # take as much N as provide with litter
        u_PlantNmax0 => i_L0/β_Ni0, 
    ]
    ODESystem(eqs; name, defaults)    
end


