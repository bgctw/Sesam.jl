function plant_const(;name)
    @parameters t 
    D = Differential(t)
    @parameters i_L0 [unit = uQ/uT] i_IN0 [unit = uQ/uT] β_Ni0
    @parameters u_PlantNmax0 [unit = uQ/uT] k_PlantN0 [unit = uT^-1]
    @variables (begin
        i_L(t),[unit = uQ/uT], β_Ni(t), i_IN(t),[unit = uQ/uT],
        u_PlantNmax(t),[unit = uQ/uT], k_PlantN(t),[unit = uT^-1]
    end) 
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


