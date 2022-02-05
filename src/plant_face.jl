"""
model i_L and β_Ni as as step function increased by factor fac_int at t in [t1,t2)
"""
function plant_face(;name, t1=20, t2=20+50, fac_inc=1.2)
    @parameters t 
    @parameters i_L0  i_IN0  β_Ni0 t1=t1 t2=t2 fac_inc=fac_inc
    @parameters u_PlantNmax0  k_PlantN0 
    @variables (begin
        i_L(t), β_Ni(t), i_IN(t),
        u_PlantNmax(t), k_PlantN(t),
        i_L_annual(t)
    end) 
    eqs = [
        i_L ~ i_L_annual,
        i_L_annual ~ IfElse.ifelse((t1 <= t) & (t < t2), fac_inc*i_L0, i_L0), 
        β_Ni ~ IfElse.ifelse((t1 <= t) & (t < t2), fac_inc*β_Ni0, β_Ni0),
        i_IN ~ i_IN0, 
        u_PlantNmax ~ u_PlantNmax0,
        k_PlantN ~ k_PlantN0,
    ]
    defaults=Pair{Num, Any}[
        # rate not constrained
        k_PlantN0 => 1e3, 
        # take as much N as provide with litter
        u_PlantNmax0 => i_L0/β_Ni0, 
    ]
    continuous_events = [
        t ~ t1,
        t ~ t2,
    ]
    ODESystem(eqs,t; name, defaults, continuous_events)    
end


