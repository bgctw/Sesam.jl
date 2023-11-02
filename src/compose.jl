function plant_sesam_system(s, pl; name, simplify = true)
    @parameters t
    eqCN = [
        s.i_L ~ pl.i_L,
        s.i_IN ~ pl.i_IN,
        s.β_Ni ~ pl.β_Ni,
        s.u_PlantNmax ~ pl.u_PlantNmax,
        s.k_PlantN ~ pl.k_PlantN,
    ]
    eqP = :i_IP ∈ propertynames(s) ?
          [
        s.i_IP ~ pl.i_IP,
        s.β_Pi ~ pl.β_Pi,
        s.u_PlantPmax ~ pl.u_PlantPmax,
        s.k_PlantP ~ pl.k_PlantP,
        s.s_EP ~ pl.s_EP,
    ] : Equation[]
    sp = compose(ODESystem(vcat(eqCN, eqP), t; name), s, pl)
    simplify ? structural_simplify(sp) : sp
end
