function plant_sesam_system(s,pl; name, simplify=true)
    @parameters t 
    sp = compose(ODESystem([
        s.i_L ~ pl.i_L,
        s.i_IN ~ pl.i_IN,
        s.β_Ni ~ pl.β_Ni,
        s.u_PlantNmax ~ pl.u_PlantNmax,
        s.k_PlantN ~ pl.k_PlantN,
      ], t; name), s, pl)
    simplify ? structural_simplify(sp) : sp
end
