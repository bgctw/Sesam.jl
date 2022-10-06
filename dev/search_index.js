var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Sesam","category":"page"},{"location":"#Sesam","page":"Home","title":"Sesam","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Sesam.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Sesam]","category":"page"},{"location":"#Sesam.calculate_β_NR_sesam3-Tuple{Any, Any}","page":"Home","title":"Sesam.calculate_β_NR_sesam3","text":"R pool is a mixture of microbial turnover and enzymes. Here the stable C/N ratio is computed based on given parameterization.\n\n\n\n\n\n","category":"method"},{"location":"#Sesam.calculate_β_PR_sesam3-Tuple{Any, Any}","page":"Home","title":"Sesam.calculate_β_PR_sesam3","text":"R pool is a mixture of microbial turnover and enzymes. Here the stable C/N ratio is computed based on given parameterization.\n\n\n\n\n\n","category":"method"},{"location":"#Sesam.plant_const-Tuple{}","page":"Home","title":"Sesam.plant_const","text":"System to be coupled with Sesam3 that provides constant inputs, i.e. drivers. Need to specify parameters\n\ni_L0: litter input rate\nβ_Ni0, β_Pi0: C:N and C:P ratio of litter input\ni_IN0, i_IP0: external inputs of inorganic N and P (N deposition, weathering)\n\noptional:\n\nu_PlantNmax0, u_PlantPmax0: maximum plant uptake of inorganic N and P\nk_PlantN0, k_PlantP0: plant uptake rate per inorganic pool of N and P \nu_Pup: plant inorganic P uplift (to potentially modify u_PlantPmax0)\n\nDefaults of plant uptake are set so that uptake rate is very high (1e3) so that the maximum is taken up, i.e. wins competition with microbes. The maximum plan uptake is set to match the input by litter.\n\n\n\n\n\n","category":"method"},{"location":"#Sesam.plant_face-Tuple{}","page":"Home","title":"Sesam.plant_face","text":"model iL and βNi as as step function increased by factor fac_int at t in [t1,t2)\n\n\n\n\n\n","category":"method"},{"location":"#Sesam.plant_face_fluct-Tuple{}","page":"Home","title":"Sesam.plant_face_fluct","text":"model iL and βNi as as step function increased by factor fac_int at t in [t1,t2). Superimpose a seasonal pattern on the litter input.\n\n\n\n\n\n","category":"method"}]
}
