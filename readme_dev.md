The SESAM model is implemented using ModelingToolkit.

The CN-model resides in src/sesam3.jl 
The C fluxes are defined in function sesam3C,
the N fluxes are defined in function sesam3N, taking a C-system to extend,

They are combined in function sesam3CN, where
- syn_B is defined as the minimum of C_synBC and β_NB*N_synBN
- and quantities that are derived:
  - elemental limitation factors ω_Enz
  - α_L as (1- sum_opther_alpha)
  - dα_R via calling get_revenue_x
- equations for computing the revenue depending on 
  - seam3: (in seam3_revMM.jl) combines revenues of different elements by weighting factors (Wutzler17)
  - sesam3CN: combines investments and returns separately across elements
  - sesam3CN_deriv: compute dalpha by Wutzler22 Appendix D

The CNP model is based on the CN model but combines in sesam3CNP
- define syn_B, ω_Enz, ω_Z,  and to include P limitation
- TODO use unnormalized weights in omega

