# FDMIC_Geodynamics
Basic Finite difference / Marker-In-Cell code for teaching purposes.

# Stokes solver
The script [Mechanics_v11_SolVI_ViscousInclusionTest.jl](./Mechanics_v11_SolVI_ViscousInclusionTest.jl) produces first order L2 convergence of velocity and pressure using the SOLVI tests:

<!-- ![](/images/SOLVI_Julia.png) -->
<center><img src="/images/SOLVI_Julia.png" alt="drawing" width="500"/></center>

# Power-law Stokes flow + advection

The script [Mechanics_v11_MultiLayerExt.jl](./Mechanics_v11_MultiLayerExt.jl) allows to model multi-layer necking instabilities:

![](/images/MLPS_Julia.png)

<!-- ![](/images/MultiLayerExtension.gif) -->
<center><img src="/images/MultiLayerExtension.gif" alt="drawing" width="600"/></center>

# Periodic simple shear Stokes flow + advection

With [Mechanics_v11_PeriodicSimpleShear.jl](./Mechanics_v11_PeriodicSimpleShear.jl) one can model periodic shear deformation.

<!-- #![](/images/Periodic_Julia.png) -->

<!-- ![](/images/PeriodicSimpleShear.gif) -->
<center><img src="/images/PeriodicSimpleShear.gif" alt="drawing" width="600"/></center>

# Compressible layer under layer parallel compression: pressure build-up

<!-- ![](/images/CompressibleLayer.gif) -->
