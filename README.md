# FDMIC_Geodynamics
Basic Finite difference / Marker-In-Cell code for teaching purposes.

# Stokes solver
The script [Mechanics_v11_SolVI_ViscousInclusionTest.jl](./Mechanics_v11_SolVI_ViscousInclusionTest.jl) produces first order L2 convergence of velocity and pressure using the SOLVI tests. See [e.g. Duretz et al., 2011](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2011GC003567) for some details. The original analytical solution is given in [Schmid & Podladchikov, 2003](https://academic.oup.com/gji/article/155/1/269/713923).

<!-- ![](/images/SOLVI_Julia.png) -->
<center><img src="/images/SOLVI_Julia.png" alt="drawing" width="500"/></center>

# Power-law Stokes flow + advection

The script [Mechanics_v11_MultiLayerExt.jl](./Mechanics_v11_MultiLayerExt.jl) allows to model multi-layer necking instabilities. See [Duretz & Schmalholz, 2015](https://pubs.geoscienceworld.org/gsa/geology/article-abstract/43/8/711/131951/From-symmetric-necking-to-localized-asymmetric?redirectedFrom=fulltext) for more details.

![](/images/MLPS_Julia.png)

<!-- ![](/images/MultiLayerExtension.gif) -->
<center><img src="/images/MultiLayerExtension.gif" alt="drawing" width="600"/></center>

# Periodic simple shear Stokes flow + advection

With [Mechanics_v11_PeriodicSimpleShear.jl](./Mechanics_v11_PeriodicSimpleShear.jl) one can model periodic shear deformation.

<!-- #![](/images/Periodic_Julia.png) -->

<!-- ![](/images/PeriodicSimpleShear.gif) -->
<center><img src="/images/PeriodicSimpleShear.gif" alt="drawing" width="600"/></center>

# Compressible layer under layer parallel compression

&#8594; Pressure build-up

This can be reproduced with [Mechanics_v11_CompressibleLayer.jl](./Mechanics_v11_CompressibleLayer.jl).
See the paper of [Moulas et al., 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/jmg.12446) for more details.

<!-- ![](/images/CompressibleLayer.gif) -->
<center><img src="/images/CompressibleLayer.gif" alt="drawing" width="600"/></center>

