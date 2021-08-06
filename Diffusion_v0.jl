using Revise
using LoopVectorization
import Plots
using LinearAlgebra, SparseArrays 
import UnicodePlots
using Base.Threads
##############
include("DataStructures.jl")
include("ThermalRoutines.jl")
##############
@views function main()
# Domain
xmin      = -1.0
xmax      =  1.0
ymin      = -1.0
ymax      =  1.0
# Scales
Lc        = 1.0
Tc        = 1.0
# Numerics
ncx       = 10
ncy       = 10
dx, dy    = (xmax-xmin)/ncx, (ymax-ymin)/ncy
dt        = 1
# Initialize coordinates
xc        = LinRange( xmin+dx/2, xmax-dx/2, ncx  )
yc        = LinRange( ymin+dy/2, ymax-dy/2, ncy  )
xce       = LinRange( xmin-dx/2, xmax+dx/2, ncx+2)
yce       = LinRange( ymin-dy/2, ymax+dy/2, ncy+2)
xv        = LinRange( xmin     , xmax     , ncx+1)
yv        = LinRange( ymin     , ymax     , ncy+1)
# Parameters
Tamp      = 1.0
sig       = 0.3
k1        = 1.0
k2        = 10.0
steady    = 1
# Allocate tables
T         = zeros( ncx  , ncy  )
T0        = zeros( ncx  , ncy  )
kc        = zeros( ncx  , ncy  )
rho       =  ones( ncx  , ncy  )
Cp        =  ones( ncx  , ncy  )
H         = zeros( ncx  , ncy  )
Vx        = zeros( ncx+1, ncy  )
Vy        = zeros( ncx  , ncy+1)
# Initial configuration
xc2  = repeat(xc, 1, length(yc))
yc2  = repeat(yc, 1, length(xc))'
kc   = k1.*ones( ncx  , ncy  )
kc[ xc2.^2 .+ yc2.^2 .< sig.^2 ] .= k2
xvx2  = repeat(xv, 1, length(yc))
yvx2  = repeat(yc, 1, length(xv))'
xvy2  = repeat(xc, 1, length(yv))
yvy2  = repeat(yv, 1, length(xc))'
Lx, Ly = xmax - xmin, ymax - ymin
ax = 1/2.0; ay = 1/2
Vx    .=  -sin.(ax*2.0 .*pi.*xvx2./Lx) .* cos.(ay*2.0.*pi .* yvx2./Ly)
Vy    .=   cos.(ax*2.0 .*pi.*xvy2./Lx) .* sin.(ay*2.0.*pi .* yvy2./Ly)
Vxc    = 0.5.*(Vx[1:end-1,:] .+ Vx[2:end-0,:])
Vyc    = 0.5.*(Vy[:,1:end-1] .+ Vy[:,2:end-0])
# BC definition
BC = BoundaryConditions()
BC.T.type_W = 0;     BC.T.Dir_W = 0.0*ones(ncy)
BC.T.type_E = 0;     BC.T.Dir_E = 0.0*ones(ncy)       
BC.T.type_S = 1;     BC.T.Dir_S = 1.0*ones(ncx)
BC.T.type_N = 1;     BC.T.Dir_N = 0.0*ones(ncx) 
# Initialize
kx, ky = AverageConductivity( kc, ncx, ncy, BC.T )
# Residual (-> RHS)
F =  ThermalResidual( T, T0, rho, Cp, H, kx, ky, dx, dy, dt, steady, BC.T, ncx, ncy )
println("Initial residual = ", norm(F))
# Matrix assembly
K = ThermalAssembly( Cp, rho, kx, ky, dx, dy, dt, steady, BC.T, ncx, ncy )
# Solve
ThermalSolve!( T, K, F, ncx, ncy )
# Residual (Check)
F =  ThermalResidual( T, T0, rho, Cp, H, kx, ky, dx, dy, dt, steady, BC.T, ncx, ncy )
println("Solve   residual = ", norm(F))
# Visualize
p1 = Plots.heatmap(xc*Lc, yc*Lc, Array((T*Tc.-0*273.15))', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="T")
# p1 = Plots.heatmap(xc*Lc, yc*Lc, Array((kc*Tc.-0*273.15))', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="T")
display(Plots.plot( p1, dpi=200 ) ) 
end

for it=1:3
    @time main()
end