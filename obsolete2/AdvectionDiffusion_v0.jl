using Revise, Printf
using LoopVectorization
import Plots
using LinearAlgebra, SparseArrays 
import UnicodePlots
using Base.Threads
##############
include("DataStructures.jl")
include("ThermalRoutines.jl")
include("MarkerRoutines.jl")
include("GridRoutines.jl")
##############
function SetMarkers!( p, R )
@tturbo for k=1:p.nmark
    in         = (p.y[k]^2) < R^2  #p.x[k]^2 + 
    p.phase[k] = (in==0)* 1.0 + (in==1)*2.0
end
end
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
ncx       = 40
ncy       = 40
nt        = 100
nout      = 50
dx, dy    = (xmax-xmin)/ncx, (ymax-ymin)/ncy
dt        = 1
nmx       = 4            # 2 marker per cell in x
nmy       = 4            # 2 marker per cell in y
nmark     = ncx*ncy*nmx*nmy; # total initial number of marker in grid
# RK4 weights
rkw = 1.0/6.0*[1.0 2.0 2.0 1.0] # for averaging
rkv = 1.0/2.0*[1.0 1.0 2.0 2.0] # for time stepping
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
k2        = 10
steady    = 1
# Allocate tables
T         = zeros( ncx  , ncy  )
T0        = zeros( ncx  , ncy  )
kc        = zeros( ncx  , ncy  )
rho       =  ones( ncx  , ncy  )
Cp        =  ones( ncx  , ncy  )
H         = zeros( ncx  , ncy  )
Vx        = zeros( ncx+1, ncy+2) # !!! GHOST ROWS
Vy        = zeros( ncx+2, ncy+1) # !!! GHOST COLUMNS
# Initialise markers
nmark0    = ncx*ncy*nmx*nmy; # total initial number of marker in grid
dx,  dy   = (xmax-xmin)/ncx, (ymax-ymin)/ncy
dxm, dym  = dx/nmx, dy/nmy 
xm1d      =  LinRange(xmin+dxm/2, xmax-dxm/2, ncx*nmx)
ym1d      =  LinRange(ymin+dym/2, ymax-dym/2, ncy*nmy)
(xmi,ymi) = ([x for x=xm1d,y=ym1d], [y for x=xm1d,y=ym1d])
xc        =  LinRange(xmin+dx/2, xmax-dx/2, ncx)
yc        =  LinRange(ymin+dy/2, ymax-dy/2, ncy)
# Over allocate markers
nmark_max = 3*nmark0;
phm    = zeros(Float64, nmark_max)
xm     = zeros(Float64, nmark_max) 
ym     = zeros(Float64, nmark_max)
Tm     = zeros(Float64, nmark_max)
cellxm = zeros(Int64,   nmark_max)
cellym = zeros(Int64,   nmark_max)
xm[1:nmark0]     = vec(xmi)
ym[1:nmark0]     = vec(ymi)
phm[1:nmark0]    = zeros(Float64, size(xmi))
cellxm[1:nmark0] = zeros(Int64,   size(xmi)) #zeros(CartesianIndex{2}, size(xm))
cellym[1:nmark0] = zeros(Int64,   size(xmi))
p      = Markers( xm, ym, Tm, phm, cellxm, cellym, nmark0, nmark_max )
mpc         = zeros(Float64,(ncx  ,ncy  )) # markers per cell
mpc_th      = zeros(Float64,(nthreads(), ncx  ,ncy   )) # markers per cell
# Initial configuration
SetMarkers!( p, sig )  # Define phase on markers
xc2   = repeat(xc, 1, length(yc))
yc2   = repeat(yc, 1, length(xc))'
xvx2  = repeat(xv, 1, length(yce))
yvx2  = repeat(yce, 1, length(xv))'
xvy2  = repeat(xce, 1, length(yv))
yvy2  = repeat(yv, 1, length(xce))'
Lx, Ly = xmax - xmin, ymax - ymin
ax   = 1.0/2.0; ay = 1.0/2.0
Vx    .=   cos.(ax*2.0 .*pi.*xvx2./Lx) .* sin.(ay*2.0.*pi .* yvx2./Ly)
Vy    .=  -sin.(ax*2.0 .*pi.*xvy2./Lx) .* cos.(ay*2.0.*pi .* yvy2./Ly)
Vxc    = 0.5.*(Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
Vyc    = 0.5.*(Vy[2:end-1,1:end-1] .+ Vy[2:end-1,2:end-0])
# BC definition
BC = BoundaryConditions()
BC.T.type_W = 2;     BC.T.Dir_W = 0.0*ones(ncy)
BC.T.type_E = 2;     BC.T.Dir_E = 0.0*ones(ncy)       
BC.T.type_S = 1;     BC.T.Dir_S = 1.0*ones(ncx)
BC.T.type_N = 1;     BC.T.Dir_N = 0.0*ones(ncx) 
# TIME LOOP
for it=1:nt
    @printf("Time step #%04d\n", it)
    # Count markers and reseed
    @time CountMarkers!(p,mpc,mpc_th,nmark,xmin,xmax,ymin,ymax,xc,yc,dx,dy,ncx,ncy)
    @time ReseedMarkers!( p, mpc, xc, yc, dx, dy, ncx, ncy)
    nmark = p.nmark
    @time CountMarkers!(p,mpc,mpc_th,nmark,xmin,xmax,ymin,ymax,xc,yc,dx,dy,ncx,ncy)
    # Define dt
    maxVx = maximum(abs.(Vx))
    maxVy = maximum(abs.(Vy))
    dt    = 0.25*(maximum([dx,dy])/maximum([maxVx,maxVy]))
    @printf("dt = %2.2e\n", dt)
    # Update cell info on markers
    LocateMarkers(p,dx,dy,xc,yc,xmin,xmax,ymin,ymax)
    # Interpolate k from markers
    kc = zeros( ncx  , ncy  )
    @time Markers2Cells3!(p,kc,xc,yc,dx,dy,ncx,ncy,[k1,k2],0,0)
    @time Markers2Cells3!(p,T0,xc,yc,dx,dy,ncx,ncy,    p.T,0,1)
    # Initialize
    kx, ky = AverageConductivity( kc, ncx, ncy, BC.T )
    # Residual (-> RHS)
    @time F =  ThermalResidual( T, T0, rho, Cp, H, kx, ky, dx, dy, dt, steady, BC.T, ncx, ncy )
    println("Initial residual = ", norm(F))
    # Matrix assembly
    @time  K = ThermalAssembly( Cp, rho, kx, ky, dx, dy, dt, steady, BC.T, ncx, ncy )
    # Solve
    @time ThermalSolve!( T, K, F, ncx, ncy )
    # Residual (Check)
    F =  ThermalResidual( T, T0, rho, Cp, H, kx, ky, dx, dy, dt, steady, BC.T, ncx, ncy )
    println("Solve   residual = ", norm(F))
    @time RungeKutta!(p, nmark, rkv, rkw, BC, dt, Vx, Vy, xv, yv, xce, yce, dx, dy, ncx, ncy)

    if mod(it, nout)==0
    # Visualize
    # p1 = Plots.heatmap(xv*Lc, yce*Lc, Array(Vx)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="T")
    # p1 = Plots.heatmap(xce*Lc, yv*Lc, Array(Vy)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="T")
    p1 = Plots.heatmap(xc*Lc, yc*Lc, Array((T*Tc.-0*273.15))', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="T")
    # stp = 1
    # Vxc2 = Vxc
    # Vyc2 = Vyc
    # println(size(Vx))
    # println(size(Vy))
    # p1 = Plots.quiver!(xc2[1:stp:end]*Lc, yc2[1:stp:end]*Lc, quiver=(Vxc2[1:stp:end], Vyc2[1:stp:end]))
    # p1 = Plots.heatmap(xc*Lc, yc*Lc, Array((T*Tc.-0*273.15))', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="T")
    # p1 = Plots.scatter!(p.x[p.phase.==1], p.y[p.phase.==1], c=:white,  markersize=0.5, alpha=0.8, legend=false)
    # p1 = Plots.scatter!(p.x[p.phase.==2], p.y[p.phase.==2], c=:green,  markersize=1.0, alpha=0.8, legend=false)
    # p1 = Plots.scatter!(p.x[p.phase.==1], p.y[p.phase.==1], c=:white, markershape=:circle,  markersize=2.0, alpha=0.5, legend=false, markerstrokewidth = 0)
    # p1 = Plots.scatter!(p.x[p.phase.==2], p.y[p.phase.==2], c=:white, markershape=:hexagon,  markersize=2.0, alpha=1.0, legend=false, markerstrokewidth = 0)
    display(Plots.plot( p1, dpi=200 ) ) 
    end
end
end

for it=1:1
    @time main()
end