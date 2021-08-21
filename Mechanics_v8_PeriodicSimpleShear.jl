##############
using Revise
using FDMIC_Geodynamics
using LoopVectorization, Printf, Base.Threads, Plots, Revise, LinearAlgebra, Statistics, SparseArrays
##############
function SetMarkers!( p, R, xmax, xmin )
    @tturbo for k=1:p.nmark
        in         = ( ( (p.x[k])^2 + p.y[k]^2) < R^2 )
        p.phase[k] = (in==0)* 1.0 + (in==1)*2.0
    end
end
##############
@views function main( )
    println("#####################")
    println("###### M2Di.jl ######")
    println("#####################")
    # Parameters
    gamma     = 1e4
    solver    = 1   # 0: coupled --- 1: decoupled 
    Dir_scale = 1.0
    show_figs = 1
    # Domain
    xmin      = -1.0
    xmax      =  1.0
    ymin      = -1.0
    ymax      =  1.0
    Ebg       = -1.0
    # Scales
    Lc        = 1.0
    Tc        = 1.0
    tc        = 1.0/abs(Ebg)       
    # Numerics
    ncx       = 50
    ncy       = 50
    nt        = 1
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
    rad       = 1.0/3.0 
    eta1      = 1.0 # matrix viscosity
    eta2      = 1e3 # inclusion viscosity
    # BC definition
    BC = BoundaryConditions()
    ## Pure shear - full Dirichlet
    BC.periodix  = 1
    BC.Vx.type_W = 0;     BC.Vx.Dir_W = 0.0*ones(ncy+2)
    BC.Vx.type_E = 0;     BC.Vx.Dir_E =-0.0*ones(ncy+2)      
    BC.Vx.type_S =11;     BC.Vx.Dir_S =-1.0*ones(ncx+1)
    BC.Vx.type_N =11;     BC.Vx.Dir_N = 1.0*ones(ncx+1)
    BC.Vy.type_W = 0;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
    BC.Vy.type_E = 0;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
    BC.Vy.type_S = 1;     BC.Vy.Dir_S =-0.0*ones(ncx+2)
    BC.Vy.type_N = 1;     BC.Vy.Dir_N = 0.0*ones(ncx+2)
    # Allocate tables
    etac      = zeros(Float64, ncx  , ncy  )
    etav      = zeros(Float64, ncx+1, ncy+1)
    rhoc      = zeros(Float64, ncx  , ncy  )
    Pc        = zeros(Float64, ncx+0, ncy+0)
    Vx        = zeros(Float64, ncx+1, ncy+2) # !!! GHOST ROWS
    Vy        = zeros(Float64, ncx+2, ncy+1) # !!! GHOST COLUMNS
    div       = zeros(Float64, ncx+0, ncy+0)
    Tau       = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Eps       = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    if BC.periodix==0 
        Fx        = zeros(Float64, ncx+1, ncy+0)
    else
        Fx        = zeros(Float64, ncx+0, ncy+0)
    end
    Fy        = zeros(Float64, ncx+0, ncy+1)
    Fp        = zeros(Float64, ncx  , ncy  )
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
    cellxm[1:nmark0] = zeros(Int64,   size(xmi)) 
    cellym[1:nmark0] = zeros(Int64,   size(xmi))
    p                = Markers( xm, ym, Tm, phm, cellxm, cellym, nmark0, nmark_max )
    mpc              = zeros(Float64,(ncx  ,ncy  ))              # markers per cell
    mpc_th           = zeros(Float64,(nthreads(), ncx  ,ncy   )) # markers per cell per thread
    # Initial configuration
    SetMarkers!( p, rad, xmax, xmin )  # Define phase on markers
    xc2   = repeat(xc, 1, length(yc))
    yc2   = repeat(yc, 1, length(xc))'
    xvx2  = repeat(xv, 1, length(yce))
    yvx2  = repeat(yce, 1, length(xv))'
    xvy2  = repeat(xce, 1, length(yv))
    yvy2  = repeat(yv, 1, length(xce))'
    Lx, Ly = xmax - xmin, ymax - ymin
    ax   = 1.0/2.0; ay = 1.0/2.0
    # Initial guess - Vx
    Vx[1    ,:] .= BC.Vx.Dir_W
    Vx[end  ,:] .= BC.Vx.Dir_E
    Vx[:,    2] .= BC.Vx.Dir_S
    Vx[:,end-1] .= BC.Vx.Dir_N
    Vy[2    ,:] .= BC.Vy.Dir_W
    Vy[end-1,:] .= BC.Vy.Dir_E
    Vy[:,    1] .= BC.Vy.Dir_S
    Vy[:,  end] .= BC.Vy.Dir_N

    Vx_ini = zeros(ncx+1,ncy+2)
    Vx_ini .= Vx
    Vy_ini = zeros(ncx+2,ncy+1)
    Vy_ini .= Vy

    for i=1:size(Vx,1)
        for j=1:size(Vx,2)
            wW = 1.0 - (xv[i]-xmin)/(xmax-xmin)
            Vx[i,j] = wW * BC.Vx.Dir_W[j] + (1.0-wW) * BC.Vx.Dir_E[j]
            # if i>1 & i<size(Vx,1)
            #     wS = 1.0 - (yce[j]-ymin)/(ymax-ymin)
            #     Vx[i,j] += wS * BC.Vx.Dir_S[i] + (1.0-wS) * BC.Vx.Dir_N[i]
            # end
        end
    end
    
    # Initial guess - Vy
    for i=1:size(Vy,1)
        for j=1:size(Vy,2)
            wS = 1.0 - (yv[j]-ymin)/(ymax-ymin)
            Vy[i,j] = wS * BC.Vy.Dir_S[i] + (1.0-wS) * BC.Vy.Dir_N[i]
            # if j>1 & j<size(Vy,2)
            #     wW = 1.0 - (xce[i]-xmin)/(xmax-xmin)
            #     Vy[i,j] += wW * BC.Vy.Dir_W[j] + (1.0-wW) * BC.Vy.Dir_E[j]
            # end
        end
    end
    Vxc    = 0.5.*(Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
    Vyc    = 0.5.*(Vy[2:end-1,1:end-1] .+ Vy[2:end-1,2:end-0])
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
        @time LocateMarkers(p,dx,dy,xc,yc,xmin,xmax,ymin,ymax)
        # Interpolate k from markers
        @time Markers2Cells3!(p,etac,xc,yc,dx,dy,ncx,ncy,[eta1,eta2],0,0)
        # Initialize
        @time  CentroidsToVertices!( etav, etac, ncx, ncy, BC )

        # Evaluate residuals
        println("Residuals")
        @time SetBCs( Vx, Vy, BC )
        @time StrainRate!(div, Eps, Vx, Vy, dx, dy)
        @time Stress!(Tau, Eps, etac, etav)
        @time Residuals!(Fx, Fy, Fp, Tau, Pc, div, BC, ncx , ncy, dx, dy)
        println("|Fx| = ", norm(Fx)/length(Fx))
        println("|Fy| = ", norm(Fy)/length(Fy))
        println("|Fp| = ", norm(Fp)/length(Fp))
    
        # Number equations
        println("Numbering")
        @time NumVx, NumVy, NumP = NumberingStokes(BC, ncx, ncy)
    
        # Assemble Stokes matrices
        println("Assembly")
        @time Kuu, Kup, Kpu = StokesAssembly( BC, NumVx, NumVy, NumP, etac, etav, Dir_scale, dx, dy )
    
        # Call solver
        println("Solver")
        @time StokesSolver!(Vx,Vy,Pc,NumVx,NumVy,NumP, Fx, Fy, Fp,Kuu,Kup,Kpu,etac,gamma,solver)
    
        # Evaluate residuals
        println("Residuals")
        @time SetBCs( Vx, Vy, BC )
        @time StrainRate!(div, Eps, Vx, Vy, dx, dy)
        @time Stress!(Tau, Eps, etac, etav)
        @time Residuals!(Fx, Fy, Fp, Tau, Pc, div, BC, ncx , ncy, dx, dy)
        println("|Fx| = ", norm(Fx)/length(Fx))
        println("|Fy| = ", norm(Fy)/length(Fy))
        println("|Fp| = ", norm(Fp)/length(Fp))
    end
    # Visualize
    if show_figs==1
        p1 = Plots.heatmap(xv*Lc, yce*Lc, Array(Vx)',   aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="Vx")
        p2 = Plots.heatmap(xce*Lc, yv*Lc, Array(Vy)',   aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="Vy")
        p3 = Plots.heatmap(xc*Lc,  yc*Lc, Array(Pc)',   aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="P")
        p4 = Plots.heatmap(xv*Lc,  yv*Lc, Array(etav)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="etav")
        display(Plots.plot( p1, p2, p3, p4, dpi=200 ) ) 
    end
    return
end
####################################################
####################################################
####################################################
main()
