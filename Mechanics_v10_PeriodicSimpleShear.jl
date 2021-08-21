##############
using Revise
using FDMIC_Geodynamics
using LoopVectorization, Printf, Base.Threads, Plots, Revise, LinearAlgebra, Statistics, SparseArrays
##############
function SetMarkers!( p::Markers, R::Float64, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, dx::Float64, dy::Float64 )
    # Use this function to set up the model geometry
     for k=1:p.nmark #@tturbo <---------rand does not work in @tturbo!!!!
        in         = ( ( (p.x[k])^2 + p.y[k]^2) < R^2 )
        p.phase[k] = (in==0)* 1.0 + (in==1)*2.0
    end
end
##############
function Rheology!( etac::Matrix{Float64}, etav::Matrix{Float64}, Eps::Tensor2D, eta0::Vector{Float64}, n::Vector{Float64}, phase_perc::Array{Float64, 3}, nphase::Int64, BC::BoundaryConditions, ncx::Int64, ncy::Int64 )
    # Compute effective viscosity for each phase
    etac .= 0.0
    for m=1:nphase
        eta    = eta0[m] * Eps.II.^(1.0/n[m] - 1.0)
        etac .+= phase_perc[m,:,:] .* eta
    end
    @time  CentroidsToVertices!( etav, etac, ncx, ncy, BC )
    println("min Eii: ", minimum(Eps.II), " --- max Eii: ", maximum(Eps.II))
    println("min eta: ", minimum(etac), " --- max eta: ", maximum(etac))
end
##############
@views function main( )
    # Main routine
    println("#####################")
    println("###### M2Di.jl ######")
    println("#####################")
    # Domain
    xmin          = -1.0
    xmax          =  1.0
    ymin          = -1.0
    ymax          =  1.0
    Ebg           = -1.0      # reference strain rate
    # Scales
    Lc            = 1.0
    Tc            = 1.0
    tc            = 1.0/abs(Ebg)       
    # Spatial discretisation
    ncx           = 100
    ncy           = 100
    nmx           = 4            # 2 marker per cell in x
    nmy           = 4            # 2 marker per cell in y
    nmark         = ncx*ncy*nmx*nmy; # total initial number of marker in grid
    # Time discretisation
    nt            = 500
    dt            = 1
    Courant       = 0.4  # Courant number
    # Boundary conditions
    BC_type       = 3     # 1: Pure shear / 2: Simple shear / 3: Simple shear periodic
    PureShear_ALE = 0     # Deform box 
    # Solver
    solver        = 1     # 0: coupled --- 1: decoupled 
    gamma         = 1e4   # penalty factor
    Dir_scale     = 1.0   # Dirichlet scaling factor
    # Non-linear iterations 
    niter_nl      = 5     # max. number of non-linear iterations
    tol_nl        = 1e-3  # non-linear tolerance
    # Visualisation
    show_figs     = 1     # activates visualisation...
    nout          = 10    # ... every nout
    experiment    = "PeriodicSimpleShear"
    # RK4 weights
    rkw = 1.0/6.0*[1.0 2.0 2.0 1.0] # for averaging
    rkv = 1.0/2.0*[1.0 1.0 2.0 2.0] # for time stepping
    # Initialize coordinates
    L0            = xmax-xmin
    dx, dy, xc, yc, xce, yce, xv, yv = GenerateMesh( xmin, xmax, ymin, ymax, ncx, ncy )
    # Material parameters
    rad       = 1.0/3.0 
    nphase    = 2
    eta0      = zeros(Float64, nphase)
    n         = zeros(Float64, nphase)
    # Material 1
    Tref1     = 2000.0  # reference flow stress
    n[1]      = 1.0
    eta0[1]   = Tref1 * (1.0/2.0) * abs(Ebg)^(-1.0/n[1])  # matrix viscosity
    # Material 2
    Tref2     = 2 # reference flow stress
    n[2]      = 1.0
    eta0[2]   = Tref2 * (1.0/2.0) * abs(Ebg)^(-1.0/n[2]) # inclusion viscosity
    # BC definition
    BC = BoundaryConditions()
    if BC_type==1
        ## Pure shear - full Dirichlet
        BC.Vx.type_W = 1;     BC.Vx.Dir_W =-Ebg*xmin*ones(ncy+2)
        BC.Vx.type_E = 1;     BC.Vx.Dir_E =-Ebg*xmax*ones(ncy+2)      
        BC.Vx.type_S =22;     BC.Vx.Dir_S = 0.0*ones(ncx+1)
        BC.Vx.type_N =22;     BC.Vx.Dir_N = 0.0*ones(ncx+1)
        BC.Vy.type_W =22;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
        BC.Vy.type_E =22;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
        BC.Vy.type_S = 1;     BC.Vy.Dir_S = Ebg*ymin*ones(ncx+2)
        BC.Vy.type_N = 1;     BC.Vy.Dir_N = Ebg*ymax*ones(ncx+2)
    elseif BC_type==2
        ## Simple shear 
        BC.Vx.type_W = 1;     BC.Vx.Dir_W = 0.0*ones(ncy+2)
        BC.Vx.type_E = 1;     BC.Vx.Dir_E =-0.0*ones(ncy+2)      
        BC.Vx.type_S =11;     BC.Vx.Dir_S =-1.0*ones(ncx+1)
        BC.Vx.type_N =11;     BC.Vx.Dir_N = 1.0*ones(ncx+1)
        BC.Vy.type_W =22;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
        BC.Vy.type_E =22;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
        BC.Vy.type_S = 1;     BC.Vy.Dir_S =-0.0*ones(ncx+2)
        BC.Vy.type_N = 1;     BC.Vy.Dir_N = 0.0*ones(ncx+2)
    elseif BC_type==3
        ## Simple shear periodic in x
        BC.periodix  = 1
        BC.Vx.type_W = 0;     BC.Vx.Dir_W = 0.0*ones(ncy+2)
        BC.Vx.type_E = 0;     BC.Vx.Dir_E =-0.0*ones(ncy+2)      
        BC.Vx.type_S =11;     BC.Vx.Dir_S =-1.0*ones(ncx+1)
        BC.Vx.type_N =11;     BC.Vx.Dir_N = 1.0*ones(ncx+1)
        BC.Vy.type_W = 0;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
        BC.Vy.type_E = 0;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
        BC.Vy.type_S = 1;     BC.Vy.Dir_S =-0.0*ones(ncx+2)
        BC.Vy.type_N = 1;     BC.Vy.Dir_N = 0.0*ones(ncx+2)
    end
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
    cellxm[1:nmark0] = zeros(Int64,   size(xmi)) #zeros(CartesianIndex{2}, size(xm))
    cellym[1:nmark0] = zeros(Int64,   size(xmi))
    p      = Markers( xm, ym, Tm, phm, cellxm, cellym, nmark0, nmark_max )
    mpc           = zeros(Float64,(ncx  ,ncy  )) # markers per cell
    mpc_th        = zeros(Float64,(nthreads(), ncx  ,ncy   )) # markers per cell
    phase_perc    = zeros(Float64,(nphase, ncx  ,ncy  )) # markers per cell
    phase_perc_th = zeros(Float64,(nthreads(), nphase, ncx  ,ncy   )) # markers per cell
    # Initial configuration
    SetMarkers!( p, rad, xmin, xmax, ymin, ymax, dx, dy )  # Define phase on markers
    xc2   = repeat(xc, 1, length(yc))
    yc2   = repeat(yc, 1, length(xc))'
    xvx2  = repeat(xv, 1, length(yce))
    yvx2  = repeat(yce, 1, length(xv))'
    xvy2  = repeat(xce, 1, length(yv))
    yvy2  = repeat(yv, 1, length(xce))'
    Lx, Ly = xmax - xmin, ymax - ymin
    ax   = 1.0/2.0; ay = 1.0/2.0
    SetInitialVelocity!( Vx, Vy, BC, xv, yv, xce, yce, xmin, xmax, ymin, ymax, ncx, ncy )
    Vxc    = 0.5.*(Vx[1:end-1,2:end-1] .+ Vx[2:end-0,2:end-1])
    Vyc    = 0.5.*(Vy[2:end-1,1:end-1] .+ Vy[2:end-1,2:end-0])
    # Visualisation
    viz_directory = string( "Figures_", experiment )
    if isdir( viz_directory ) == false 
        mkdir( viz_directory ) 
    end 
    path = string( "./", viz_directory, "/" ) 
    anim = Plots.Animation( path, String[] ) 
    # TIME LOOP
    for it=1:nt
        @printf("-------------------------------------\n")
        @printf("---------- Time step #%04d ----------\n", it)
        @printf("-------------------------------------\n")
        # Count markers and reseed
        @time CountMarkers2!(p,mpc,mpc_th,phase_perc,phase_perc_th,nphase,nmark,xmin,xmax,ymin,ymax,xc,yc,dx,dy,ncx,ncy)
        @time ReseedMarkers!( p, mpc, xc, yc, dx, dy, ncx, ncy)
        nmark = p.nmark
        @time CountMarkers2!(p,mpc,mpc_th,phase_perc,phase_perc_th,nphase,nmark,xmin,xmax,ymin,ymax,xc,yc,dx,dy,ncx,ncy)
        # Define dt
        maxVx = maximum(abs.(Vx))
        maxVy = maximum(abs.(Vy))
        dt    = Courant*(maximum([dx,dy])/maximum([maxVx,maxVy]))
        @printf("max V = %2.2e --> dt = %2.2e\n", maximum([maxVx,maxVy]),dt)
        # Update cell info on markers
        @time LocateMarkers(p,dx,dy,xc,yc,xmin,xmax,ymin,ymax)
        # Interpolate k from markers
        # etac = zeros( ncx  , ncy  )
        # @time Markers2Cells3!(p,etac,xc,yc,dx,dy,ncx,ncy,[eta1,eta2],0,0)
        # Initialize
        # @time  etav = CentroidsToVertices( etac, ncx, ncy, BC )
        # Non-linear iterations
        for iter=1:niter_nl

            @printf("-------------------------------------\n")
            @printf("---------- Iteration #%04d ----------\n", iter)
            @printf("-------------------------------------\n")
            
            # Evaluate residuals
            println("Residuals")
            @time SetBCs( Vx, Vy, BC )
            @time StrainRate!(div, Eps, Vx, Vy, dx, dy)
            @time Rheology!( etac, etav, Eps, eta0, n, phase_perc, nphase, BC, ncx, ncy )
            @time Stress!(Tau, Eps, etac, etav)
            @time Residuals!(Fx, Fy, Fp, Tau, Pc, div, BC, ncx , ncy, dx, dy)
            println("|Fx| = ", norm(Fx)/length(Fx))
            println("|Fy| = ", norm(Fy)/length(Fy))
            println("|Fp| = ", norm(Fp)/length(Fp))
            if ( (norm(Fx)/length(Fx)<tol_nl) & (norm(Fy)/length(Fy)<tol_nl) )
                println("Mechanical solver converged")
                break
            end

            # Number equations
            println("Numbering")
            @time NumVx, NumVy, NumP = NumberingStokes(BC, ncx, ncy)
            # Assemble Stokes matrices
            println("Assembly")
            @time Kuu, Kup, Kpu = StokesAssembly( BC, NumVx, NumVy, NumP, etac, etav, Dir_scale, dx, dy )
            # Call solver
            println("Solver")
            @time StokesSolver!(Vx,Vy,Pc,NumVx,NumVy,NumP, Fx, Fy, Fp,Kuu,Kup,Kpu,etac,gamma,solver)
        
            # # Evaluate residuals
            # println("Residuals")
            # @time SetBCs( Vx, Vy, BC )
            # @time StrainRate!(div, Eps, Vx, Vy, dx, dy)
            # @time Stress!(Tau, Eps, etac, etav)
            # @time Residuals!(Fx, Fy, Fp, Tau, Pc, div, BC, ncx , ncy, dx, dy)
            # println("|Fx| = ", norm(Fx)/length(Fx))
            # println("|Fy| = ", norm(Fy)/length(Fy))
            # println("|Fp| = ", norm(Fp)/length(Fp))
        end

        # Advect particles
        @time RungeKutta!(p, nmark, rkv, rkw, BC, dt, Vx, Vy, xv, yv, xce, yce, dx, dy, ncx, ncy)
        # Deform box
        if PureShear_ALE==1
            xmin        += mean(BC.Vx.Dir_W) * dt
            xmax        += mean(BC.Vx.Dir_E) * dt
            ymin        += mean(BC.Vy.Dir_S) * dt
            ymax        += mean(BC.Vy.Dir_N) * dt            
            BC.Vx.Dir_W .=-Ebg*xmin*ones(ncy+2)
            BC.Vx.Dir_E .=-Ebg*xmax*ones(ncy+2)            
            BC.Vy.Dir_S .= Ebg*ymin*ones(ncx+2)
            BC.Vy.Dir_N .= Ebg*ymax*ones(ncx+2)
            L = xmax - xmin
            println("Box pure shear deformation: ", (L-L0)/L0*100, "%" )
            dx, dy, xc, yc, xce, yce, xv, yv = GenerateMesh( xmin, xmax, ymin, ymax, ncx, ncy )
            SetInitialVelocity!( Vx, Vy, BC, xv, yv, xce, yce, xmin, xmax, ymin, ymax, ncx, ncy )
        end
        # Visualisation
        if show_figs==1 && ( mod(it,nout)==0 || it==1 )
            # Visualize
            # p4 = Plots.heatmap(xv*Lc, yce*Lc, Array(Vx)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="Vx")
            # p2 = Plots.heatmap(xce*Lc, yv*Lc, Array(Vy)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="Vy")
            # p3 = Plots.heatmap(xc*Lc,  yc*Lc, Array(Pc)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="P")
            # p4 = Plots.heatmap(xv*Lc,  yv*Lc, Array(etav)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="etav")
            p4 = Plots.heatmap(xc*Lc,  yc*Lc, Array(Tau.II)', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="Tii")
            # p4 = Plots.heatmap(xc*Lc,  yc*Lc, Array(phase_perc[1,:,:])', aspect_ratio=1, xlims=(minimum(xv*Lc), maximum(xv)*Lc), ylims=(minimum(yv)*Lc, maximum(yv)*Lc), c=Plots.cgrad(:roma, rev = true), title="phase %")
            display(Plots.plot( p4, dpi=200 ) ); Plots.frame(anim) 
        end
    end
    Plots.gif(anim, string( path, experiment, ".gif" ), fps = 15)
    return 
    end

   @time main() 