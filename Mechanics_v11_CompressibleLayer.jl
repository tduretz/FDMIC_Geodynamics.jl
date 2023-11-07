##############
using Revise
using FDMIC_Geodynamics
using Printf, Base.Threads, Plots, Revise, LinearAlgebra, Statistics, SparseArrays
##############
function SetMarkers!( p, params, domain )    # Use this function to set up the model geometry
    R  = params.user[1]
    L  = domain.xmax - domain.xmin
    H  = domain.ymax - domain.ymin
    x1 = -1.0
    x2 =  1.0
    for k=1:p.nmark # <---------rand does not work in !!!!
        p.phase[k] = 1
        if abs(p.y[k]) < R
            p.phase[k] = 2.0
        end
    end
end
##############
function Rheology!( f, Eps, Tau, params, materials, BC, ncx, ncy )
    # Compute effective viscosity for each phase
    f.etac .= 0.0
    f.Kc   .= 0.0
    for m=1:materials.nphase
        eta    = 2.0* materials.eta0[m] * Eps.II.^(1.0/materials.n[m] - 1.0)
        f.etac .+= f.phase_perc[m,:,:] .* eta
        f.Kc   .+= f.phase_perc[m,:,:] .* materials.K[m]
    end
    @time CentroidsToVertices!( f.etav, f.etac, ncx, ncy, BC )
    println("min Eii: ", minimum(Eps.II), " --- max Eii: ", maximum(Eps.II))
    println("min eta: ", minimum(f.etac), " --- max eta: ", maximum(f.etac))
end
export Rheology!
##############
@views function main( )
    # Main routine
    println("#####################")
    println("###### M2Di.jl ######")
    println("#####################")
    # Name of the experiement
    experiment       = "CompressibleLayer" 
    # Domain
    domain           = ModelDomain()
    domain.xmin      = -1.0/2.0
    domain.xmax      =  1.0/2.0
    domain.ymin      = -1.0/2.0
    domain.ymax      =  1.0/2.0
    # parameters
    params           = ModelParameters()
    params.Ebg       =  1.0      # reference strain rate (>0 horizontal compression)
    # Scales
    scale            = CharScales()
    scale.L          = 1.0
    scale.T          = 1.0
    scale.t          = 1.0/abs(params.Ebg)       
    # Spatial discretisation
    ncx              = 50; domain.ncx = ncx
    ncy              = 50; domain.ncy = ncy
    # Time discretisation
    params.advection = 0
    params.nt        = 100
    params.dt        = 1e4
    params.Courant   = 0.25  # Courant number
    params.dt_var    = 0     # variable dt
    params.t         = 0.0   # time
    # Boundary conditions
    params.BC_type   = 1     # 1: Pure shear / 2: Simple shear / 3: Simple shear periodic
    params.ALE       = 1     # Deform box 
    params.comp      = 1     # Compressible
    # Solver
    params.solver    = 2     # 0: coupled --- 1: decoupled 
    params.gamma     = 1e8   # Penalty factor
    params.Dir_scale = 1.0   # Dirichlet scaling factor
    # Non-linear iterations 
    params.niter_nl  = 1     # max. number of non-linear iterations
    params.tol_nl    = 1e-7  # non-linear tolerance
    params.JFNK      = 1     # Jacobian-Free Newton_Krylov
    # Visualisation
    params.show_figs = 1     # activates visualisation...
    params.nout      = 5     # ... every nout
    # Initialize coordinates
    domain.L0 = domain.xmax-domain.xmin
    GenerateMesh!( domain )
    dx, dy = domain.dx, domain.dy
    # Material setup
    params.user[1]    = 0.13333333/2 # layer half thickness
    # Material parameters
    materials         = MaterialParameters()  
    materials.nphase  = 2
    materials.Tref    = zeros(Float64, materials.nphase)
    materials.eta0    = zeros(Float64, materials.nphase)
    materials.n       = zeros(Float64, materials.nphase)
    materials.K       = zeros(Float64, materials.nphase)
    # Material 1
    materials.Tref[1] = 1e-1  # reference flow stress
    materials.n[1]    = 1.0
    materials.K[1]    = 1e-5
    materials.eta0[1] = materials.Tref[1] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[1])  # matrix viscosity
    # Material 2
    materials.Tref[2] = 1     # reference flow stress
    materials.n[2]    = 1.0
    materials.K[2]    = 1e-5
    materials.eta0[2] = materials.Tref[2] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[2]) # inclusion viscosity
    # Output arrays for monitoring
    Pana = zeros(Float64, params.nt)
    Pnum = zeros(Float64, params.nt)
    # BC definition
    BC = BoundaryConditions()
    if params.BC_type==1
        ## Pure shear
        BC.Vx.type_W = 1;     BC.Vx.Dir_W =-params.Ebg*domain.xmin*ones(ncy+2)
        BC.Vx.type_E = 1;     BC.Vx.Dir_E =-params.Ebg*domain.xmax*ones(ncy+2)      
        BC.Vx.type_S =22;     BC.Vx.Dir_S = 0.0*ones(ncx+1)
        BC.Vx.type_N =22;     BC.Vx.Dir_N = 0.0*ones(ncx+1)
        BC.Vy.type_W =22;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
        BC.Vy.type_E =22;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
        BC.Vy.type_S = 1;     BC.Vy.Dir_S = params.Ebg*domain.ymin*ones(ncx+2)
        BC.Vy.type_N = 1;     BC.Vy.Dir_N = params.Ebg*domain.ymax*ones(ncx+2)
    elseif params.BC_type==2
        ## Simple shear 
        BC.Vx.type_W = 1;     BC.Vx.Dir_W = 0.0*ones(ncy+2)
        BC.Vx.type_E = 1;     BC.Vx.Dir_E =-0.0*ones(ncy+2)      
        BC.Vx.type_S =11;     BC.Vx.Dir_S =-1.0*ones(ncx+1)
        BC.Vx.type_N =11;     BC.Vx.Dir_N = 1.0*ones(ncx+1)
        BC.Vy.type_W =22;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
        BC.Vy.type_E =22;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
        BC.Vy.type_S = 1;     BC.Vy.Dir_S =-0.0*ones(ncx+2)
        BC.Vy.type_N = 1;     BC.Vy.Dir_N = 0.0*ones(ncx+2)
    elseif params.BC_type==3
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
    fields           = Fields2D()
    fields.etac      = zeros(Float64, ncx  , ncy  )
    fields.etav      = zeros(Float64, ncx+1, ncy+1)
    fields.Kc        = zeros(Float64, ncx  , ncy  )
    fields.rhoc      = zeros(Float64, ncx  , ncy  )
    fields.Pc0       = zeros(Float64, ncx+0, ncy+0)
    fields.Pc        = zeros(Float64, ncx+0, ncy+0)
    fields.Vx        = zeros(Float64, ncx+1, ncy+2) # !!! GHOST ROWS
    fields.Vy        = zeros(Float64, ncx+2, ncy+1) # !!! GHOST COLUMNS
    fields.Damc      = ones(Float64, ncx  , ncy  )
    fields.Damv      = ones(Float64, ncx+1, ncy+1)
    Tau              = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Eps              = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    if BC.periodix==0 
        fields.Fx    = zeros(Float64, ncx+1, ncy+0)
    else
        fields.Fx    = zeros(Float64, ncx+0, ncy+0)
    end
    fields.Fy        = zeros(Float64, ncx+0, ncy+1)
    fields.Fp        = zeros(Float64, ncx  , ncy  )
    # Initialise markers
    p                = Markers()
    p.nmx            = 4            # 2 marker per cell in x
    p.nmy            = 4            # 2 marker per cell in y
    InitialiseMarkerPositions!( p, domain )
    # Over allocate markers
    p.phase          = zeros(Float64, p.nmark_max)
    p.T              = zeros(Float64, p.nmark_max)
    p.cellx          = zeros(Int64,   p.nmark_max)
    p.celly          = zeros(Int64,   p.nmark_max)
    # Marker 2 cell mapping
    fields.mpc           = zeros(Float64,(ncx  ,ncy  )) # markers per cell
    fields.mpc_th        = zeros(Float64,(nthreads(), ncx  ,ncy   )) # markers per cell
    fields.phase_perc    = zeros(Float64,(materials.nphase, ncx  ,ncy  )) # markers per cell
    fields.phase_perc_th = zeros(Float64,(nthreads(), materials.nphase, ncx  ,ncy   )) # markers per cell
    # Initial configuration
    SetMarkers!( p, params, domain )  # Define phase on markers
    SetInitialVelocity2!( fields.Vx, fields.Vy, BC, domain )
    fields.Vxc    = 0.5.*(fields.Vx[1:end-1,2:end-1] .+ fields.Vx[2:end-0,2:end-1])
    fields.Vyc    = 0.5.*(fields.Vy[2:end-1,1:end-1] .+ fields.Vy[2:end-1,2:end-0])
    # Visualisation
    viz_directory = string( "Figures_", experiment )
    if isdir( viz_directory ) == false 
        mkdir( viz_directory ) 
    end 
    path = string( "./", viz_directory, "/" ) 
    anim = Plots.Animation( path, String[] ) 
    # TIME LOOP
    for it=1:params.nt
        @printf("-------------------------------------\n")
        @printf("---------- Time step #%04d ----------\n", it)
        @printf("-------------------------------------\n")
        # History
        fields.Pc0 .= fields.Pc
        # Count markers and reseed
        @time CountMarkers3!( p, fields, materials, domain )
        @time ReseedMarkers2!( p, fields.mpc, domain)
        nmark = p.nmark
        @time CountMarkers3!( p, fields, materials, domain )
        # Define dt
        if params.dt_var == 1
            maxVx     = maximum(abs.(fields.Vx))
            maxVy     = maximum(abs.(fields.Vy))
            params.dt = params.Courant*(maximum([domain.dx,domain.dy])/maximum([maxVx,maxVy]))
            @printf("dt = %2.2e\n", params.dt)
        end
        # Update cell info on markers
        @time LocateMarkers2!(p, domain)
        # Number equations
        println("Numbering")
        @time NumberingStokes!( fields, BC, domain )
        # Non-linear iterations
        for iter=1:params.niter_nl

            @printf("-------------------------------------\n")
            @printf("---------- Iteration #%04d ----------\n", iter)
            @printf("-------------------------------------\n")
            
            # Evaluate residuals
            println("Residuals")
            @time SetBCs( fields.Vx, fields.Vy, BC )
            @time StrainRate!( Eps, fields.Vx, fields.Vy, domain )
            @time Rheology!( fields, Eps, Tau, params, materials, BC, ncx, ncy )
            @time Stress!( Tau, Eps, fields )
            @time ResidualsComp!( fields, Tau, Eps, BC, domain, params )
            println("|Fx| = ", norm(fields.Fx)/length(fields.Fx))
            println("|Fy| = ", norm(fields.Fy)/length(fields.Fy))
            println("|Fp| = ", norm(fields.Fp)/length(fields.Fp))
            if ( (norm(fields.Fx)/length(fields.Fx) < params.tol_nl) & (norm(fields.Fy)/length(fields.Fy) < params.tol_nl) & (norm(fields.Fp)/length(fields.Fp) < params.tol_nl) )
                println("Mechanical solver converged")
                break
            end

            # Assemble Stokes matrices
            println("Assembly")
            @time Kuu, Kup, Kpu = StokesAssembly( BC, fields, params.Dir_scale, domain )
            # Call solver
            println("Solver")
            @time StokesSolver!( fields, params, domain, materials, BC, Kuu, Kup, Kpu )
        
            # Evaluate residuals
            println("Residuals")
            @time SetBCs( fields.Vx, fields.Vy, BC )
            @time StrainRate!( Eps, fields.Vx, fields.Vy, domain )
            @time Stress!( Tau, Eps, fields )
            @time ResidualsComp!( fields, Tau, Eps, BC, domain, params )
            println("|Fx| = ", norm(fields.Fx)/length(fields.Fx))
            println("|Fy| = ", norm(fields.Fy)/length(fields.Fy))
            println("|Fp| = ", norm(fields.Fp)/length(fields.Fp))
        end
        params.t += params.dt
        # Advect particles
        if params.advection == 1
            @time RungeKutta2!( p, fields, params, BC, domain )
            # Deform box
            if params.ALE==1
                dx, dy = PureShearBoxUpdate!( domain, BC, params )
                SetInitialVelocity2!( fields.Vx, fields.Vy, BC, domain )
            end
        end
        # Postprocessing
        tmax     = 2materials.eta0[2]/materials.K[1]
        dmus     = 2*(materials.eta0[2] - materials.eta0[1]) 
        exx_tot  = diff( fields.Vx[:,2:end-1], dims=1 )./dx
        ffs      = abs( 2exx_tot[1,Int64(floor(ncy/2))]* dmus)
        Pana[it] = 1.0 - exp(-3.0/4.0*params.t/tmax) 
        Pnum[it] = ( fields.Pc[1,Int64(floor(ncy/2))] - fields.Pc[1,1] ) / ffs
        println("Pressure in layer num. :", Pnum[it] )
        println("Pressure in layer ana. :", Pana[it] )
        # Visualisation
        if params.show_figs==1 && ( mod( it, params.nout )==0 || it==1 )
            p1 = Plots.heatmap( domain.xc, domain.yc, Array(fields.Pc)'./ffs, aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="P",  clim=(0.0,1.0) )
            # p1 = Plots.heatmap( domain.xc, domain.yc, Array(fields.etac)'./ffs, aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="ηc" )
            p2 = Plots.plot( (1:it)*params.dt/tmax, Pnum[1:it], markershape=:dtriangle, label="P num.", legend=:bottomright)
            p2 = Plots.plot!((1:it)*params.dt/tmax, Pana[1:it], markershape=:circle,    label="P ana.")
            display(Plots.plot( p1, p2, dpi=200 ) ); Plots.frame(anim) 
        end
    end
    if params.show_figs==1 Plots.gif(anim, string( path, experiment, ".gif" ), fps = 15) end
    return 
end

for it=1:1
   @time main()
end