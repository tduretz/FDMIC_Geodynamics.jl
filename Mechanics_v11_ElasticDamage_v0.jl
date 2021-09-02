##############
using Revise
using FDMIC_Geodynamics
using LoopVectorization, Printf, Base.Threads, Plots, Revise, LinearAlgebra, Statistics, SparseArrays
##############
function SetMarkers!( p::Markers, params::ModelParameters, domain::ModelDomain )    # Use this function to set up the model geometry
    R  = params.user[1]
    L  = domain.xmax - domain.xmin
    H  = domain.ymax - domain.ymin
    x1 = -1.0
    x2 =  1.0
    for k=1:p.nmark #@tturbo <---------rand does not work in @tturbo!!!!
        notch      = ( p.x[k] < 0.0) && (abs(p.y[k] - H/2.0) < R)
        in         = notch
        p.phase[k] = (in==0)* 1.0 + (in==1)*2.0
    end
    return
end
##############
function Rheology!( f, Eps, Tau, params, materials, BC, ncx, ncy )
    # Compute effective viscosity for each phase
    f.etac .= 0.0
    f.Kc   .= 0.0
    f.Gc   .= 0.0
    f.GDam .= 0.0
    f.lDam .= 0.0
    for m=1:materials.nphase
        @tturbo eta    = materials.G[m]*params.dt
        @tturbo f.etac .+= f.phase_perc[m,:,:] .* eta
        @tturbo f.Kc   .+= f.phase_perc[m,:,:] .* materials.K[m]
        @tturbo f.Gc   .+= f.phase_perc[m,:,:] .* materials.G[m]
        @tturbo f.GDam .+= f.phase_perc[m,:,:] .* materials.Gc[m]
        @tturbo f.lDam .+= f.phase_perc[m,:,:] .* materials.l0[m]
    end
    # Damage
    f.etac .*= f.Damc
    f.Kc   .*= f.Damc
    @time CentroidsToVertices!( f.etav, f.etac, ncx, ncy, BC )
    @time CentroidsToVertices!( f.Gv,   f.Gc,   ncx, ncy, BC )
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
    domain.xmin      = -0.5e-3
    domain.xmax      =  0.5e-3
    domain.ymin      = -0.0e-3
    domain.ymax      =  1.0e-3
    # Parameters
    params           = ModelParameters()
    params.Ebg       =  2.0e-5      # reference strain rate (>0 horizontal compression)
    # Scales
    scale            = CharScales()
    scale.L          = 1e-7
    scale.T          = 1.0
    scale.t          = 1.0/abs(params.Ebg)  
    scale.S          = 1e6     
    # Spatial discretisation
    ncx              = 103; domain.ncx = ncx
    ncy              = 100;  domain.ncy = ncy
    # Time discretisation
    params.advection = 0
    params.nt        = 100
    params.dt        = 5e0
    params.Courant   = 0.01  # Courant number
    params.dt_var    = 0     # variable dt
    params.t         = 0.0   # time
    # Boundary conditions
    params.BC_type   = 1     # 1: Pure shear / 2: Simple shear / 3: Simple shear periodic
    params.ALE       = 1     # Deform box 
    params.comp      = 1     # Compressible
    # Solver
    params.solver    = 1     # 0: coupled --- 1: decoupled 
    params.gamma     = 1e8   # Penalty factor
    params.Dir_scale = 1.0   # Dirichlet scaling factor
    # Non-linear iterations 
    params.niter_nl  = 3     # max. number of non-linear iterations
    params.tol_nl    = 1e-7  # non-linear tolerance
    params.JFNK      = 0     # Jacobian-Free Newton_Krylov
    # Visualisation
    params.show_figs = 1     # activates visualisation...
    params.nout      = 10    # ... every nout
    # Material setup
    params.user[1]    = 15e-6/2.0 # layer half thickness
    params.user[2]    = 1e-3 # numerical parameter for phase-field algorithm
    # Material parameters
    materials         = MaterialParameters()  
    materials.nphase  = 2
    materials.Tref    = zeros(Float64, materials.nphase)
    materials.eta0    = zeros(Float64, materials.nphase)
    materials.n       = zeros(Float64, materials.nphase)
    materials.G       = zeros(Float64, materials.nphase)
    materials.K       = zeros(Float64, materials.nphase)
    materials.Gc      = zeros(Float64, materials.nphase)
    materials.l0      = zeros(Float64, materials.nphase)
    # Material 1
    materials.Tref[1] = 1e-1  # reference flow stress
    materials.n[1]    = 1.0
    materials.K[1]    = 174.9e9
    materials.G[1]    = 80.77e9
    materials.Gc[1]   = 2700.0
    materials.l0[1]   = 15e-6
    materials.eta0[1] = materials.Tref[1] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[1])  # matrix viscosity
    # Material 2
    materials.Tref[2] = 1     # reference flow stress
    materials.n[2]    = 1.0
    materials.K[2]    = 174.9e8
    materials.G[2]    = 80.77e8
    materials.Gc[2]   = 2700.0
    materials.l0[2]   = 15e-6
    materials.eta0[2] = materials.Tref[2] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[2]) # inclusion viscosity
   
   domain.xmin /= scale.L
   domain.xmax /= scale.L
   domain.ymin /= scale.L
   domain.ymax /= scale.L
   materials.Gc ./= (scale.S*scale.L) 
   materials.G  ./= scale.S 
   materials.K  ./= scale.S 
   materials.l0 ./=  scale.L
   params.Ebg  /= (1.0/scale.t)
   params.dt   /= scale.t
   params.user[1] /=  scale.L

   displ = 0
   
    # Initialize coordinates
    domain.L0 = domain.xmax-domain.xmin
    GenerateMesh!( domain )
    dx, dy = domain.dx, domain.dy   
    # BC definition
    BC = BoundaryConditions()
    if params.BC_type==1
        ## Pure shear
        BC.Vx.type_W = 1;     BC.Vx.Dir_W =-0.0*params.Ebg*domain.xmin*ones(ncy+2)
        BC.Vx.type_E = 1;     BC.Vx.Dir_E =-0.0*params.Ebg*domain.xmax*ones(ncy+2)      
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
    BC_Dam = Thermal_BC()
    BC_Dam.type_W = 2;     BC_Dam.Dir_W = 0.0*ones(ncy)
    BC_Dam.type_E = 2;     BC_Dam.Dir_E = 0.0*ones(ncy)       
    BC_Dam.type_S = 2;     BC_Dam.Dir_S = 0.0*ones(ncx)
    BC_Dam.type_N = 2;     BC_Dam.Dir_N = 0.0*ones(ncx) 
    # Allocate tables
    fields           = Fields2D()
    fields.etac      = zeros(Float64, ncx  , ncy  )
    fields.etav      = zeros(Float64, ncx+1, ncy+1)
    fields.Kc        = zeros(Float64, ncx  , ncy  )
    fields.Gc        = zeros(Float64, ncx  , ncy  )
    fields.Gv        = zeros(Float64, ncx+1, ncy+1)
    fields.We        = zeros(Float64, ncx  , ncy  )
    fields.rhoc      = zeros(Float64, ncx  , ncy  )
    fields.Pc0       = zeros(Float64, ncx+0, ncy+0)
    fields.Pc        = zeros(Float64, ncx+0, ncy+0)
    fields.Vx        = zeros(Float64, ncx+1, ncy+2) # !!! GHOST ROWS
    fields.Vy        = zeros(Float64, ncx+2, ncy+1) # !!! GHOST COLUMNS
    Tau0             = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Tau              = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Eps              = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Str0             = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    Str              = Tensor2D( zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+1, ncy+1), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0), zeros(Float64, ncx+0, ncy+0) ) 
    if BC.periodix==0 
        fields.Fx    = zeros(Float64, ncx+1, ncy+0)
    else
        fields.Fx    = zeros(Float64, ncx+0, ncy+0)
    end
    fields.Fy        = zeros(Float64, ncx+0, ncy+1)
    fields.Fp        = zeros(Float64, ncx  , ncy  )
    fields.Damc      =  ones(Float64, ncx  , ncy  )
    fields.phiDam    = zeros(Float64, ncx, ncy)
    fields.GDam      = zeros(Float64, ncx  , ncy  )
    fields.lDam      = zeros(Float64, ncx  , ncy  )
    fields.FDam      =  ones(Float64, ncx, ncy)
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
    # Storage
    Wev    = zeros(Float64, params.nt)
    Force  = zeros(Float64, params.nt)
    Displ  = zeros(Float64, params.nt)
    # TIME LOOP
    for it=1:params.nt
        @printf("-------------------------------------\n")
        @printf("---------- Time step #%04d ----------\n", it)
        @printf("-------------------------------------\n")
        # History
        fields.Pc0 .= fields.Pc
        Tau0.xx    .= Tau.xx
        Tau0.yy    .= Tau.yy
        Tau0.zz    .= Tau.zz
        Tau0.xy    .= Tau.xy
        Str0.xx    .= Str.xx
        Str0.yy    .= Str.yy
        Str0.zz    .= Str.zz
        Str0.xy    .= Str.xy
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

            @printf("-------------------------------------------------\n")
            @printf("---------- Iteration #%04d (step #%04d)----------\n", iter, it)
            @printf("-------------------------------------------------\n")
            
            # Evaluate residuals
            println("Residuals")
            @time SetBCs( fields.Vx, fields.Vy, BC )
            @time StrainRate!( Eps, fields.Vx, fields.Vy, domain )
            fields.Damc .= (1.0 .- fields.phiDam ).^2 .+ 1e-6
            @time Rheology!( fields, Eps, Tau, params, materials, BC, ncx, ncy )
            @time StressVE!( Tau, Tau0, Eps, fields, params.dt )
            @time StrainEnergy!( fields.We, Tau, Eps, Str0, Str, fields.Pc, params.dt )

            @time ResidualsComp!( fields, Tau, Eps, BC, domain, params )
            @time DamageResidual!( fields, BC_Dam, domain )
            println("|Fx|   = ", norm(fields.Fx)/length(fields.Fx))
            println("|Fy|   = ", norm(fields.Fy)/length(fields.Fy))
            println("|Fp|   = ", norm(fields.Fp)/length(fields.Fp))
            println("|FDam| = ", norm(fields.FDam)/length(fields.FDam))
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

            K_Dam = DamageAssembly( fields, BC_Dam, domain )
            # dd = K_Dam\fields.FDam[:]
            # fields.phiDam .-= dd[fields.NumP]
            # display(Plots.spy(K_Dam))
            Kc = cholesky(K_Dam)
            dd = Kc\fields.FDam[:]
            fields.phiDam .-= dd[fields.NumP]
            # @show K_Dam - K_Dam'

            # @show fields.FDam
            # @show norm(fields.FDam)
        
            # Evaluate residuals
            println("Residuals")
            @time SetBCs( fields.Vx, fields.Vy, BC )
            @time StrainRate!( Eps, fields.Vx, fields.Vy, domain )
            @time StressVE!( Tau, Tau0, Eps, fields, params.dt )
            @time ResidualsComp!( fields, Tau, Eps, BC, domain, params )
            @time DamageResidual!( fields, BC_Dam, domain )

            println("|Fx|   = ", norm(fields.Fx)/length(fields.Fx))
            println("|Fy|   = ", norm(fields.Fy)/length(fields.Fy))
            println("|Fp|   = ", norm(fields.Fp)/length(fields.Fp))
            println("|FDam| = ", norm(fields.FDam)/length(fields.FDam))
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
        Wev[it]   = fields.Kc[1,1]#mean(Tau.II)
        Wev[it]   = mean(Tau.II)
        Force[it] = sum( Tau.yy[:,end]*dx .- fields.Pc[:,end]*dx) * (domain.ymax-domain.ymin)
        displ    += maximum(abs.(fields.Vy).*params.dt)
        Displ[it] = displ
        println("Displacement: ", displ * scale.L * 1e3, " --- Max incremental displacement: ", maximum((fields.Vy).*params.dt * scale.L * 1e3 ),  " mm")


        # Visualisation
        if params.show_figs==1 && ( mod( it, params.nout )==0 || it==1 )
            l = @layout [ [a; b] c]
            # p1 = Plots.heatmap( domain.xv, domain.yce, Array(fields.Vx*scale.L/scale.t)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="K" )
            # p2 = Plots.heatmap( domain.xce, domain.yv, Array(fields.Vy*scale.L/scale.t)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="phi" )
            # p1 = Plots.heatmap( domain.xc, domain.yc, Array(fields.We)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="We" )
            p1 = Plots.heatmap( domain.xc, domain.yc, Array(Tau.II)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="str" )
            p2 = Plots.heatmap( domain.xc, domain.yc, Array(fields.phiDam)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="phi" )
            # p3 = Plots.heatmap( domain.xc, domain.yc, Array(Tau.II*scale.S)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Tii" )
            # p3 = Plots.plot( Displ[1:it]*scale.L*1e3, Force[1:it]*(scale.S*scale.L^2)/1e3, markershape=:circle)
            p3 = Plots.plot( Displ[1:it]*scale.L*1e3, Wev[1:it]*scale.S/1e9, markershape=:circle)

            display( Plots.plot( p1, p2, p3, dpi=200, layout=l ) ); Plots.frame(anim) 
        end
    end
    if params.show_figs==1 Plots.gif(anim, string( path, experiment, ".gif" ), fps = 15) end
    return 
end

for it=1:1
   @time main()
end
