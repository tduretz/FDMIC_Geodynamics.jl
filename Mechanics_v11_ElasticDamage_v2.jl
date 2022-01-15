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
        p.phase[k] = 1.0
    end
    # for k=1:p.nmark #@tturbo <---------rand does not work in @tturbo!!!!
    #     notch      = ( p.x[k] < 0.0) && (abs(p.y[k] - H/2.0) < R)
    #     in         = notch
    #     p.phase[k] = (in==0)* 1.0 + (in==1)*2.0
    # end
    return
end
##############
function Rheology!( f, Eps, Tau, params, materials, BC, ncx, ncy )
    # Compute effective viscosity for each phase
    f.etac .= 0.0
    f.Kc   .= 0.0
    f.Gc   .= 0.0
    # f.GDam .= 0.0
    # f.lDam .= 0.0
    for m=1:materials.nphase
        @tturbo eta    = materials.G[m]*params.dt
        @tturbo f.etac .+= f.phase_perc[m,:,:] .* eta           
        @tturbo f.Kc   .+= f.phase_perc[m,:,:] .* materials.K[m]
        @tturbo f.Gc   .+= f.phase_perc[m,:,:] .* materials.G[m] 
        # @tturbo f.GDam .+= f.phase_perc[m,:,:] .* materials.Gc[m] 
        # @tturbo f.lDam .+= f.phase_perc[m,:,:] .* materials.l0[m]
    end
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
    # Name of the experiment
    experiment       = "ElasticDamageDyke" 
    # Domain
    domain           = ModelDomain()
    domain.xmin      = -2.5e3
    domain.xmax      =  2.5e3
    domain.ymin      = -5.0e3
    domain.ymax      =  0.0e-3
    # Parameters
    params           = ModelParameters()
    params.Ebg       = -2.0e-14      # reference strain rate (>0 horizontal compression)
    # Scales
    scale            = CharScales()
    scale.L          = 1e3
    scale.T          = 1.0
    scale.t          = 1.0/abs(params.Ebg)  
    scale.S          = 1e8   
    scale.m          = scale.S*scale.L*scale.t^2 
    scale.rho        = scale.S*scale.L^(-2)*scale.t^2
    # Spatial discretisation
    ncx              = 150;  domain.ncx = ncx
    ncy              = 150;  domain.ncy = ncy
    # Time discretisation
    params.advection = 0
    params.nt        = 2#6*150
    params.dt        = 1.25e8/4
    params.Courant   = 0.01  # Courant number
    params.dt_var    = 0     # variable dt
    params.t         = 0.0   # time
    # Boundary conditions
    params.BC_type   = 1     # 1: Pure shear / 2: Simple shear / 3: Simple shear periodic
    params.ALE       = 1     # Deform box 
    params.comp      = 1     # Compressible
    # Solver
    params.solver    = 1     # 0: coupled --- 1: decoupled 
    params.gamma     = 1e7   # Penalty factor
    params.Dir_scale = 1.0   # Dirichlet scaling factor
    # Non-linear iterations 
    params.niter_nl  = 15    # max. number of non-linear iterations
    params.tol_nl    = 1e-4  # non-linear tolerance
    params.JFNK      = 0     # Jacobian-Free Newton_Krylov
    # Visualisation
    params.show_figs = 1     # activates visualisation...
    params.nout      = 10    # ... every nout
    # Material setup
    params.user[1]    = 5e2 # layer half thickness
    params.user[2]    = 1e-3 # numerical parameter for phase-field algorithm
    # Material parameters
    materials         = MaterialParameters()  
    materials.nphase  = 2
    materials.Tref    = zeros(Float64, materials.nphase)
    materials.eta0    = zeros(Float64, materials.nphase)
    materials.n       = zeros(Float64, materials.nphase)
    materials.G       = zeros(Float64, materials.nphase)
    materials.K       = zeros(Float64, materials.nphase)
    # Constant material parameters (all phases)
    materials.gDam    = 2.70e3#zeros(Float64, materials.nphase)
    materials.lDam    = 50.0#zeros(Float64, materials.nphase)
    materials.eDam    = 1.0e20#zeros(Float64, materials.nphase)
    # Material 1
    materials.Tref[1] = 1e-1  # reference flow stress
    materials.n[1]    = 1.0
    materials.K[1]    = 174.9e9
    materials.G[1]    = 80.77e9
    materials.eta0[1] = materials.Tref[1] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[1])  # matrix viscosity
    # Material 2
    materials.Tref[2] = 1     # reference flow stress
    materials.n[2]    = 1.0
    materials.K[2]    = 174.9e8
    materials.G[2]    = 80.77e8
    materials.eta0[2] = materials.Tref[2] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[2]) # inclusion viscosity
    # Scaling
    domain.xmin     /= scale.L
    domain.xmax     /= scale.L
    domain.ymin     /= scale.L
    domain.ymax     /= scale.L
    materials.gDam /= (scale.S*scale.L) 
    materials.eDam /= (scale.S*scale.t) 
    materials.lDam /=  scale.L
    materials.G    ./= scale.S 
    materials.K    ./= scale.S 
    params.Ebg      /= (1.0/scale.t)
    params.dt       /= scale.t
    params.user[1]  /=  scale.L
    # Initialize coordinates
    domain.L0 = domain.xmax-domain.xmin
    GenerateMesh!( domain )
    dx, dy = domain.dx, domain.dy   
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
        BC.Vy.type_S = 1;     BC.Vy.Dir_S = 0.0*params.Ebg*domain.ymin*ones(ncx+2)
        BC.Vy.type_N = 1;     BC.Vy.Dir_N = 0.0*params.Ebg*domain.ymax*ones(ncx+2)
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
    BC_Dam.Pinned = zeros(Int64,ncx,ncy) 
    # Allocate tables
    fields           = Fields2D()
    fields.etac      = zeros(Float64, ncx  , ncy  )
    fields.etav      = zeros(Float64, ncx+1, ncy+1)
    fields.Kc        = zeros(Float64, ncx  , ncy  )
    fields.Gc        = zeros(Float64, ncx  , ncy  )
    fields.Gv        = zeros(Float64, ncx+1, ncy+1)
    fields.We0       = zeros(Float64, ncx  , ncy  )
    fields.We        = zeros(Float64, ncx  , ncy  ) 
    fields.rhoc      = zeros(Float64, ncx  , ncy  )
    fields.rhov      = zeros(Float64, ncx  , ncy  )
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
    fields.Damc      = zeros(Float64, ncx  , ncy  )
    fields.Damv      = zeros(Float64, ncx+1, ncy+1)
    fields.phiDam    = zeros(Float64, ncx,   ncy  )
    fields.phiDam0   = zeros(Float64, ncx,   ncy  )
    fields.FDam      = zeros(Float64, ncx,   ncy  )
    # Number equations
    println("Numbering")
    @time NumberingStokes!( fields, BC, domain )
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
    
    drho         = 200.0  / scale.rho
    rho0         = 2700.0 / scale.rho
    fields.rhoc .= rho0 # mc      = tauc*Lc*tc^2;  
    params.gy    = -0*9.81   / (scale.L/scale.t^2)

    # Initial phase field
    BC_Dam.Pinned = zeros(Int64,ncx,ncy) 
    x0 = (domain.xmax+domain.xmin)/2.0
    y0 = domain.ymin+0.2*(domain.ymax-domain.ymin)
    for i=1:ncx
        for j=1:ncy
            if  ( (domain.xc[i] - x0)^2 + (domain.yc[j] - y0)^2 ) < params.user[1]^2 
                BC_Dam.Pinned[i,j] = 1
            end
        end
    end
    fields.phiDam[BC_Dam.Pinned.==1] .= 1.0
    @time CountMarkers3!( p, fields, materials, domain )
    @time Rheology!( fields, Eps, Tau, params, materials, BC, ncx, ncy )
    InitialDamageResidual!( fields, BC_Dam, domain, materials )
    println("|FDam| = ", norm(fields.FDam)/length(fields.FDam))
    K_Dam = InitialDamageAssembly( fields, BC_Dam, domain, materials )
    Kc    = cholesky(K_Dam)
    dd    = Kc\fields.FDam[:]
    fields.phiDam .-= dd[fields.NumP]
    InitialDamageResidual!( fields, BC_Dam, domain, materials )
    println("|FDam| = ", norm(fields.FDam)/length(fields.FDam))

    fields.rhoc .-= drho.*fields.phiDam
    @time CentroidsToVertices!( fields.rhov, fields.rhoc, ncx, ncy, BC )

    # Visualisation
    viz_directory = string( "Figures_", experiment )
    if isdir( viz_directory ) == false 
        mkdir( viz_directory ) 
    end 
    path = string( "./", viz_directory, "/" ) 
    anim = Plots.Animation( path, String[] ) 
    # Storage
    displ  = 0
    Wev    = zeros(Float64, params.nt)
    Force  = zeros(Float64, params.nt)
    Displ  = zeros(Float64, params.nt)
    resxv  = zeros(Float64, params.niter_nl)
    # TIME LOOP
    for it=1:params.nt
        @printf("-------------------------------------\n")
        @printf("---------- Time step #%04d ----------\n", it)
        @printf("-------------------------------------\n")
        # History
        fields.We0 .= fields.We
        if it>1 
            fields.phiDam0 .= fields.phiDam 
        end
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
        # Non-linear iterations
        iterations = 0
        for iter=1:params.niter_nl
            iterations+=1
            @printf("-------------------------------------------------\n")
            @printf("---------- Iteration #%04d (step #%04d)----------\n", iter, it)
            @printf("-------------------------------------------------\n")
            # Evaluate residuals
            println("Residuals")
            @time SetBCs( fields.Vx, fields.Vy, BC )
            @time StrainRate!( Eps, fields.Vx, fields.Vy, domain )
            fields.Damc .= (1.0 .- fields.phiDam ).^2 .+ 1e-3
            @time CentroidsToVertices!( fields.Damv, fields.Damc, ncx, ncy, BC )
            fields.rhoc .= rho0 .- drho.*fields.phiDam
            @time CentroidsToVertices!( fields.rhov, fields.rhoc, ncx, ncy, BC )
        
            @time Rheology!( fields, Eps, Tau, params, materials, BC, ncx, ncy )
            @time StressVE!( Tau, Tau0, Eps, fields, params.dt )
            @time StrainEnergy!( fields.We, Tau, Eps, Str0, Str, fields.Pc, params.dt, fields, dx, dy )
            @time ResidualsCompDam!( fields, Tau, Eps, BC, domain, params )
            @time DamageResidual!( fields, BC_Dam, domain, materials, params.dt )
            println("|Fx|   = ", norm(fields.Fx)/length(fields.Fx))
            println("|Fy|   = ", norm(fields.Fy)/length(fields.Fy))
            println("|Fp|   = ", norm(fields.Fp)/length(fields.Fp))
            println("|FDam| = ", norm(fields.FDam)/length(fields.FDam))
            if ( (norm(fields.Fx)/length(fields.Fx) < params.tol_nl) & (norm(fields.Fy)/length(fields.Fy) < params.tol_nl) & (norm(fields.Fp)/length(fields.Fp) < params.tol_nl)  && (norm(fields.FDam)/length(fields.FDam) < params.tol_nl))
                println("Mechanical solver converged")
                break
            end
            resxv[iter] = norm(fields.Fx)/length(fields.Fx)

            # Assemble Stokes matrices
            println("Assembly")
            @time Kuu, Kup, Kpu = StokesAssembly( BC, fields, params.Dir_scale, domain )
            # Call solver
            println("Solver")
            @time StokesSolver!( fields, params, domain, materials, BC, Kuu, Kup, Kpu )

            K_Dam = DamageAssembly( fields, BC_Dam, domain, materials, params.dt )
            # dd = K_Dam\fields.FDam[:]
            # fields.phiDam .-= dd[fields.NumP]
            # display(Plots.spy(K_Dam))
            Kc = cholesky(K_Dam)
            dd = Kc\fields.FDam[:]
            fields.phiDam .-= dd[fields.NumP]

            # Evaluate residuals
            println("Residuals")
            @time SetBCs( fields.Vx, fields.Vy, BC )
            @time StrainRate!( Eps, fields.Vx, fields.Vy, domain )
            @time StressVE!( Tau, Tau0, Eps, fields, params.dt )
            @time ResidualsCompDam!( fields, Tau, Eps, BC, domain, params )
            @time DamageResidual!( fields, BC_Dam, domain, materials, params.dt )

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
        # Wev[it]   = fields.Kc[1,1]#mean(Tau.II)
        Wev[it]   = mean(fields.We)
        Force[it] = sum( Tau.yy[:,end] .- fields.Pc[:,end])
        Force[it]+= sum( Tau.yy[:,1] .- fields.Pc[:,end])
        Force[it]+= sum( Tau.xx[1,:] .- fields.Pc[1,:])
        Force[it]+= sum( Tau.xx[end,:] .- fields.Pc[end,:])
        displ    += maximum(abs.(fields.Vy[:,end]).*params.dt)
        Displ[it] = displ
        
        DIV = maximum(fields.Vy)/(domain.ymax-domain.ymin)
        P   = -mean(fields.Kc)*DIV*params.t
        EXXd = 0.0 - 1.0/3.0*DIV
        EYYd = DIV - 1.0/3.0*DIV
        TXX = 2.0*mean(fields.Gc) * EXXd*params.t
        TYY = 2.0*mean(fields.Gc) * EYYd*params.t
        SXX = -P + TXX
        SYY = -P + TYY
        EYY = DIV*params.t
        WE  = EYY*SYY
        println("Displacement: ", displ * scale.L * 1e3, " --- Max incremental displacement: ", maximum((fields.Vy).*params.dt * scale.L * 1e3 ),  " mm")
        println("dx = ", domain.dx*scale.L*1e3, " mm --- l0 = ",materials.lDam*scale.L*1e3, " mm")
        println("min/max div: ", minimum(Eps.div)*(1.0/scale.t), " ", maximum(Eps.div)*(1.0/scale.t))
        println("bulk divergence: ",  DIV*(1.0/scale.t) )
        println("min/max P: ", minimum(fields.Pc)*scale.S, " ", maximum(fields.Pc)*scale.S)
        println("numerical pressure: ", mean(fields.Pc)*scale.S, " --- expected: ", P*scale.S)
        println("numerical Txx:      ", mean(Tau.xx)*scale.S, " --- expected: ", TXX*scale.S)
        println("numerical Tyy:      ", mean(Tau.yy)*scale.S, " --- expected: ", TYY*scale.S)
        println("numerical Sxx:      ", mean(Tau.xx.-fields.Pc)*scale.S, " --- expected: ", SXX*scale.S)
        println("numerical Syy:      ", mean(Tau.yy.-fields.Pc)*scale.S, " --- expected: ", SYY*scale.S)
        println("numerical We:       ", mean(fields.We)*scale.S, " --- expected: ", WE*scale.S)

        # Visualisation
        if params.show_figs==1 && ( mod( it, params.nout )==0 || it==1 )
            l = @layout [ [a; b] c]
            cmy = 3600*24*365*100
            p2 = Plots.heatmap( domain.xc, domain.yc, (fields.rhoc*scale.rho)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="rho [kg/m^3]",  clim=(2500,2700)  )
            # p2 = Plots.heatmap( domain.xc, domain.yc, (Tau.xy_c .- 0.0*fields.Pc)*scale.S)
            # p2 = Plots.heatmap( domain.xc, domain.yc, Eps.div*(1.0/scale.t) )
            p1 = Plots.heatmap( domain.xv, domain.yce, Array(fields.Vx*cmy*scale.L/scale.t)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Vx [cm/y]" ,  clim=(-0.15,0.15) )
            # p2 = Plots.heatmap( domain.xce, domain.yv, Array(fields.Vy*cmy*scale.L/scale.t)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Vy [cm/y]" )
            # p1 = Plots.heatmap( domain.xc, domain.yc, Array(fields.We)'*scale.S, aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="We" )
            p3 = Plots.heatmap( domain.xc, domain.yc, Array(fields.phiDam)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="phi, niter = $iterations",  clim=(0.0,1.0) )
            p4 = Plots.heatmap( domain.xc, domain.yc, Array(Tau.II)'./1e6*scale.S, aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Sii [MPa], nt = $it",  clim=(0.0,40)  )
            # p4 = Plots.heatmap( domain.xc, domain.yc, Array(Str.II)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Strii" )

            # p3 = Plots.heatmap( domain.xc, domain.yc, Array(Tau.II*scale.S)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Tii" )
            # p3 = Plots.plot( Displ[1:it]*scale.L*1e3, Force[1:it]*(scale.S), markershape=:circle)
            # p3 = Plots.plot( Displ[1:it]*scale.L*1e3, Wev[1:it]*scale.S, markershape=:circle)
            # p3 = Plots.plot( [1:iterations], log10.(resxv[1:iterations]), markershape=:circle)
            display( Plots.plot( p1, p2, p3, p4, dpi=200 ) ); Plots.frame(anim)
            # display( Plots.plot( p1, p2, p3, dpi=200, layout=l ) ); Plots.frame(anim) 
        end
    end
    if params.show_figs==1 Plots.gif(anim, string( path, experiment, ".gif" ), fps = 15) end
    return 
end

for it=1:1
   @time main()
end
