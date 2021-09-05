##############
using Revise
using FDMIC_Geodynamics
using LoopVectorization, Printf, Base.Threads, Plots, Revise, LinearAlgebra, Statistics, SparseArrays
include("./src/EvalAnalDani_v2.jl")
##############
function SetMarkers!( p::Markers, params::ModelParameters, dom::ModelDomain )
    R = params.user[1]
    # Use this function to set up the model geometry
    for k=1:p.nmark 
        in         = (p.x[k]^2 + p.y[k]^2) < R^2 
        p.phase[k] = (in==0)* 1.0 + (in==1)*2.0
    end
end
##############
function Rheology!( f, Eps, Tau, params, materials, BC, ncx, ncy )
    # Compute effective viscosity for each phase
    f.etac .= 0.0
    for m=1:materials.nphase
        eta      = materials.eta0[m] * Eps.II.^(1.0/materials.n[m] - 1.0)
        f.etac .+= f.phase_perc[m,:,:] .* eta
    end
    @time CentroidsToVertices!( f.etav, f.etac, ncx, ncy, BC )
    println("min Eii: ", minimum(Eps.II), " --- max Eii: ", maximum(Eps.II))
    println("min eta: ", minimum(f.etac), " --- max eta: ", maximum(f.etac))
end
export Rheology!
##############
@views function main( N )
    # Main routine
    println("#####################")
    println("###### M2Di.jl ######")
    println("#####################")
    # Name of the experiement
    experiment       = "PeriodicSimpleShear" 
    # Domain
    domain           = ModelDomain()
    domain.xmin      = -1.0
    domain.xmax      =  1.0
    domain.ymin      = -1.0
    domain.ymax      =  1.0
    # parameters
    params           = ModelParameters()
    params.Ebg       = -1.0      # reference strain rate (>0 horizontal compression)
    # Scales
    scale            = CharScales()
    scale.L          = 1.0
    scale.T          = 1.0
    scale.t          = 1.0/abs(params.Ebg)       
    # Spatial discretisation
    ncx              = N; domain.ncx = ncx
    ncy              = N; domain.ncy = ncy
    # Time discretisation
    params.advection = 0
    params.nt        = 1
    params.dt        = 1e4
    params.Courant   = 0.40  # Courant number
    params.dt_var    = 1     # variable dt
    params.t         = 0.0   # time
    # Boundary conditions
    params.BC_type   = 3     # 1: Pure shear / 2: Simple shear / 3: Simple shear periodic
    params.ALE       = 0     # Deform box 
    params.comp      = 0     # Compressible
    # Solver
    #------------------------ PINNED PRESSURE
    params.solver    =-1     # -1: coupled with pinned P, 0: coupled with slight compressibility --- 1: decoupled Powell-Hestenes --- 2: decoupled KSP GCR
    #------------------------ PINNED PRESSURE
    params.gamma     = 1e4   # Penalty factor
    params.Dir_scale = 1.0   # Dirichlet scaling factor
    # Non-linear iterations 
    params.niter_nl  = 5     # max. number of non-linear iterations
    params.tol_nl    = 1e-7  # non-linear tolerance
    params.JFNK      = 0     # Jacobian-Free Newton_Krylov
    # Visualisation
    params.show_figs = 1     # activates visualisation...
    params.nout      = 10    # ... every nout
    # Initialize coordinates
    domain.L0 = domain.xmax-domain.xmin
    GenerateMesh!( domain )
    dx, dy = domain.dx, domain.dy
    # Material setup
    params.user[1]        = 1.0/3.0  # layer half thickness
    # Material parameters
    materials        = MaterialParameters()  
    materials.nphase = 2  
    materials.Tref   = zeros(Float64, materials.nphase)
    materials.eta0   = zeros(Float64, materials.nphase)
    materials.n      = zeros(Float64, materials.nphase)
    # Material 1
    materials.Tref[1]   = 2.0  # reference flow stress
    materials.n[1]      = 1.0
    materials.eta0[1]   = materials.Tref[1] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[1])  # matrix viscosity
    # Material 2
    materials.Tref[2]   = 2000 # reference flow stress
    materials.n[2]      = 1.0
    materials.eta0[2]   = materials.Tref[2] * (1.0/2.0) * abs(params.Ebg)^(-1.0/materials.n[2]) # inclusion viscosity
    # Output arrays for monitoring
    Pana = zeros(Float64, params.nt)
    Pnum = zeros(Float64, params.nt)
    # BC definition
    BC = BoundaryConditions()
    # Pure shear - full Dirichlet
    BC.Vx.type_W = 1;     BC.Vx.Dir_W = 1.0*ones(ncy+2)
    BC.Vx.type_E = 1;     BC.Vx.Dir_E =-1.0*ones(ncy+2)      
    BC.Vx.type_S =11;     BC.Vx.Dir_S = 0.0*ones(ncx+1)
    BC.Vx.type_N =11;     BC.Vx.Dir_N = 0.0*ones(ncx+1)
    BC.Vy.type_W =11;     BC.Vy.Dir_W = 0.0*ones(ncy+1)
    BC.Vy.type_E =11;     BC.Vy.Dir_E = 0.0*ones(ncy+1)       
    BC.Vy.type_S = 1;     BC.Vy.Dir_S =-1.0*ones(ncx+2)
    BC.Vy.type_N = 1;     BC.Vy.Dir_N = 1.0*ones(ncx+2)
    # Evaluate BC's: West/East Vx
    for j=1:ncy+2
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xmin, domain.yce[j], 0)
        BC.Vx.Dir_W[j] = vxa
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xmax, domain.yce[j], 0)
        BC.Vx.Dir_E[j] = vxa
    end
    # Evaluate BC's: South/North Vx
    for i=1:ncx+1
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xv[i], domain.ymin, 0)
        BC.Vx.Dir_S[i] = vxa
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xv[i], domain.ymax, 0)
        BC.Vx.Dir_N[i] = vxa
    end
    # Evaluate BC's: West/East Vy
    for j=1:ncy+1
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xmin, domain.yv[j], 0)
        BC.Vy.Dir_W[j] = vya
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xmax, domain.yv[j], 0)
        BC.Vy.Dir_E[j] = vya
    end
    # Evaluate BC's: South/North Vy
    for i=1:ncx+2
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xce[i], domain.ymin, 0)
        BC.Vy.Dir_S[i] = vya
        vxa, vya = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, domain.xce[i], domain.ymax, 0)
        BC.Vy.Dir_N[i] = vya
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
    # Postprocessing
    Pca       = zeros(Float64, ncx+0, ncy+0)
    Vxca      = zeros(Float64, ncx+0, ncy+0)
    Vyca      = zeros(Float64, ncx+0, ncy+0)
    #------------------------ PINNED PRESSURE
    # Evaluate analytical solution on centroids (moved it before such that pinned pressure can be taken from Pca)
    for i=1:size(Pca,1)
        for j=1:size(Pca,2)
            x  = domain.xc[i]
            y  = domain.yc[j]
            in = sqrt(x^2.0 + y^2.0) <= (params.user[1])
            pa        = SolutionFields_p(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, x, y, in)
            vxa, vya  = SolutionFields_v(materials.eta0[1], materials.eta0[2], params.user[1], params.Ebg, 0.0, x, y, in)
            Pca[i,j]  = pa
            Vxca[i,j] = vxa
            Vyca[i,j] = vya
        end
    end
    # Set pinned pressure
    fields.Pc[1,1] = Pca[1,1]
    #------------------------ PINNED PRESSURE
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
    end
    # Errors
    Vxc   = 0.5*(fields.Vx[1:end-1,:2:end-1] .+ fields.Vx[2:end-0,:2:end-1])
    Vyc   = 0.5*(fields.Vy[2:end-1,:1:end-1] .+ fields.Vy[2:end-1,:2:end-0])
    errV  = sqrt.((Vxc .- Vxca).^2 .+ (Vyc .- Vyca).^2) #.- sqrt.(Vxca.^2 .+ Vyca.^2)
    errP  = abs.(fields.Pc.-Pca)
    println("L1 error V :", mean(errV))
    println("L1 error P :", mean(errP ))
    # Visualize
    if params.show_figs==1
        p1 = Plots.heatmap(domain.xc, domain.yc, Array(Vxc .- Vxca)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Vx")
        p2 = Plots.heatmap(domain.xc, domain.yc, Array(Vyc .- Vyca)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Vy")
        p3 = Plots.heatmap(domain.xc, domain.yc, Array(errP)',        aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="P")
        p4 = Plots.heatmap(domain.xv, domain.yv, Array(fields.etav)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="etav")
        display(Plots.plot( p1, p2, p3, p4, dpi=200 ) ) 
    end
    return mean(errV), mean(errP)
end

####################################################
####################################################
####################################################
N        = [50; 100; 200; 400;]
L1_errV  = zeros(length(N))
L1_errP  = zeros(length(N))
h        = 2.0 ./ N
# Call solver
for ires=1:length(N)
    @time eV, eP = main( N[ires] )
    L1_errV[ires]     = eV
    L1_errP[ires]     = eP
end
# Visualise
p = Plots.plot(  log10.(1.0./h),  log10.(L1_errV), markershape=:rect, color=:blue, label="V", xlabel="1/h" )
p = Plots.plot!( log10.(1.0./h),  log10.(L1_errP), markershape=:rect, color=:red,  label="P", ylabel="L1 error" )
display( p )
display(L1_errV[1:end-1]./L1_errV[2:end-0])
display(L1_errP[1:end-1]./L1_errP[2:end-0])