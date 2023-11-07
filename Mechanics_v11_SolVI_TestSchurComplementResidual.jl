##############
using Revise
using FDMIC_Geodynamics
using Printf, Base.Threads, Plots, Revise, LinearAlgebra, Statistics, SparseArrays
include("./src/EvalAnalDani_v2.jl")
##############
function Residual_Stokes_SchurComplement!( Fx, Fy, Vx, Vy, p, γ, bcv, BC )
    Δx=p.Δx; Δy=p.Δy; BCx=p.BCx; BCy=p.BCy; ηv=p.ηv; ηc=p.ηc; ncx=p.ncx; ncy=p.ncy; bx=p.bx; by=p.by; bp=p.bp; nvx=ncx+1; nvy=ncy+1
    
    for j=1:ncy, i=2:nvx-1
        if j==1
            VxS = bcv*2*BC.Vx.Dir_S[i] - Vx[i,j]
        else
            VxS = Vx[i,j-1]
        end
        if j==ncy  
            VxN = bcv*2*BC.Vx.Dir_N[i] - Vx[i,j] 
        else
            VxN = Vx[i,j+1] 
        end
        ∂Vx∂xW  = (Vx[i,j]   - Vx[i-1,j]) / Δx
        ∂Vx∂xE  = (Vx[i+1,j] - Vx[i,j]  ) / Δx
        ∂Vy∂yW  = (Vy[i-1,j+1] - Vy[i-1,j]) / Δy
        ∂Vy∂yE  = (Vy[i,j+1]   - Vy[i,j]  ) / Δy
        ∂Vx∂yS  = (Vx[i,j] - VxS) / Δy
        ∂Vx∂yN  = (VxN - Vx[i,j]) / Δy
        ∂Vy∂xS  = (Vy[i,j]   - Vy[i-1,j]) / Δx
        ∂Vy∂xN  = (Vy[i,j+1] - Vy[i-1,j+1]) / Δx
        ∇vW     = ∂Vx∂xW + ∂Vy∂yW
        ∇vE     = ∂Vx∂xE + ∂Vy∂yE
        εxxW    = ∂Vx∂xW - 1//3*∇vW
        εxxE    = ∂Vx∂xE - 1//3*∇vE
        εyxS    = 0.5*(∂Vx∂yS + ∂Vy∂xS)
        εyxN    = 0.5*(∂Vx∂yN + ∂Vy∂xN)
        τxxW    = 2.0*ηc[i-1,j]*εxxW + γ*∇vW
        τxxE    = 2.0*ηc[i  ,j]*εxxE + γ*∇vE
        τyxS    = 2.0*ηv[i,j]*εyxS
        τyxN    = 2.0*ηv[i,j+1]*εyxN
        Fx[i,j] = (τxxE - τxxW)/Δx + (τyxN - τyxS)/Δy + bx[i,j]
    end

    for j=2:nvy-1, i=1:ncx # avoid Dirichlets
        if i==1
            VyW = bcv*2*BC.Vy.Dir_W[j] - Vy[i,j]
        else
            VyW = Vy[i-1,j]
        end
        if i==ncx  
            VyE = bcv*2*BC.Vy.Dir_E[j] - Vy[i,j] 
        else
            VyE = Vy[i+1,j] 
        end
        ∂Vx∂xS  = (Vx[i+1,j-1] - Vx[i,j-1]) / Δx
        ∂Vx∂xN  = (Vx[i+1,j-0] - Vx[i,j-0]) / Δx
        ∂Vy∂yS  = (Vy[i,j]   - Vy[i,j-1]) / Δy
        ∂Vy∂yN  = (Vy[i,j+1] - Vy[i,j]  ) / Δy
        ∂Vx∂yW  = (Vx[i,j]   - Vx[i,j-1]  ) / Δy
        ∂Vx∂yE  = (Vx[i+1,j] - Vx[i+1,j-1]) / Δy
        ∂Vy∂xW  = (Vy[i,j]   - VyW) / Δx
        ∂Vy∂xE  = (VyE - Vy[i,j]  ) / Δx
        ∇vS     = ∂Vx∂xS + ∂Vy∂yS
        ∇vN     = ∂Vx∂xN + ∂Vy∂yN
        εyyS    = ∂Vy∂yS - 1//3*∇vS
        εyyN    = ∂Vy∂yN - 1//3*∇vN
        εxyW    = 0.5*(∂Vy∂xW + ∂Vx∂yW) 
        εxyE    = 0.5*(∂Vy∂xE + ∂Vx∂yE)
        τyyS    = 2.0*ηc[i,j-1]*εyyS + γ*∇vS
        τyyN    = 2.0*ηc[i,j]*εyyN + γ*∇vN
        τxyW    = 2.0*ηv[i,j]*εxyW
        τxyE    = 2.0*ηv[i+1,j]*εxyE
        Fy[i,j] = (τyyN - τyyS)/Δy + (τxyE - τxyW)/Δx + by[i,j]
    end
end
##############
function Residual_Stokes!( Fx, Fy, Fp, Vx, Vy, P, p, γ, bcv, BC )
    Δx=p.Δx; Δy=p.Δy; BCx=p.BCx; BCy=p.BCy; ηv=p.ηv; ηc=p.ηc; ncx=p.ncx; ncy=p.ncy; bx=p.bx; by=p.by; bp=p.bp; nvx=ncx+1; nvy=ncy+1
    
    for j=1:ncy, i=2:nvx-1
        if j==1
            VxS = 2*BC.Vx.Dir_S[i] - Vx[i,j]
        else
            VxS = Vx[i,j-1]
        end
        if j==ncy  
            VxN = 2*BC.Vx.Dir_N[i] - Vx[i,j] 
        else
            VxN = Vx[i,j+1] 
        end
        ∂Vx∂xW  = (Vx[i,j]   - Vx[i-1,j]) / Δx
        ∂Vx∂xE  = (Vx[i+1,j] - Vx[i,j]  ) / Δx
        ∂Vy∂yW  = (Vy[i-1,j+1] - Vy[i-1,j]) / Δy
        ∂Vy∂yE  = (Vy[i,j+1]   - Vy[i,j]  ) / Δy
        ∂Vx∂yS  = (Vx[i,j] - VxS) / Δy
        ∂Vx∂yN  = (VxN - Vx[i,j]) / Δy
        ∂Vy∂xS  = (Vy[i,j]   - Vy[i-1,j]) / Δx
        ∂Vy∂xN  = (Vy[i,j+1] - Vy[i-1,j+1]) / Δx
        ∇vW     = ∂Vx∂xW + ∂Vy∂yW
        ∇vE     = ∂Vx∂xE + ∂Vy∂yE
        εxxW    = ∂Vx∂xW - 1//3*∇vW
        εxxE    = ∂Vx∂xE - 1//3*∇vE
        εyxS    = 0.5*(∂Vx∂yS + ∂Vy∂xS)
        εyxN    = 0.5*(∂Vx∂yN + ∂Vy∂xN)
        τxxW    = 2.0*ηc[i-1,j]*εxxW 
        τxxE    = 2.0*ηc[i  ,j]*εxxE 
        τyxS    = 2.0*ηv[i,j]*εyxS
        τyxN    = 2.0*ηv[i,j+1]*εyxN
        Fx[i,j] = (τxxE - τxxW)/Δx + (τyxN - τyxS)/Δy - (P[i,j] - P[i-1,j])/Δx
    end

    for j=2:nvy-1, i=1:ncx # avoid Dirichlets
        if i==1
            VyW = 2*BC.Vy.Dir_W[j] - Vy[i,j]
        else
            VyW = Vy[i-1,j]
        end
        if i==ncx  
            VyE = 2*BC.Vy.Dir_E[j] - Vy[i,j] 
        else
            VyE = Vy[i+1,j] 
        end
        ∂Vx∂xS  = (Vx[i+1,j-1] - Vx[i,j-1]) / Δx
        ∂Vx∂xN  = (Vx[i+1,j-0] - Vx[i,j-0]) / Δx
        ∂Vy∂yS  = (Vy[i,j]   - Vy[i,j-1]) / Δy
        ∂Vy∂yN  = (Vy[i,j+1] - Vy[i,j]  ) / Δy
        ∂Vx∂yW  = (Vx[i,j]   - Vx[i,j-1]  ) / Δy
        ∂Vx∂yE  = (Vx[i+1,j] - Vx[i+1,j-1]) / Δy
        ∂Vy∂xW  = (Vy[i,j]   - VyW) / Δx
        ∂Vy∂xE  = (VyE - Vy[i,j]  ) / Δx
        ∇vS     = ∂Vx∂xS + ∂Vy∂yS
        ∇vN     = ∂Vx∂xN + ∂Vy∂yN
        εyyS    = ∂Vy∂yS - 1//3*∇vS
        εyyN    = ∂Vy∂yN - 1//3*∇vN
        εxyW    = 0.5*(∂Vy∂xW + ∂Vx∂yW) 
        εxyE    = 0.5*(∂Vy∂xE + ∂Vx∂yE)
        τyyS    = 2.0*ηc[i,j-1]*εyyS 
        τyyN    = 2.0*ηc[i,j]*εyyN 
        τxyW    = 2.0*ηv[i,j]*εxyW
        τxyE    = 2.0*ηv[i+1,j]*εxyE
        Fy[i,j] = (τyyN - τyyS)/Δy + (τxyE - τxyW)/Δx - (P[i,j] - P[i,j-1])/Δy 
    end
end
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
function main( N )
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
    params.solver    = 1     # 1: coupled with pinned P, 0: coupled with slight compressibility --- 1: decoupled Powell-Hestenes --- 2: decoupled KSP GCR
    #------------------------ PINNED PRESSURE
    params.gamma     = 1e4   # Penalty factor
    params.Dir_scale = 1.0   # Dirichlet scaling factor
    # Non-linear iterations 
    params.niter_nl  = 1     # max. number of non-linear iterations
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

    Fx = zeros(ncx+1,ncy)
    Fy = zeros(ncx,ncy+1)
    Fp = zeros(ncx,ncy)
    Vx = zeros(ncx+1,ncy)
    Vy = zeros(ncx,ncy+1)
    bx = zeros(ncx+1,ncy)
    by = zeros(ncx,ncy+1)

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
            #-------------------------------------------------
            # @time StokesSolver!( fields, params, domain, materials, BC, Kuu, Kup, Kpu )
            

            fu = zeros(Float64, length(fields.Fx) + length(fields.Fy))
            fp = zeros(Float64, length(fields.Fp))
            fu[1:length(fields.Fx)]     .= fields.Fx[:]
            fu[length(fields.Fx)+1:end] .= fields.Fy[:]
            fp .= fields.Fp[:]

            # Decoupled solve
            ndofu = size(Kup,1)
            ndofp = size(Kup,2)
            if params.comp==1
                coef      = 1.0./(f.Kc[:].*params.dt).*ones(length(f.etac))#.*etac[:]
                coef_inv  = f.Kc[:].*params.dt.*ones(length(f.etac))#.*etac[:]
            else
                coef      =    0.0*ones(length(fields.etac))#.*etac[:]
                coef_inv  = params.gamma*ones(length(fields.etac))#.*etac[:]
            end
            Kpp   = spdiagm(coef)
            Kppi  = spdiagm(coef_inv)
            Kuusc = Kuu - Kup*(Kppi*Kpu) # OK
            
            println()
            show(stdout, "text/plain", Kuu[45,:])
            println() 
            show(stdout, "text/plain", Kuusc[45,:])
            # show(stdout, "text/plain", Kuusc[44,:])
            # show(stdout, "text/plain", Kuusc[45,:])
            
            PC    =  0.5*(Kuusc + Kuusc') 
            t = @elapsed Kf    = cholesky(Hermitian(PC), check = false)
            @printf("Cholesky took = %02.2e s\n", t)
            u     = zeros(Float64,ndofu)
            ru    = zeros(Float64,ndofu)
            fusc  = zeros(Float64,ndofu)
            p     = zeros(Float64,ndofp)
            rp    = zeros(Float64,ndofp)
            # Iterations
            for rit=1:10
                ru   .= fu .- Kuu*u .- Kup*p;
                rp   .= fp .- Kpu*u .- Kpp*p;
                @printf("  --> Powell-Hestenes Iteration %02d\n  Momentum res.   = %2.2e\n  Continuity res. = %2.2e\n", rit, norm(ru)/sqrt(length(ru)), norm(rp)/sqrt(length(rp)))
                if norm(ru)/(length(ru)) < 1e-13 && norm(rp)/(length(ru)) < 1e-13
                    break
                end
                fusc .=  fu .- Kup*(Kppi*fp .+ p)
                u    .= Kf\fusc
                p   .+= Kppi*(fp .- Kpu*u .- Kpp*p)

                # Residual
                # f = fusc .- Kuusc*u
                # Fx .= f[fields.NumVx]
                # Fy .= f[fields.NumVy]

                Rx = fu[fields.NumVx]
                Ry = fu[fields.NumVy]
                Rp = fp[fields.NumP]

                # Matrix free Schur-complement RHS
                p_dum = p[fields.NumP] .+ params.gamma.*Rp
                bx[2:end-1,:] .= Rx[2:end-1,:] .- diff(p_dum,dims=1) ./ dx  
                by[:,2:end-1] .= Ry[:,2:end-1] .- diff(p_dum,dims=2) ./ dy 

                # Matrix free Schur-complement residual: 🎫
                # bx = fusc[fields.NumVx]
                # by = fusc[fields.NumVy]
                Vx.= u[fields.NumVx]
                Vy.= u[fields.NumVy]
                bp = zeros(ncx,ncy)
                pm = (Δx=dx, Δy=dy, BCx=0, BCy=0, ηv=fields.etav, ηc=fields.etac, ncx=ncx, ncy=ncy, bx=bx, by=by, bp=bp )
                Residual_Stokes_SchurComplement!( Fx, Fy, Vx, Vy, pm, params.gamma, 0.0, BC )

                dP = params.gamma.*(Rp .- (diff(Vx,dims=1) ./ dx + diff(Vy,dims=2) ./ dy))
                p .+= dP[:]
            end
            fields.Vx[:,2:end-1] .= fields.Vx[:,2:end-1] .+ u[fields.NumVx]
            fields.Vy[2:end-1,:] .= fields.Vy[2:end-1,:] .+ u[fields.NumVy]
            fields.Pc            .= fields.Pc            .+ p[fields.NumP]
            
            #-------------------------------------------------
        
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
    println("L1 error P :", mean(errP))
    # Visualize
    if params.show_figs==1

        # display(fields.NumVx')

        # bx = zeros(ncx+1,ncy)
        # by = zeros(ncx,ncy+1)
        # bp = zeros(ncx,ncy)
        # pm = (Δx=dx, Δy=dy, BCx=0, BCy=0, ηv=fields.etav, ηc=fields.etac, ncx=ncx, ncy=ncy, bx=bx, by=by, bp=bp )
        # Residual_Stokes!( Fx, Fy, Fp, fields.Vx[:,2:end-1], fields.Vy[2:end-1,:], fields.Pc, pm, params.gamma, 1.0, BC )


        p1 = Plots.heatmap(domain.xv, domain.yc, Array(Fx)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Fx")
        p2 = Plots.heatmap(domain.xc, domain.yv, Array(Fy)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Fy")
        # p1 = Plots.heatmap(domain.xc, domain.yc, Array(Vxc .- Vxca)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Vx")
        # p2 = Plots.heatmap(domain.xc, domain.yc, Array(Vyc .- Vyca)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="Vy")
        p3 = Plots.heatmap(domain.xc, domain.yc, Array(errP)',        aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="P")
        p4 = Plots.heatmap(domain.xv, domain.yv, Array(fields.etav)', aspect_ratio=1, xlims=(domain.xmin, domain.xmax), ylims=(domain.ymin, domain.ymax), c=Plots.cgrad(:roma, rev = true), title="etav")
        display(Plots.plot( p1, p2, p3, p4, dpi=200 ) ) 
    end
    return mean(errV), mean(errP)
end

####################################################
####################################################
####################################################
N        = [20;]
# Call solver
for ires in eachindex(N)
    @time eV, eP = main( N[ires] )
end
