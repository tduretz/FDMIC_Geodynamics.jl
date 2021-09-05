##############

Base.@kwdef mutable struct Markers
    """ 
        Data structure for marker fields 
    """
    x         ::  Union{Array{Float64,1}, Missing} = missing
    y         ::  Union{Array{Float64,1}, Missing} = missing
    T         ::  Union{Array{Float64,1}, Missing} = missing
    phase     ::  Union{Array{Float64,1}, Missing} = missing
    cellx     ::  Union{Array{Int64,1}, Missing} = missing
    celly     ::  Union{Array{Int64,1}, Missing} = missing
    nmark     ::  Union{Int64, Missing} = missing
    nmark_ini ::  Union{Int64, Missing} = missing
    nmark_max ::  Union{Int64, Missing} = missing
    nmx       ::  Union{Int64, Missing} = missing
    nmy       ::  Union{Int64, Missing} = missing
end
export Markers

##############

Base.@kwdef mutable struct Thermal_BC
    type_W  ::Union{Int64, Missing} = 1 # Default Dirichlet
    type_E  ::Union{Int64, Missing} = 1 # Default Dirichlet
    type_S  ::Union{Int64, Missing} = 1 # Default Dirichlet
    type_N  ::Union{Int64, Missing} = 1 # Default Dirichlet
    Dir_W   ::Union{Vector{Float64}, Missing} = missing
    Dir_E   ::Union{Vector{Float64}, Missing} = missing
    Dir_S   ::Union{Vector{Float64}, Missing} = missing
    Dir_N   ::Union{Vector{Float64}, Missing} = missing
    Pinned  ::Union{Matrix{Int64}, Missing}   = missing
end
export Thermal_BC

##############

Base.@kwdef mutable struct BoundaryConditions
    T   :: Union{Thermal_BC, Missing} = Thermal_BC()
    Vx  :: Union{Thermal_BC, Missing} = Thermal_BC()
    Vy  :: Union{Thermal_BC, Missing} = Thermal_BC()
    periodix  ::Union{Int64, Missing} = 0
end
export BoundaryConditions

##############

Base.@kwdef mutable struct Tensor2D 
    xx   :: Union{Matrix{Float64}, Missing} = missing
    yy   :: Union{Matrix{Float64}, Missing} = missing
    zz   :: Union{Matrix{Float64}, Missing} = missing
    xy   :: Union{Matrix{Float64}, Missing} = missing
    xy_c :: Union{Matrix{Float64}, Missing} = missing
    II   :: Union{Matrix{Float64}, Missing} = missing
    div  :: Union{Matrix{Float64}, Missing} = missing
end
export Tensor2D

##############

Base.@kwdef mutable struct ModelDomain 
    xc   :: Union{LinRange{Float64}, Missing} = missing
    yc   :: Union{LinRange{Float64}, Missing} = missing 
    xv   :: Union{LinRange{Float64}, Missing} = missing 
    yv   :: Union{LinRange{Float64}, Missing} = missing 
    xce  :: Union{LinRange{Float64}, Missing} = missing 
    yce  :: Union{LinRange{Float64}, Missing} = missing 
    xmin :: Union{Float64, Missing} = missing
    xmax :: Union{Float64, Missing} = missing
    ymin :: Union{Float64, Missing} = missing
    ymax :: Union{Float64, Missing} = missing
    dx   :: Union{Float64, Missing} = missing
    dy   :: Union{Float64, Missing} = missing
    ncx  :: Union{Int64, Missing}   = missing
    ncy  :: Union{Int64, Missing}   = missing
    L    :: Union{Float64, Missing} = missing
    L0   :: Union{Float64, Missing} = missing
end
export ModelDomain

##############

Base.@kwdef mutable struct MaterialParameters 
    nphase       :: Union{Int64, Missing}   = missing
    # Viscosity
    n    :: Union{Vector{Float64}, Missing} = missing
    eta0 :: Union{Vector{Float64}, Missing} = missing 
    Tref :: Union{Vector{Float64}, Missing} = missing 
    # Elasticity
    K    :: Union{Vector{Float64}, Missing} = missing
    G    :: Union{Vector{Float64}, Missing} = missing
    # Phase field damage
    # lDam :: Union{Vector{Float64}, Missing} = missing
    # gDam :: Union{Vector{Float64}, Missing} = missing
    # eDam :: Union{Vector{Float64}, Missing} = missing
    lDam :: Union{Float64, Missing} = missing
    gDam :: Union{Float64, Missing} = missing
    eDam :: Union{Float64, Missing} = missing
end
export MaterialParameters

##############

Base.@kwdef mutable struct ModelParameters 
    # User interface
    user          = zeros(Float64, 10) 
    # Switches
    BC_type       :: Union{Int64, Missing}     = missing    # 1: Pure shear / 2: Simple shear / 3: Simple shear periodic
    comp          :: Union{Int64, Missing}     = missing    # Compressible
    # Time integration
    advection     :: Union{Int64, Missing}     = missing
    ALE           :: Union{Int64, Missing}     = missing    # Deform box 
    Ebg           :: Union{Float64, Missing}   = missing 
    nt            :: Union{Int64, Missing}     = missing
    dt            :: Union{Float64, Missing}   = missing
    Courant       :: Union{Float64, Missing}   = missing
    dt_var        :: Union{Int64, Missing}     = missing
    t             :: Union{Float64, Missing}   = missing
    rkw           = 1.0/6.0*[1.0 2.0 2.0 1.0] # RK4 weights for averaging
    rkv           = 1.0/2.0*[1.0 1.0 2.0 2.0] # RK4 weights for time stepping
    # Linear solver
    solver        :: Union{Int64, Missing}     = missing    # 0: coupled --- 1: decoupled 
    gamma         :: Union{Float64, Missing}   = missing    # Penalty factor
    Dir_scale     :: Union{Float64, Missing}   = missing    # Dirichlet scaling factor
    # Non-linear solver 
    niter_nl      :: Union{Int64, Missing}     = missing    # max. number of non-linear iterations
    tol_nl        :: Union{Float64, Missing}   = missing    # non-linear tolerance
    JFNK          :: Union{Int64, Missing}     = missing    # Jacobian-Free Newton_Krylov
    # Visualisation
    show_figs     :: Union{Int64, Missing}     = missing    # activates visualisation...
    nout          :: Union{Int64, Missing}     = missing    # ... every nout
end
export ModelParameters

##############

Base.@kwdef mutable struct Fields2D 
    etac      :: Union{Matrix{Float64}, Missing} = missing
    etav      :: Union{Matrix{Float64}, Missing} = missing
    Kc        :: Union{Matrix{Float64}, Missing} = missing
    Gc        :: Union{Matrix{Float64}, Missing} = missing
    Gv        :: Union{Matrix{Float64}, Missing} = missing
    We        :: Union{Matrix{Float64}, Missing} = missing
    We0       :: Union{Matrix{Float64}, Missing} = missing
    rhoc      :: Union{Matrix{Float64}, Missing} = missing
    Pc0       :: Union{Matrix{Float64}, Missing} = missing
    Pc        :: Union{Matrix{Float64}, Missing} = missing
    Vx        :: Union{Matrix{Float64}, Missing} = missing
    Vy        :: Union{Matrix{Float64}, Missing} = missing
    Vxc       :: Union{Matrix{Float64}, Missing} = missing
    Vyc       :: Union{Matrix{Float64}, Missing} = missing
    div       :: Union{Matrix{Float64}, Missing} = missing
    Fx        :: Union{Matrix{Float64}, Missing} = missing
    Fy        :: Union{Matrix{Float64}, Missing} = missing
    Fp        :: Union{Matrix{Float64}, Missing} = missing
    FDam      :: Union{Matrix{Float64}, Missing} = missing
    Damc      :: Union{Matrix{Float64}, Missing} = missing
    Damv      :: Union{Matrix{Float64}, Missing} = missing
    phiDam    :: Union{Matrix{Float64}, Missing} = missing
    phiDam0   :: Union{Matrix{Float64}, Missing} = missing
    GDam      :: Union{Matrix{Float64}, Missing} = missing
    lDam      :: Union{Matrix{Float64}, Missing} = missing
    #####
    NumVx     :: Union{Matrix{Int64}, Missing} = missing
    NumVy     :: Union{Matrix{Int64}, Missing} = missing
    NumP      :: Union{Matrix{Int64}, Missing} = missing
    #####
    mpc           :: Union{Matrix{Float64}, Missing} = missing
    mpc_th        :: Union{Array{Float64, 3}, Missing} = missing
    phase_perc    :: Union{Array{Float64, 3}, Missing} = missing
    phase_perc_th :: Union{Array{Float64, 4}, Missing} = missing
end
export Fields2D

##############

Base.@kwdef mutable struct CharScales 
    L      :: Union{Float64, Missing}   = missing
    t      :: Union{Float64, Missing}   = missing
    S      :: Union{Float64, Missing}   = missing
    T      :: Union{Float64, Missing}   = missing
end
export CharScales
