##############

mutable struct Markers
    x         ::  Array{Float64,1}
    y         ::  Array{Float64,1}
    T         ::  Array{Float64,1}
    phase     ::  Array{Float64,1}
    cellx     ::  Array{Int64,1}#Vector{CartesianIndex{2}}
    celly     ::  Array{Int64,1}
    nmark     ::  Int64
    nmark_max ::  Int64
end

##############

Base.@kwdef mutable struct Thermal_BC
    type_W  ::Union{Int64, Missing} = 1 # Defautl Dirichlet
    type_E  ::Union{Int64, Missing} = 1 # Defautl Dirichlet
    type_S  ::Union{Int64, Missing} = 1 # Defautl Dirichlet
    type_N  ::Union{Int64, Missing} = 1 # Defautl Dirichlet
    Dir_W   ::Union{Vector{Float64}, Missing} = missing
    Dir_E   ::Union{Vector{Float64}, Missing} = missing
    Dir_S   ::Union{Vector{Float64}, Missing} = missing
    Dir_N   ::Union{Vector{Float64}, Missing} = missing
end

##############

Base.@kwdef mutable struct BoundaryConditions
    T   :: Union{Thermal_BC, Missing} = Thermal_BC()
    Vx  :: Union{Thermal_BC, Missing} = Thermal_BC()
    Vy  :: Union{Thermal_BC, Missing} = Thermal_BC()
    periodix  ::Union{Int64, Missing} = 0
end

##############

Base.@kwdef mutable struct Tensor2D 
    xx   :: Union{Matrix{Float64}, Missing} 
    yy   :: Union{Matrix{Float64}, Missing} 
    zz   :: Union{Matrix{Float64}, Missing} 
    xy   :: Union{Matrix{Float64}, Missing} 
    xy_c :: Union{Matrix{Float64}, Missing} 
    II   :: Union{Matrix{Float64}, Missing} 
end

##############