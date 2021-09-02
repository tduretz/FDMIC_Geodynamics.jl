module FDMIC_Geodynamics
using Printf, Statistics
using LoopVectorization
import Plots
using LinearAlgebra, SparseArrays 
import UnicodePlots
using Base.Threads
include("DataStructures.jl")
include("ThermalRoutines.jl")
include("MechanicsRoutines.jl")
include("MechanicsSolvers.jl")
include("MarkerRoutines.jl")
include("GridRoutines.jl")
include("SparseRoutines.jl")
include("DamageRoutines.jl")
end
