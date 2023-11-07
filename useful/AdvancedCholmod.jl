using SparseArrays
using LinearAlgebra
using SuiteSparse
using BenchmarkTools
using UnicodePlots

import SuiteSparse.CHOLMOD: Sparse, Dense, Factor, C_Dense, C_Factor, C_Sparse, SuiteSparse_long
import SuiteSparse.CHOLMOD: common_struct, common_supernodal, common_nmethods, change_factor!, set_print_level, defaults
import SuiteSparse.CHOLMOD: VTypes, @cholmod_name, fact_, cholesky!, CHOLMOD_A
import SparseArrays: getcolptr, getrowval

function supercholesky(A::Sparse; shift::Real=0.0, check::Bool = true,
  perm::Union{Nothing,AbstractVector{SuiteSparse_long}}=nothing)
  
  cm = defaults(common_struct[Threads.threadid()])
  set_print_level(cm, 0)
  
  # Force a supernodal solution (eliminates alloc on solve)
  unsafe_store!(common_supernodal[Threads.threadid()], 2)
  
  # Compute the symbolic factorization
  @time F = fact_(A, cm; perm = perm)
  
  # Compute the numerical factorization
  @time cholesky!(F, A; shift = shift, check = check)
  
  return F
  
end

function _solve!(X::Dense{Tv}, Yref::Ref{Ptr{C_Dense{Tv}}}, Eref::Ref{Ptr{C_Dense{Tv}}}, F::Factor{Tv}, B::Dense{Tv}) where Tv<:VTypes
  
  # Pointer to pre-allocated dense matrix
  Xref = Ptr{C_Dense{Tv}}(pointer(X))
  sys = CHOLMOD_A # Solve type
  Bset = C_NULL   # Unused parameter
  Xset = C_NULL   # Unused parameter
  
  if size(F,1) != size(B,1)
    throw(DimensionMismatch("LHS and RHS should have the same number of rows. " *
    "LHS has $(size(F,1)) rows, but RHS has $(size(B,1)) rows."))
  end
  
  if !issuccess(F)
    s = unsafe_load(pointer(F))
    if s.is_ll == 1
      throw(LinearAlgebra.PosDefException(s.minor))
    else
      throw(LinearAlgebra.ZeroPivotException(s.minor))
    end
  end
  
  res = ccall((@cholmod_name("solve2"), :libcholmod), Cint,
  # Input                                                        # Output                                         # Workspace
  (Cint, Ptr{C_Factor{Tv}}, Ptr{C_Dense{Tv}}, Ptr{C_Sparse{Tv}}, Ref{Ptr{C_Dense{Tv}}},  Ref{Ptr{C_Sparse{Tv}}},  Ref{Ptr{C_Dense{Tv}}},  Ref{Ptr{C_Dense{Tv}}}, Ptr{UInt8}),
  sys,   F,                   B,              Bset,              Xref,                   Xset,                    Yref,                   Eref,                  SuiteSparse.CHOLMOD.common_struct[Threads.threadid()])
  
  if(res != 1)
    throw(ErrorException("CHOLMOD solve failure"))
  end

  return 
  
end


T = Float64
N = 8000

X = Dense(Array{T}(undef, N))
B = Dense(rand(N))
Yref = Ref(Ptr{C_Dense{T}}(C_NULL))
Eref = Ref(Ptr{C_Dense{T}}(C_NULL))

A = sprand(N, N, 0.005)
A = A + A' + N*I
As = Sparse(A)

F_simple = cholesky(A)
F_super =  supercholesky(As);

_solve!(X, Yref, Eref, F_simple, B);
@show F_simple
A\B ≈ X ? println("Simplical pass") : error("Simplical fail")

_solve!(X, Yref, Eref, F_super, B);
@show F_super
A\B ≈ X ? println("Supernodal pass") : error("Supernodal fail")

# @benchmark _solve!($X, $Yref, $Eref, $F_simple, $B)

# @benchmark _solve!($X, $Yref, $Eref, $F_super, $B)

# @benchmark $F_super\B

# Factorization + solve
cm = defaults(common_struct[Threads.threadid()])
@time Asymb = fact_(As, cm)
@time cholesky!(Asymb,As)


println("solves...")

@time _solve!(X, Yref, Eref, F_simple, B);

@time _solve!(X, Yref, Eref, F_super, B);

@time F_super\B;