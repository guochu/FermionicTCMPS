using Logging: @warn
using Base.Threads: @spawn
using LinearAlgebra: Hermitian, eigen, Diagonal, mul!, rmul!, axpy!, diagm, tr
using Random, QuadGK, Reexport, TensorKit
@reexport using DMRG, GeneralHamiltonians, TEBD, ImpurityModelBase
using GeneralHamiltonians: AbstractCoefficient


# auxiliary
include("auxiliary/freegf.jl")
include("auxiliary/rungekutta.jl")
include("auxiliary/discretization.jl")

# fermionic hamiltonian wrapper
include("operators/fterm.jl")
include("operators/abstractdefs.jl")
include("operators/hamiltonian.jl")
include("operators/freehamiltonian.jl")


# TCMPS algorithm

# utility
include("tcmps/star_to_chain.jl")
include("tcmps/quenchscheme.jl")
include("tcmps/transformations.jl")
include("tcmps/abstractmodels.jl")
include("tcmps/current.jl")
include("tcmps/separablestate.jl")

# # TCMPS algorithm
# single impurity
include("tcmps/onebath/abstractmodel.jl")
include("tcmps/onebath/oneimpurity.jl")
include("tcmps/onebath/thermalstate.jl")
include("tcmps/onebath/freefermions.jl")
# two impurity
include("tcmps/twobath/abstractmodel.jl")
include("tcmps/twobath/twoimpurity.jl")
include("tcmps/twobath/thermalstate.jl")
include("tcmps/twobath/freefermions.jl")
# green's functions
include("tcmps/gf/gf.jl")