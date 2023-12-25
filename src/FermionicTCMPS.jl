module FermionicTCMPS

# auxiliary
export FermionicCommutator, RK2Stepper, RK4Stepper
export DiscretizationScheme, LinearDiscretization, frequencies, spectrum_couplings

# fermionic hamiltonian wrapper
export AbstractFermionicTerm, TwoBodyTerm, FourBodyTerm, consolidate, changepositions, AbstractFermionicHamiltonian, FermionicHamiltonian, hamiltonian
export FreeFermionicHamiltonian, free, coefficient_matrix


# TCMPS backend
# utilities
export QuenchScheme, KohnSantoroQuench, quench_time
export StarConfiguration, StarThermofieldConfiguration, ChainConfiguration, ChainThermofieldConfiguration
export star, thermofield, chainmapping, couplings


# single impurity model
export AbstractDIM, OneBathDIM, TwoBathDIM 
export default_sys_sites, default_env_sites
export hybridization_hamiltonian, quench_hamiltonian, environment_hamiltonian, sysbath_tunneling
export SISBD, separable_initial_state
export FreeSISBD, quench
export separable_state, thermal_state, ground_state


# double impurity model
export SIDBD, GIDBD, FreeSIDBD, FreeGIDBD, left_sysbath_tunneling, right_sysbath_tunneling
export IRLMD, FreeIRLMD, add_sys!

# green's functions
export gf_greater_t, gf_greater_τ, gf_lesser_t, gf_greater_lesser_t
export gs_gf_greater_t, gs_gf_greater_τ, gs_gf_lesser_t, gs_gf_greater_lesser_t


using Logging: @warn
using Base.Threads: @spawn
using LinearAlgebra: Hermitian, eigen, Diagonal, mul!, rmul!, axpy!, diagm, tr
using Random, QuadGK, Reexport
@reexport using DMRG, Hamiltonians, TEBD, ImpurityModelBase
using Hamiltonians: AbstractCoefficient


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


# utilities for real spectrum functions
include("utilities/utilities.jl")


end