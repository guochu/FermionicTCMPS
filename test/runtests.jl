using Test, Random
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum,
                         frequencies, spectrumvalues, spectrumcouplings, fermidirac,
                         freefermions_greater_lesser

Random.seed!(20261119)

include("g2_dense_anchor.jl")
include("g3_dynamics.jl")
include("gtempo_crosscheck.jl")
