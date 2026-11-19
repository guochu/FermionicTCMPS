using Test, Random
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum,
                         frequencies, spectrumvalues

Random.seed!(20261119)

include("g2_dense_anchor.jl")
