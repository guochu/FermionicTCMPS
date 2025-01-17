push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/DMRG/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/InfiniteDMRG/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/GeneralHamiltonians/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/TEBD/src")

using Test, Random
using TensorKit, DMRG, GeneralHamiltonians, TEBD

push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")
using FermionicTCMPS

# include("../src/includes.jl")

Random.seed!(12354)

include("util.jl")

### tcmps
include("tcmps/configurations.jl")
# ground state and thermal state quality for single impurity model
include("tcmps/eq_onebath.jl")
include("tcmps/neq_onebath.jl")

# thermal state quality check for double impurity model
include("tcmps/neq_twobath.jl")

include("tcmps/eq_twobath.jl")
