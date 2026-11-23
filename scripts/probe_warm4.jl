# Decisive: does the W application populate the hop amplitude?
using Random, TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath
using LinearAlgebra
const F = FermionicTCMPS
Random.seed!(31)
ws = [-0.9, 0.5]; fs = [0.35, 0.55]
x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
h = mpo(x)
δw = 0.0025
Wm = DMRG.MPO(DMRG.timeevompo(h, -im * δw, DMRG.WI()))
ψ = fiducial_state(x)
c = ψ.comps[1]
occ = occupation_op(x; spin = :up)
env0 = DMRG.environments(h, DMRG.ExactMPS(c))
println("⟨n⟩ before = ", round(real(DMRG.expectation(env0.state, occ, env0.state)), digits = 12))
capp = Wm * c
println("after apply (uncompressed): bonds = ", [dim(space(capp.data[b], 3)) for b in 1:5])
env1 = DMRG.environments(h, DMRG.ExactMPS(capp))
println("⟨n⟩ after apply = ", round(real(DMRG.expectation(env1.state, occ, env1.state)), digits = 12))
canonicalize!(capp, alg = Orthogonalize(alg = DMRG.SVD(), trunc = truncdimcutoff(D = 256, ϵ = 1e-12, add_back = 0)))
println("after compress: bonds = ", [dim(space(capp.data[b], 3)) for b in 1:5])
env2 = DMRG.environments(h, DMRG.ExactMPS(capp))
println("⟨n⟩ after compress = ", round(real(DMRG.expectation(env2.state, occ, env2.state)), digits = 12))
