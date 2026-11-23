# Is MPO(W) == exp(-iH δw) as an operator? Full 64×64 comparison on K=1.
using Random, TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath
using LinearAlgebra
const F = FermionicTCMPS
include(joinpath(pwd(), "test", "anchor_util.jl"))
Random.seed!(31)
ws = [-0.9]; fs = [0.35]
β = 0.8; δw = 0.0025
x = SISBD(discretefermionicbath(ws, fs; β, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
spec = F.chain_spec(x)
L = length(spec.kinds)
imp = spec.imp_site
_, v1, v2, ns = F.thermofield_couplings(x.bath)
couples = [(2, 3, ws[1], x.J * v1[1], x.J * v2[1], ns[1])]
H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
Wdense = exp(-im * δw * H)

h = mpo(x)
Wmpo = DMRG.MPO(DMRG.timeevompo(h, -im * δw, DMRG.WI()))

# dense matrix elements of the W-MPO
basis = [_basis_state(spec, _occ_of(L, w)) for w in 1:4^L]
Wm = zeros(ComplexF64, 4^L, 4^L)
for w in 1:4^L
    for v in 1:4^L
        _qtotal(spec, _occ_of(L, v)) == _qtotal(spec, _occ_of(L, w)) || continue
        Wm[v, w] = DMRG.expectation(basis[v], Wmpo, basis[w])
    end
end
err = maximum(abs.(Wm .- Wdense))
println("max |MPO(W) − exp(−iHδw)| = ", round(err, digits = 10))
# where are the largest deviations?
dev = abs.(Wm .- Wdense)
for idx in sortperm(vec(dev), rev = true)[1:5]
    (dev[idx] < 1e-10) && break
    println("   [", (idx-1) ÷ 4^L + 1, ", ", (idx-1) % 4^L + 1, "] dev = ", round(dev[idx], digits = 8),
            "  Wm = ", round(Wm[idx], digits = 6), "  dense = ", round(Wdense[idx], digits = 6))
end
