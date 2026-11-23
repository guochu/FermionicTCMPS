# Apply-path vs expectation-path: does MPO * MPS lose the hop paths?
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
basis = [_basis_state(spec, _occ_of(L, w)) for w in 1:4^L]
# pick a bath-filled basis state (config with a₂ filled): find one
idx_fill = [w for w in 1:4^L if any(o -> o == (1, 1), _occ_of(L, w))][1]
println("source basis idx = ", idx_fill, " occs = ", _occ_of(L, idx_fill))
applied = Wmpo * basis[idx_fill]
for v in 1:4^L
    dv = Wdense[v, idx_fill]
    abs(dv) > 1e-8 || continue
    ov = dot(basis[v], applied)
    println("  ⟨v=$v|apply|w⟩ = ", round(ov, digits = 8), "   dense = ", round(dv, digits = 8))
end
