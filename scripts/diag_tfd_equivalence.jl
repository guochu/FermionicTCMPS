# 正向对照 + 反向差异一次性呈现（K=1 有限温热杂质，两侧 IMB 二次型）
using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase
using ImpurityModelBase: discretefermionicbath, fermidirac, frequencies, spectrumcouplings,
                         quadratichamiltonian, tunneling, adaga, adagadag, aa
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))

ws = [0.4]; fs = [0.55]
x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
V = spectrumcouplings(x.bath.f)[1]
_, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
εd = x.εd; ω = ws[1]; nF = fermidirac(0.8, 0.0, ω)

# star: IMB normal quadratic + thermal ρ₀ (imp-hot ⊗ diag(n_F))
ham_star = quadratichamiltonian(2, [tunneling(1, 1; coeff = εd), tunneling(2, 2; coeff = ω),
                              tunneling(1, 2; coeff = x.J * V), tunneling(2, 1; coeff = x.J * V)])
h_star = ImpurityModelBase.cmatrix(ham_star)
ρ0s = diagm([1.0, nF])
n_star(t) = real((exp(-im * h_star * t) * ρ0s * exp(im * h_star * t))[1, 1])
# star (t,0)-correlator C(t) = ⟨d†(0)d(t)⟩ = [e^{−iht}ρ₀]_{11}
C_star(t) = (exp(-im * h_star * t) * ρ0s)[1, 1]

# chain: IMB generic quadratic (BdG) — εd n_d + ω n_{a₁} − ω n_{f₂} + v1(d†a₁+h.c.) + v2(d†f₂†+f₂d)
ham_chain = quadratichamiltonian(3, [tunneling(1, 1; coeff = εd), tunneling(2, 2; coeff = ω),
                              tunneling(3, 3; coeff = -ω), tunneling(1, 2; coeff = v1[1]),
                              tunneling(2, 1; coeff = v1[1]), adagadag(1, 3; coeff = v2[1]),
                              aa(3, 1; coeff = v2[1])])
# chain exact evolution: 4^3 = 64-dim ED on the G2-pinned dense H
spec = FermionicTCMPS.chain_spec(x); L = length(spec.kinds)
couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)]
H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
F = eigen(Hermitian(H))
flat(occs) = 1 + sum((occs[s][1] + 2occs[s][2]) * 4^(L - s) for s in 1:L)
occs0 = ntuple(i -> i == spec.imp_site ? (1, 0) : (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
ψ0 = zeros(ComplexF64, 4^L); ψ0[flat(occs0)] = 1.0
g0 = F.vectors' * ψ0
nddiag = [(i - 1) ÷ 4^(L - spec.imp_site) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]
println("t     star_n(t)   chain_n(t)")
for t in (0.1, 0.2, 0.3, 0.4)
    ψt = F.vectors * (exp.(-im .* F.values .* t) .* g0)
    n_ch = real(ψt' * Diagonal(nddiag) * ψt)
    println(round(t; digits=3), "   ", round(n_star(t); digits=6), "   ", round(n_ch; digits=6))
end
