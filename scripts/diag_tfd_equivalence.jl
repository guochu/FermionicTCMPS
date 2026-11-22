# TFD-equivalence protocol (K=1 finite-β hot impurity, star vs chain, both dense).
#
# History: under the previous chain (KS Eq. (9) T=0 structure carried to finite T)
# this protocol EXPOSED the reported bug — chain n(t) fell below the star frame
# from O(t³) (0.898 vs 0.941 at t = 0.4). The chain now carries the EXACT
# Bogoliubov image (weighted on-sites + intra-mode pairing + thermofield-double
# fiducial), so both equal-time occupations and the (0,t) GFs reproduce the star
# frame at machine precision.
#
#   julia --project=. scripts/diag_tfd_equivalence.jl
using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase
using ImpurityModelBase: discretefermionicbath, fermidirac, frequencies, spectrumcouplings,
                         quadratichamiltonian, tunneling, freefermions_greater_lesser
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))

ws = [0.4]; fs = [0.55]
x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
V = spectrumcouplings(x.bath.f)[1]
_, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
εd = x.εd; ω = ws[1]; n = ns[1]
println("V = ", V, "  v1 = ", x.J * v1[1], "  v2 = ", x.J * v2[1], "  n_F = ", n)

# star: IMB normal quadratic + thermal ρ₀ (imp-hot ⊗ diag(n_F))
ham_star = quadratichamiltonian(2, [tunneling(1, 1; coeff = εd), tunneling(2, 2; coeff = ω),
                              tunneling(1, 2; coeff = x.J * V), tunneling(2, 1; coeff = x.J * V)])
h_star = ImpurityModelBase.cmatrix(ham_star)
ρ0s = diagm([1.0, n])
n_star(t) = real((exp(-im * h_star * t) * ρ0s * exp(im * h_star * t))[1, 1])
gt_star(t) = -im * ((exp(-im * h_star * t) * (I - ρ0s))[1, 1])
lt_star(t) = im * ((ρ0s * exp(-im * h_star * t))[1, 1])

# chain: dense mechanical-JW H (6-tuple couples with the intra-mode pairing),
# TFD fiducial, exact eigendecomposition
spec = FermionicTCMPS.chain_spec(x); L = length(spec.kinds)
couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k], ns[k]) for k in eachindex(ws)]
H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
modes = [(c[1], c[2], c[6]) for c in couples]
ψ0 = dense_tfd_vector(spec, modes, (1, 0))

F = eigen(Hermitian(H))
g0 = F.vectors' * ψ0
nddiag = [(i - 1) ÷ 4^(L - spec.imp_site) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]

times = collect(0.0:0.1:0.6)
gt_ch, lt_ch = dense_gf_sandwich(H, ψ0, spec.imp_site, times; imp_up = 1)
hu = h_star
gtr_e, = freefermions_greater_lesser(hu, diagm(vcat([0.0], n)), times, 1, 1)
_, ltr_f = freefermions_greater_lesser(hu, diagm(vcat([1.0], n)), times, 1, 1)

println("\nt     star_n(t)   chain_n(t)   |Δn|")
global allΔ = 0.0
for (j, t) in enumerate(times)
    ψtt = F.vectors * (exp.(-im .* F.values .* t) .* g0)
    n_ch = real(sum(conj(ψtt) .* nddiag .* ψtt))
    Δ = abs(n_star(t) - n_ch)
    global allΔ = max(allΔ, Δ)
    println(round(t; digits=3), "   ", round(n_star(t); digits=6), "      ", round(n_ch; digits=6), "      ", round(Δ, digits=12))
end

println("\nGF layer (per-spin star resolvent vs spin-summed chain sandwich):")
# the hot impurity: only the EMPTY down spin contributes to the chain's
# greater, and only the OCCUPIED up spin to its lesser — compare per spin
maxgt = maximum(abs(gt_ch[j] - gtr_e[j]) for j in eachindex(times))
maxlt = maximum(abs(lt_ch[j] - ltr_f[j]) for j in eachindex(times))
println("   max |Δ greater| = ", maxgt, "   max |Δ lesser| = ", maxlt)
println("\nmax |Δ n(t)| = ", allΔ)
