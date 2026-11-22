# One-step error scaling of the engine sweep vs exact evolution:
#   err(dt) = |<chi|chi_sweep(dt)> - <chi|chi_exact(dt)>|
# second-order scheme => ratio 8 under dt halving; first-order => ratio 4.
using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))
Random.seed!(20261119)
xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
spec = FermionicTCMPS.chain_spec(x)
L = length(spec.kinds); imp = spec.imp_site
wsl, vl1, vl2, nsl = FermionicTCMPS.thermofield_couplings(x.bathl)
wsr, vr1, vr2, nsr = FermionicTCMPS.thermofield_couplings(x.bathr)
couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k]) for k in eachindex(wsl)],
               [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k]) for k in eachindex(wsr)])
H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = -0.3, U = 1.7)
h = mpo(x)
TRUNC = truncdimcutoff(D = 256, ϵ = 1e-10, add_back = 10)
modes = vcat([(2k - 1, 2k, nsl[k]) for k in eachindex(wsl)],
             [(imp + 2k - 1, imp + 2k, nsr[k]) for k in eachindex(wsr)])
ψ0v = dense_tfd_vector(spec, modes, (0, 0))
χv = zeros(ComplexF64, 4^L)
for i in eachindex(ψ0v)
    iszero(ψ0v[i]) && continue
    d = (i - 1) ÷ 4^(L - imp) % 4
    q_in = (d == 1 || d == 3 ? 1 : 0, d >= 2 ? 1 : 0)
    (q_in[1] == 0) || continue
    amp = FermionicTCMPS._amplitude(FermionicTCMPS._adag_up, (1, q_in[2]), q_in)
    amp === nothing && continue
    χv[i + 4^(L - imp)] += amp * ψ0v[i]
end
Feig = eigen(Hermitian(H))
gv = Feig.vectors' * χv
χ0m = FermionicTCMPS._shifted_state(x, 0, 0, :up, 1)
for δt in (0.02, 0.01, 0.005)
    exact = sum(abs2(gv[k]) * exp(-im * Feig.values[k] * δt) for k in eachindex(gv))
    cache = environments(h, deepcopy(χ0m))
    alg = FermionicTCMPS._tdvp2(δt, TRUNC)
    sweep!(cache, alg)
    swp = dot(χ0m, cache.mps)
    println("δt=", δt, ":  one-step err = ", abs(swp - exact),
            "   (sweep ", round(swp, digits = 9), "  exact ", round(exact, digits = 9), ")")
end
