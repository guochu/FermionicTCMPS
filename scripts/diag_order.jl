# Order analysis of the double-bath sweep error: is the TDVP2 trajectory O(dt)?
# ⚠️ BROKEN probe — do not quote its numbers. Two independent defects:
#   1. the evolved chain `_shifted_state` is the sharp VACUUM product chain,
#      while the reference χv is the finite-β thermofield-double insertion —
#      different initial states, so the "error" is dominated by the mismatch
#      (measured: 0.19/0.30/0.44 at t=0.3/0.5/1.0, flat in δt);
#   2. `copy(χ0m)` shares the site tensors with the evolved MPS (sweep! replaces
#      .data entries in place), so the raw overlap needs a deepcopy snapshot of
#      the initial state.
# The corrected probe is scripts/diag_order_plain.jl (proper per-component
# thermofield insertion + deepcopy snapshot): err 2.7e-4/3.9e-4/4.0e-4 at
# t=0.3/0.5/1.0, δt-independent (projection floor, not Trotter).
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

# reference: raw overlap <chi_up|chi(t)> from the dense eigendecomposition
modes = vcat([(2k - 1, 2k, nsl[k]) for k in eachindex(wsl)],
             [(imp + 2k - 1, imp + 2k, nsr[k]) for k in eachindex(wsr)])
ψ0 = dense_tfd_vector(spec, modes, (0, 0))
χv = zeros(ComplexF64, 4^L)
for i in eachindex(ψ0)
    iszero(ψ0[i]) && continue
    d = (i - 1) ÷ 4^(L - imp) % 4
    q_in = (d == 1 || d == 3 ? 1 : 0, d >= 2 ? 1 : 0)
    (q_in[1] == 0) || continue
    amp = FermionicTCMPS._amplitude(FermionicTCMPS._adag_up, (1, q_in[2]), q_in)
    amp === nothing && continue
    χv[i + 4^(L - imp)] += amp * ψ0[i]
end
F = eigen(Hermitian(H))
g = F.vectors' * χv
SAMPLE_T = [0.3, 0.5, 1.0]
ovl_ref = [sum(abs2(g[k]) * exp(-im * F.values[k] * t) for k in eachindex(g)) for t in SAMPLE_T]

χ0m = FermionicTCMPS._shifted_state(x, 0, 0, :up, 1)

function run_sweep(δt; fresh_cache::Bool)
    ψ = copy(χ0m)
    cache = environments(h, ψ)
    alg = FermionicTCMPS._tdvp2(δt, TRUNC)
    nsteps = round(Int, 1.0 / δt)
    out = Dict{Float64,ComplexF64}()
    for n in 1:nsteps
        fresh_cache && (cache = environments(h, ψ))
        sweep!(cache, alg)
        t = n * δt
        if any(abs(t - st) < 1e-9 for st in SAMPLE_T)
            out[round(t; digits = 6)] = dot(χ0m, ψ)
        end
    end
    return out
end

println("δt      fresh?   err(0.3)      err(0.5)       err(1.0)")
for δt in (0.02, 0.01, 0.005)
    for fresh in (false, true)
        res = run_sweep(δt; fresh_cache = fresh)
        errs = [abs(res[round(t; digits = 6)] - ovl_ref[j]) for (j, t) in enumerate(SAMPLE_T)]
        println(δt, "  ", fresh ? "yes" : "no ", "   ",
                [round(e, digits = 10) for e in errs])
    end
end
