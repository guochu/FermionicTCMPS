# Plain-route raw-trajectory accuracy probe, CORRECTED insertion side. The
# older (removed) diag_order.jl / diag_onestep.jl probes evolved `_shifted_state`
# — the sharp VACUUM product chain — against a finite-β TFD insertion reference;
# the two initial states differ, so their "error" was dominated by the
# initial-state mismatch, not the sweep (and a shallow `copy` shares tensors
# with the in-place sweep, so the raw overlap needs a deepcopy snapshot).
# Here both sides use the thermofield-double insertion:
#
#   sweep side: χ = _inserted_state(x, :up, 1; ψ = fiducial_state(x)) (the exact
#   per-component mechanical d†-apply the GF pipeline itself uses), evolved with
#   the plain two-site TDVP; raw overlap ovl(t) = Σ_j w_j ⟨χ_j(0)|χ_j(t)⟩.
#
#   reference side: dense eigendecomposition of the mechanical-JW Hamiltonian
#   applied to the dense TFD insertion vector (test/anchor_util.jl assembly).
#
# Diagnostics: error vs t (shape), error vs δt (Trotter vs projection), and the
# fresh-cache control (environment staleness).
#
#   julia --project=. scripts/diag_order_plain.jl
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

# reference: raw overlap ⟨χ_up|χ(t)⟩ from the dense eigendecomposition of the
# thermofield-double insertion vector
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

println("reference overlaps: ", [round(v, digits = 6) for v in ovl_ref])

"plain-route insertion chain, evolved with the production sweep; per-component
raw overlap against the un-evolved chain"
function run_sweep(δt; fresh_cache::Bool)
    ψ = FermionicTCMPS.fiducial_state(x)
    χ = FermionicTCMPS._inserted_state(x, :up, 1; ψ)
    # snapshot BEFORE the caches attach: sweep! mutates χ.comps in place
    χ0 = [deepcopy(c) for c in χ.comps]
    cache = [environments(h, c) for c in χ.comps]
    alg = FermionicTCMPS._tdvp2(δt, TRUNC)
    nsteps = round(Int, 1.0 / δt)
    out = Dict{Float64,ComplexF64}()
    for n in 1:nsteps
        fresh_cache && (cache = [environments(h, c) for c in χ.comps])
        foreach(c -> sweep!(c, alg), cache)
        t = n * δt
        if any(abs(t - st) < 1e-9 for st in SAMPLE_T)
            v = sum(w * dot(c0, c.state) for (w, c0, c) in zip(χ.weights, χ0, cache))
            out[round(t; digits = 6)] = χ.phase * ψ.phase * v
        end
    end
    return out
end

println("δt      fresh?   err(0.3)      err(0.5)       err(1.0)      ovl(0.3)")
for δt in (0.02, 0.01, 0.005)
    for fresh in (false, true)
        res = run_sweep(δt; fresh_cache = fresh)
        errs = [abs(res[round(t; digits = 6)] - ovl_ref[j]) for (j, t) in enumerate(SAMPLE_T)]
        println(δt, "  ", fresh ? "yes" : "no ", "   ",
                [round(e, digits = 10) for e in errs], "   ", round(real(res[0.3]), digits = 8))
        flush(stdout)
    end
end
