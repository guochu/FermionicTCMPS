# K=40 compact-expansion validation under a TIGHT bond cap (D ≤ 100):
#   1. construction + state identity (expanded ≡ raw fiducial)
#   2. short-time evolution through the PUBLIC gf_greater_lesser_t with
#      trunc = D(100) — the sweeps pull the chain down to the cap immediately;
#   3. stability: drift vs the exact star evolution stays bounded, bond ≤ 100.
# Usage: julia --project=. scripts/diag_expand_k40.jl [K] [β] [D] [tmax]
using FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, freefermions_greater_lesser, fermidirac,
                         frequencies, spectrumcouplings
using FermionicTCMPS: SISBD, chain_spec, mpo, expanded_fiducial, _fiducial_mps, _maxbond
using DMRG: truncdimcutoff, LinearAlgebra
using LinearAlgebra: dot, diagm
using Printf

K = parse(Int, (length(ARGS) ≥ 1 ? ARGS[1] : "40"))
βv = parse(Float64, (length(ARGS) ≥ 2 ? ARGS[2] : "1.0"))
D = parse(Int, (length(ARGS) ≥ 3 ? ARGS[3] : "100"))
tmax = parse(Float64, (length(ARGS) ≥ 4 ? ARGS[4] : "0.15"))

ws = collect(range(-3.0, 3.0; length = K))
fs = [0.5 / (1 + ω^2) for ω in ws]
x = SISBD(discretefermionicbath(ws, fs; β = βv, μ = 0.0); U = 0.0, εd = -0.3, J = 1.0)
ns = [fermidirac(βv, 0.0, ω) for ω in ws]
spec = chain_spec(x)
L = length(spec.kinds)
@printf("K=%d β=%.2f D=%d L=%d  n_k ∈ [%.3f, %.3f]\n", K, βv, D, L, minimum(ns), maximum(ns))
flush(stdout)

# ---------- 1) construction + identity ----------
t0 = time()
ψe = expanded_fiducial(x)
t_exp = time() - t0
@printf("expanded_fiducial: %.1f s, comps=%d, maxbond=%d\n", t_exp, length(ψe.comps), maximum(_maxbond.(ψe.comps)))
flush(stdout)

raw = _fiducial_mps(x; imp_up = 0, imp_down = 0)
nrm_raw = sqrt(real(dot(raw, raw)))
ovl = dot(raw, ψe.comps[1]) / (nrm_raw * sqrt(real(dot(ψe.comps[1], ψe.comps[1]))))
@printf("overlap(expanded, raw fiducial) = %.12f\n", real(ovl))
@assert abs(ovl - 1) < 1e-10
flush(stdout)

# ---------- 2) short-time GF through the public pipeline, bond capped ----------
δt = 0.025
times = collect(0.0:0.05:tmax)
trunc = truncdimcutoff(D = D, ϵ = 1.0e-10, add_back = 0)
t0 = time()
gt, lt = gf_greater_lesser_t(x, times; δt, expand = true, trunc)
@printf("gf pipeline: %.1f s\n", time() - t0)
flush(stdout)

# star exact (per spin; the pipeline GT is spin-summed = 2×)
wsr = frequencies(x.bath.f)
V = spectrumcouplings(x.bath.f)
hu = zeros(1 + K, 1 + K)
hu[1, 1] = x.εd
for k in 1:K
    hu[1, 1+k] = x.J * V[k]
    hu[1+k, 1] = hu[1, 1+k]
    hu[1+k, 1+k] = wsr[k]
end
ρ0 = diagm(vcat([0.0], [fermidirac(βv, 0.0, ω) for ω in wsr]))
gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)

worst = 0.0
for (i, t) in enumerate(times)
    d = abs(gt[i] - 2gtr[i])
    worst = max(worst, d)
    @printf("t=%.2f  gt=%+.6f%+.6fim  star2=%+.6f%+.6fim  Δ=%.2e\n",
            t, real(gt[i]), imag(gt[i]), real(2gtr[i]), imag(2gtr[i]), d)
    flush(stdout)
end
@printf("K=%d D=%d: worst |gt − 2·star| = %.3e\n", K, D, worst)
@printf("live memory ≈ %.2f GB\n", Base.summarysize(Base) / 2^30)
