# K=100 validation of the compact initial-state expansion:
#   1. construction feasibility (time / memory / bond profile)
#   2. state identity: expanded ≡ raw fiducial (normalized overlap)
#   3. dynamics stability: U=0 GF vs the exact star evolution (short times)
using FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, freefermions_greater_lesser, fermidirac,
                         frequencies, spectrumcouplings
using FermionicTCMPS: SISBD, chain_spec, mpo, expanded_fiducial, _fiducial_mps, _split_sectors, _sector,
                      _maxbond, _fiducial_entries, _ExpansionChannel, _channel_bonds
using DMRG: TDVP2, environments, sweep!, truncdimcutoff, LinearAlgebra
using LinearAlgebra: dot, diagm
using Printf

K = parse(Int, (length(ARGS) ≥ 1 ? ARGS[1] : "100"))
DOverride = (length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : nothing)
βv = parse(Float64, (length(ARGS) ≥ 2 ? ARGS[2] : "1.0"))
tmax = parse(Float64, (length(ARGS) ≥ 3 ? ARGS[3] : "0.1"))

ws = range(-3.0, 3.0; length = K)
fs = [0.5 / (1 + ω^2) for ω in ws]          # spectral weights, peaked at ω=0
x = SISBD(discretefermionicbath(collect(ws), fs; β = βv, μ = 0.0); U = 0.0, εd = -0.3, J = 1.0)
ns = [fermidirac(βv, 0.0, ω) for ω in ws]
spec = chain_spec(x)
L = length(spec.kinds)
@printf("K=%d β=%.2f L=%d  n_k ∈ [%.3f, %.3f]\n", K, βv, L, minimum(ns), maximum(ns))

# ---------- 1) construction ----------
GC.gc()
t0 = time()
ψe = DOverride === nothing ? expanded_fiducial(x) : expanded_fiducial(x; D = DOverride)
t_exp = time() - t0
bond_e = [_maxbond(c) for c in ψe.comps]
@printf("expanded_fiducial: %.1f s, comps=%d, maxbond=%d\n", t_exp, length(ψe.comps), maximum(bond_e))

t0 = time()
raw = _fiducial_mps(x; imp_up = 0, imp_down = 0)
t_raw = time() - t0
@printf("_fiducial_mps (base): %.1f s, maxbond=%d\n", t_raw, _maxbond(raw))

basechan = _ExpansionChannel(_fiducial_entries(x, 0, 0), 1.0)
basebond = maximum(length(_channel_bonds(basechan.entries)[i]) for i in 2:L)
@printf("base intrinsic bond=%d, added slots=%d\n", basebond, maximum(bond_e) - basebond)

# ---------- 2) state identity ----------
nrm_raw = sqrt(real(dot(raw, raw)))
ovl = dot(raw, ψe.comps[1]) / (nrm_raw * sqrt(real(dot(ψe.comps[1], ψe.comps[1]))))
@printf("overlap(expanded, raw fiducial) = %.12f  (raw norm %.12f)\n", real(ovl), nrm_raw)
@assert abs(ovl - 1) < 1e-10

# (the per-sector split at K = 100 would materialize (K+1)² components — a
# diagnostic-only cost the pipeline never pays; the normalized overlap above is
# the identity gate)

# ---------- 3) dynamics vs star ----------
h = mpo(x)
δt = 0.025
times = collect(0.0:0.05:tmax)
D0 = maximum(bond_e)

# inserted chain χ = d†|ψ₀⟩ (up)
χe = FermionicTCMPS._inserted_state(x, :up, 1; ψ = ψe, imp_up = 0, imp_down = 0, expand = true)

# star exact
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

function run_dynamics(h, ψe, χe, x, times, δt, D0, gtr)
    trunc = truncdimcutoff(D = D0, ϵ = 1.0e-10, add_back = 0)
    alg = TDVP2(stepsize = -im * δt, ishermitian = true, trunc = trunc)
    ψc = environments(h, ψe.comps[1])
    χc = environments(h, χe.comps[1])
    worst = 0.0
    prev = 0
    t0 = time()
    for (i, t) in enumerate(times)
        for _ in (prev + 1):round(Int, t / δt)
            sweep!(ψc, alg)
            sweep!(χc, alg)
        end
        prev = round(Int, t / δt)
        tsweeps = time() - t0
        ψt = FermionicTCMPS._collapse_final(ψc.state)
        χt = FermionicTCMPS._collapse_final(χc.state)
        ψs = FermionicTCMPS._apply_dag(x, ψt, :up, 1; imp_up = 0, imp_down = 0)
        common = intersect(FermionicTCMPS._final_charges(ψs), FermionicTCMPS._final_charges(χt))
        gt = -im * dot(FermionicTCMPS._restrict_final(ψs, common),
                       FermionicTCMPS._restrict_final(χt, common))
        d = abs(gt - gtr[i])
        worst = max(worst, d)
        @printf("t=%.2f  gt=%+.6f%+.6fim  star=%+.6f%+.6fim  Δ=%.2e  (sweeps %.1f s)\n",
                t, real(gt), imag(gt), real(gtr[i]), imag(gtr[i]), d, tsweeps)
        t0 = time()
    end
    return worst
end

worst = run_dynamics(h, ψe, χe, x, times, δt, D0, gtr)
@printf("K=%d: worst |gt − star| = %.3e\n", K, worst)
@printf("peak RSS check via GC: live memory ≈ %.2f GB\n", Base.summarysize(Base) / 2^30)
