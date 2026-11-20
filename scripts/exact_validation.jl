# ==============================================================================
# Exact-vs-sweep validation in the chain frame (small discretized baths).
#
# Three independent legs per config:
#   (i)  dense exact evolution of the mechanical-JW Hamiltonian (assembled by
#        test/anchor_util.jl — fully independent of the MPO code path)
#   (ii) DMRG.exact_timeevolution on OUR MPO (ExactMPS = lossless state, Lanczos
#        exponentiation of the full effective Hamiltonian; memory = state-vector
#        size, NOT 4^L × 4^L)
#   (iii) the production sweep-TDVP pipeline
# (i) ≈ (ii) validates the MPO assembly through a different engine path;
# (iii) vs (i) quantifies the sweep error per config class.
#
# All configs are L = 5 (chain Hilbert dim 1024) — trivial memory.
#
#   julia --project=. scripts/exact_validation.jl
# ==============================================================================

using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))

Random.seed!(20261119)

const δt = 0.02
const TIMES = collect(0.0:0.1:1.0)
const TRUNC = truncdimcutoff(D = 256, ϵ = 1e-10)

# ---------- flat-vector helpers on the mechanical-JW dense Hamiltonian ----------

"flat index of a particle-occupation tuple (site 1 most significant, ↑ + 2↓ per site)"
_flat(L, occs) = 1 + sum((occs[s][1] + 2 * occs[s][2]) * 4^(L - s) for s in 1:L)

"occupation tuple of a flat basis index"
_occ_of(L, idx) = ntuple(s -> begin
    d = (idx - 1) ÷ 4^(L - s) % 4
    (d == 1 || d == 3 ? 1 : 0, d >= 2 ? 1 : 0)
end, L)

"χ = d†_σ|ψ⟩ or η = d_σ|ψ⟩ as a flat vector; nothing if the insertion kills ψ"
function _insert_vector(ψvec::Vector{ComplexF64}, spin::Symbol, delta::Int, L::Int, imp::Int)
    tbl = spin === :up ? (delta == 1 ? FermionicTCMPS._adag_up : FermionicTCMPS._a_up) :
          (delta == 1 ? FermionicTCMPS._adag_down : FermionicTCMPS._a_down)
    step = spin === :up ? 4^(L - imp) : 2 * 4^(L - imp)
    out = zeros(ComplexF64, length(ψvec))
    found = false
    for i in eachindex(ψvec)
        iszero(ψvec[i]) && continue
        q_in = _occ_of(L, i)[imp]   # impurity (n↑, n↓) pair
        q_out = spin === :up ? (q_in[1] + delta, q_in[2]) : (q_in[1], q_in[2] + delta)
        ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || continue
        amp = FermionicTCMPS._amplitude(tbl, q_out, q_in)
        amp === nothing && continue
        out[i + delta * step] += amp * ψvec[i]
        found = true
    end
    return found ? out : nothing
end

"survival overlaps ⟨v(0)|v(t)⟩ for a set of inserted vectors (eigendecomposition)"
function _exact_overlaps(H::Matrix{Float64}, inserted::Vector, times)
    F = eigen(Hermitian(H))
    result = Vector{Vector{ComplexF64}}(undef, length(inserted))
    for (i, vec) in enumerate(inserted)
        if vec === nothing
            result[i] = [zero(ComplexF64) for _ in times]
            continue
        end
        g = F.vectors' * vec
        result[i] = [sum(abs2(g[k]) * exp(-im * F.values[k] * t) for k in eachindex(g)) for t in times]
    end
    return result, F
end

"exact ⟨n_imp↑(t)⟩ — the operator is diagonal in the FLAT basis, so evolve the
flat amplitudes through the eigendecomposition before contracting"
function _nd_exact(H::Matrix{Float64}, ψ0::Vector{ComplexF64}, times, L::Int, imp::Int)
    F = eigen(Hermitian(H))
    g = F.vectors' * ψ0
    nd_diag = [((i - 1) ÷ 4^(L - imp)) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:size(H, 1)]
    out = zeros(Float64, length(times))
    for (i, t) in enumerate(times)
        ψt = F.vectors * (exp.(-im .* F.values .* t) .* g)
        out[i] = real(sum(conj(ψt) .* nd_diag .* ψt))
    end
    return out
end

# ---------- sweep side ----------

function _sweep_nd(x; imp_up, imp_down)
    h = mpo(x)
    occ = occupation_op(x; spin = :up)
    ψ = separable_state(x; imp_up, imp_down)
    cache = environments(h, ψ)
    alg = FermionicTCMPS._tdvp2(δt, TRUNC)
    out = Dict{Float64,Float64}()
    for n in 1:round(Int, last(TIMES) / δt)
        sweep!(cache, alg)
        t = n * δt
        (abs(t / 0.1 - round(t / 0.1)) < 1e-9) && (out[round(t; digits = 6)] = real(expectationvalue(occ, ψ)))
    end
    return out
end

# ---------- exact_timeevolution leg (our MPO, lossless engine) ----------

function _et_gt(x, times; imp_up = 0, imp_down = 0)
    h = mpo(x)
    ψ0m = separable_state(x; imp_up, imp_down)
    E = real(expectationvalue(h, ψ0m))
    gt = zeros(ComplexF64, length(times))
    for (spin, nσ) in ((:up, imp_up), (:down, imp_down))
        nσ == 1 && continue   # greater needs an empty spin channel
        χ0 = ExactMPS(FermionicTCMPS._shifted_state(x, imp_up, imp_down, spin, 1))
        env = environments(h, χ0)
        c = χ0.center
        snapshot = reshape(convert(Array, χ0[c]), :)
        prev = 0
        for (i, t) in enumerate(times)
            n = round(Int, t / δt)
            for _ in (prev + 1):n
                exact_timeevolution!(env, -im * δt; ishermitian = true)
            end
            prev = n
            cur = reshape(convert(Array, env.state[c]), :)
            # the two folds of one config share their side isomorphisms, so the
            # center-tensor inner product IS the state overlap (the generic MPS
            # dot short-circuits on folded ExactMPS states)
            gt[i] += -im * exp(im * E * t) * dot(snapshot, cur)
        end
    end
    return gt
end

# ---------- config runner ----------

"ET leg of the lesser GF: η = d|ψ₀⟩ chains, only occupied spin channels contribute"
function _et_lt(x, times; imp_up = 0, imp_down = 0)
    h = mpo(x)
    ψ0m = separable_state(x; imp_up, imp_down)
    E = real(expectationvalue(h, ψ0m))
    lt = zeros(ComplexF64, length(times))
    for (spin, nσ) in ((:up, imp_up), (:down, imp_down))
        nσ == 0 && continue   # lesser needs an occupied spin channel
        η0 = ExactMPS(FermionicTCMPS._shifted_state(x, imp_up, imp_down, spin, -1))
        env = environments(h, η0)
        c = η0.center
        snapshot = reshape(convert(Array, η0[c]), :)
        prev = 0
        for (i, t) in enumerate(times)
            n = round(Int, t / δt)
            for _ in (prev + 1):n
                exact_timeevolution!(env, -im * δt; ishermitian = true)
            end
            prev = n
            cur = reshape(convert(Array, env.state[c]), :)
            lt[i] += im * exp(-im * E * t) * conj(dot(snapshot, cur))
        end
    end
    return lt
end

function run_config(name, x, couples, L; εd, U, imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd, U)
    size(H, 1) == 4^L || error("dim mismatch: $(size(H, 1)) vs $(4^L)")

    occs = ntuple(i -> i == spec.imp_site ? (imp_up, imp_down) :
                  (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
    idx0 = _flat(L, occs)
    ψ0 = zeros(ComplexF64, 4^L)
    ψ0[idx0] = 1.0
    E0 = real(H[idx0, idx0])

    χ_up = _insert_vector(ψ0, :up, 1, L, spec.imp_site)
    χ_dn = _insert_vector(ψ0, :down, 1, L, spec.imp_site)
    η_up = _insert_vector(ψ0, :up, -1, L, spec.imp_site)
    η_dn = _insert_vector(ψ0, :down, -1, L, spec.imp_site)
    overlaps, = _exact_overlaps(H, Any[χ_up, χ_dn, η_up, η_dn], TIMES)

    gt_exact = [-im * exp(im * E0 * t) * (overlaps[1][i] + overlaps[2][i])
                for (i, t) in enumerate(TIMES)]
    lt_exact = [im * exp(-im * E0 * t) * conj(overlaps[3][i] + overlaps[4][i])
                for (i, t) in enumerate(TIMES)]
    nd_exact = _nd_exact(H, ψ0, TIMES, L, spec.imp_site)

    gt_et = _et_gt(x, TIMES; imp_up, imp_down)
    lt_et = _et_lt(x, TIMES; imp_up, imp_down)
    gt_sw, lt_sw = gf_greater_lesser_t(x, TIMES; δt = δt, trunc = TRUNC, imp_up, imp_down)
    nd_sw = _sweep_nd(x; imp_up, imp_down)

    fmt(v) = round(v, digits = 6)
    fmt2(v) = round(v, digits = 10)
    Δgt_et = maximum(abs.(gt_et - gt_exact))
    Δlt_et = maximum(abs.(lt_et - lt_exact))
    Δgt_sw = maximum(abs.(gt_sw - gt_exact))
    Δlt_sw = maximum(abs.(lt_sw - lt_exact))
    Δnd_sw = maximum(abs(get(nd_sw, round(t; digits = 6), 0.0) - nd_exact[i])
                     for (i, t) in enumerate(TIMES) if t > 0)
    println("== ", name, "  (imp = ", (imp_up, imp_down), ") ==")
    println("   |gt| exact at (0.1, 1.0): ", fmt(abs(gt_exact[2])), " / ", fmt(abs(gt_exact[end])))
    println("   Δgt  exact_timeevolution(our MPO) vs dense: ", fmt(Δgt_et), "  (raw ", fmt2(Δgt_et), ")")
    println("   Δlt  exact_timeevolution(our MPO) vs dense: ", fmt(Δlt_et), "  (raw ", fmt2(Δlt_et), ")")
    println("   Δgt  sweep-TDVP vs dense:                   ", fmt(Δgt_sw))
    println("   Δlt  sweep-TDVP vs dense:                   ", fmt(Δlt_sw))
    println("   Δn   sweep-TDVP vs dense:                   ", fmt(Δnd_sw))
    return (Δgt_et = Δgt_et, Δlt_et = Δlt_et, Δgt_sw = Δgt_sw, Δlt_sw = Δlt_sw, Δnd_sw = Δnd_sw)
end

function main()
    # A. single 2-mode finite β, free — G3-like sanity
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 0.0, εd = -0.3, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("A single 2mode finiteβ U=0", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 0.0)

    # B. single 2-mode finite β, interacting — the G4 failing class
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 1.7, εd = -0.3, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("B single 2mode finiteβ U=1.7", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 1.7)

    # C. single 2-mode vacuum above μ, interacting — the T0 class, multi-mode
    wsc = [0.4, 0.9]; fsc = [0.55, 0.35]
    x = SISBD(discretefermionicvacuum(wsc, fsc; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("C single 2mode vacuum-above U=1.7", x,
        [(2k, 2k + 1, wsc[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(wsc)], 5;
        εd = 0.2, U = 1.7)

    # D. single 2-mode straddling vacuum, interacting — the sharp-sea class
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.7, εd = -0.3, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("D single 2mode straddling vacuum U=1.7", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 1.7)

    # E. double bath 1+1 mode finite β, interacting
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
    imp = spec.imp_site
    couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in 1:length(wsl)],
                   [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in 1:length(wsr)])
    run_config("E double 1+1mode finiteβ U=1.7", x, couples, 5; εd = -0.3, U = 1.7)

    # F. single 2-mode vacuum above μ, interacting, HOT impurity — the ET legs
    # must stay exact here too (the sweep legs quantify the TDVP2 sharp-state
    # limitation on the down-channel greater / up-channel lesser)
    x = SISBD(discretefermionicvacuum(wsc, fsc; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
    _, v1f, v2f = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("F single 2mode vacuum-above U=1.7 hot imp", x,
        [(2k, 2k + 1, wsc[k], x.J * v1f[k], x.J * v2f[k]) for k in eachindex(wsc)], 5;
        εd = 0.2, U = 1.7, imp_up = 1)

    # G. double bath 1+1 mode finite β, interacting, HOT impurity
    x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    specg = FermionicTCMPS.chain_spec(x)
    impg = specg.imp_site
    run_config("G double 1+1mode finiteβ U=1.7 hot imp", x, couples, 5;
        εd = -0.3, U = 1.7, imp_up = 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
