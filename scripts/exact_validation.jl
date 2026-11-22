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
# All correlators use the sandwich estimator (⟨ψ(t)|d|χ(t)⟩, ⟨η(t)|d|ψ(t)⟩):
# with the thermofield-double fiducial the chain state is not a Hamiltonian
# eigenstate, so the eigenstate phase-reduction of the earlier revision is gone.
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
const TRUNC = truncdimcutoff(D = 256, ϵ = 1e-10, add_back = 10)

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
    ψ = fiducial_state(x; imp_up, imp_down)
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

"ET leg of both GFs via the sandwich estimator: the fiducial ψ and every insertion
chain evolved in lockstep with the lossless engine, sampled as
⟨ψ(t)|d_σ|χ(t)⟩ / ⟨η(t)|d_σ|ψ(t)⟩"
function _et_gfls(x, times; imp_up = 0, imp_down = 0)
    h = mpo(x)
    spec = chain_spec(x)
    dmpo = Dict(:up => site_operator_mpo(spec.spaces, spec.imp_site, _a_up),
                :down => site_operator_mpo(spec.spaces, spec.imp_site, _a_down))
    ψenv = environments(h, ExactMPS(fiducial_state(x; imp_up, imp_down)))
    inserts = Tuple{Symbol,Any,Any}[]   # (kind :gt/:lt, dmpo, env)
    for spin in (:up, :down)
        χ = FermionicTCMPS._inserted_state(x, spin, 1; imp_up, imp_down)
        χ === nothing || push!(inserts, (:gt, dmpo[spin], environments(h, ExactMPS(χ))))
        η = FermionicTCMPS._inserted_state(x, spin, -1; imp_up, imp_down)
        η === nothing || push!(inserts, (:lt, dmpo[spin], environments(h, ExactMPS(η))))
    end
    gt = zeros(ComplexF64, length(times))
    lt = zeros(ComplexF64, length(times))
    prev = 0
    for (i, t) in enumerate(times)
        n = round(Int, t / δt)
        for _ in (prev + 1):n
            exact_timeevolution!(ψenv, -im * δt; ishermitian = true)
            foreach(c -> exact_timeevolution!(c[3], -im * δt; ishermitian = true), inserts)
        end
        prev = n
        ψt = ψenv.state
        for (kind, dm, env) in inserts
            if kind === :gt
                gt[i] += -im * expectation(ψt, dm, env.state)
            else
                lt[i] += im * expectation(env.state, dm, ψt)
            end
        end
    end
    return gt, lt
end

# ---------- config runner ----------

function run_config(name, x, couples, L; εd, U, imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd, U)
    size(H, 1) == 4^L || error("dim mismatch: $(size(H, 1)) vs $(4^L)")
    modes = [(c[1], c[2], c[6]) for c in couples]
    ψ0 = dense_tfd_vector(spec, modes, (imp_up, imp_down))

    gt_exact, lt_exact = dense_gf_sandwich(H, ψ0, spec.imp_site, TIMES; imp_up, imp_down)
    nd_exact = _nd_exact(H, ψ0, TIMES, L, spec.imp_site)

    gt_et, lt_et = _et_gfls(x, TIMES; imp_up, imp_down)
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
    _, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("A single 2mode finiteβ U=0", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k], ns[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 0.0)

    # B. single 2-mode finite β, interacting
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 1.7, εd = -0.3, J = 1.0)
    _, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("B single 2mode finiteβ U=1.7", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k], ns[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 1.7)

    # C. single 2-mode vacuum above μ, interacting — the T0 class, multi-mode
    wsc = [0.4, 0.9]; fsc = [0.55, 0.35]
    x = SISBD(discretefermionicvacuum(wsc, fsc; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
    _, v1c, v2c, nsc = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("C single 2mode vacuum-above U=1.7", x,
        [(2k, 2k + 1, wsc[k], x.J * v1c[k], x.J * v2c[k], nsc[k]) for k in eachindex(wsc)], 5;
        εd = 0.2, U = 1.7)

    # D. single 2-mode straddling vacuum, interacting — the sharp-sea class
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.7, εd = -0.3, J = 1.0)
    _, v1d, v2d, nsd = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("D single 2mode straddling vacuum U=1.7", x,
        [(2k, 2k + 1, ws[k], x.J * v1d[k], x.J * v2d[k], nsd[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 1.7)

    # E. double bath 1+1 mode finite β, interacting
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    wsl, vl1, vl2, nsl = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2, nsr = FermionicTCMPS.thermofield_couplings(x.bathr)
    imp = spec.imp_site
    couplesE = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k]) for k in 1:length(wsl)],
                    [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k]) for k in 1:length(wsr)])
    run_config("E double 1+1mode finiteβ U=1.7", x, couplesE, 5; εd = -0.3, U = 1.7)

    # F. single 2-mode vacuum above μ, interacting, HOT impurity — the ET legs
    # must stay exact here too (the sweep legs quantify the TDVP2 sharp-state
    # limitation on the down-channel greater / up-channel lesser)
    x = SISBD(discretefermionicvacuum(wsc, fsc; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
    _, v1f, v2f, nsf = FermionicTCMPS.thermofield_couplings(x.bath)
    run_config("F single 2mode vacuum-above U=1.7 hot imp", x,
        [(2k, 2k + 1, wsc[k], x.J * v1f[k], x.J * v2f[k], nsf[k]) for k in eachindex(wsc)], 5;
        εd = 0.2, U = 1.7, imp_up = 1)

    # G. double bath 1+1 mode finite β, interacting, HOT impurity
    x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    run_config("G double 1+1mode finiteβ U=1.7 hot imp", x, couplesE, 5;
        εd = -0.3, U = 1.7, imp_up = 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
