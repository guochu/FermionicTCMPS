# ==============================================================================
# FermionicTCMPS side of the GTEMPO crosscheck benchmark.
#
# Run in THIS package's environment:  julia --project=. scripts/tcmps_benchmark.jl
# Writes scripts/benchmark_tcmps.csv. Compare with scripts/compare_gtempo.jl
# against the GTEMPO-side CSV produced by scripts/gtempo_benchmark.jl (which must
# run inside the GTEMPO project environment — the two stacks pin incompatible
# TensorKit versions and cannot share an environment).
#
# Semantics note: the TFD separable-state scheme computes the GREATER function of
# the fixed-impurity sector exactly. GTEMPO's thermal-impurity greater function is
# (1 − n_F(εd)) × the empty-sector greater (the occupied-sector greater vanishes),
# which is the quantity compared here. The lesser / n_d(t) / J(t) of the thermal
# impurity are NOT reproduced by the no-ramp separable scheme (they need the
# Kohn–Santoro equilibration ramp); they are recorded for reference only.
# ==============================================================================

using Random
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum, fermidirac

Random.seed!(20261119)

struct BenchConfig
    id::String
    nbaths::Int
    β::Float64                       # Inf → zero temperature
    ws::Vector{Float64}
    fs::Vector{Float64}
    εd::Float64
    U::Float64
    J::Float64
end

const CONFIGS = [
    BenchConfig("single_b05", 1, 0.5, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("single_b1", 1, 1.0, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("double_b05", 2, 0.5, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("double_b1", 2, 1.0, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("single_T0", 1, Inf, [0.4], [0.55], 0.2, 1.7, 1.0),
]

const δt = 0.02
const TIMES = collect(0.0:0.1:1.0)
const TRUNC = truncdimcutoff(D = 256, ϵ = 1e-10)

function run_config(cfg::BenchConfig)
    K = length(cfg.ws)
    if cfg.nbaths == 1
        bath = isinf(cfg.β) ?
               discretefermionicvacuum(cfg.ws, cfg.fs; μ = 0.0) :
               discretefermionicbath(cfg.ws, cfg.fs; β = cfg.β, μ = 0.0)
        x = SISBD(bath; U = cfg.U, εd = cfg.εd, J = cfg.J)
    else
        hl = cfg.ws[1:2]; hr = cfg.ws[3:3]
        fl = cfg.fs[1:2]; fr = cfg.fs[3:3]
        bathl = isinf(cfg.β) ?
                discretefermionicvacuum(hl, fl; μ = 0.0) :
                discretefermionicbath(hl, fl; β = cfg.β, μ = 0.0)
        bathr = isinf(cfg.β) ?
                discretefermionicvacuum(hr, fr; μ = 0.0) :
                discretefermionicbath(hr, fr; β = cfg.β, μ = 0.0)
        x = SIDBD(bathl, bathr; U = cfg.U, εd = cfg.εd, Jl = cfg.J, Jr = cfg.J)
    end
    gt, lt = gf_greater_lesser_t(x, TIMES; δt = δt, trunc = TRUNC)
    # Thermal-impurity greater, per spin band (GTEMPO's greater is single band):
    # the impurity thermal state of H_imp = εd(n↑+n↓) + U n↑n↓ has
    # ⟨n↑⟩ = (e^{−βεd} + e^{−β(2εd+U)})/Z; only the empty-spin sector contributes
    # to the greater (the occupied-sector greater vanishes), hence
    # G^>_band = (1 − ⟨n↑⟩) · G^>_empty-sector, and our gt is the spin sum (/2).
    nd = _thermal_impurity_occupation(cfg.εd, cfg.U, cfg.β)
    return ((1.0 - nd) .* gt ./ 2.0, lt)
end

function _thermal_impurity_occupation(εd::Real, U::Real, β::Real)
    isinf(β) && return Float64(εd < 0.0)
    w = exp.(-β .* [0.0, εd, εd, 2εd + U])   # empty, up, down, full
    return (w[2] + w[4]) / sum(w)
end

function main(outfile::String = joinpath(@__DIR__, "benchmark_tcmps.csv"))
    open(outfile, "w") do io
        println(io, "config,nbaths,beta,t,gt,lt")
        for cfg in CONFIGS
            gt, lt = run_config(cfg)
            for (i, t) in enumerate(TIMES)
                println(io, cfg.id, ",", cfg.nbaths, ",", isinf(cfg.β) ? "Inf" : cfg.β, ",",
                        t, ",", gt[i], ",", lt[i])
            end
            println(stderr, "done: ", cfg.id)
        end
    end
    println("written: ", outfile)
end

if abspath(PROGRAM_FILE) == @__FILE__; main(); end
