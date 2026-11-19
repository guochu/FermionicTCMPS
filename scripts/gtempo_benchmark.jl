# ==============================================================================
# GTEMPO side of the crosscheck benchmark.
#
# MUST run inside the GTEMPO project environment (the GTEMPO stack pins
# TensorKit 0.17 while FermionicTCMPS pins 0.13 — the two cannot share an env):
#
#   julia --project=/home/guochu/Documents/Missile/GTEMPO \
#         /home/guochu/Documents/Missile/FermionicTCMPS/scripts/gtempo_benchmark.jl
#
# Writes benchmark_gtempo.csv next to this file. The config literals below are
# duplicated from scripts/tcmps_benchmark.jl by design (双侧同源字面量).
# ==============================================================================

using GTEMPO
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum

const δt = 0.02
const NT = 50                      # times δt·(1:NT); compared against TCMPS t = 0.1:0.1:1.0

struct BenchConfig
    id::String
    β::Float64
    ws::Vector{Float64}
    fs::Vector{Float64}
    εd::Float64
    U::Float64
    J::Float64
end

# GTEMPO convention: `bands` are spin bands (AndersonIM has 2) and the bath
# spectrum is per band. The double-bath config concatenates the two baths'
# modes (the influence functional only sees the total hybridization).
const CONFIGS = [
    BenchConfig("single_b05", 0.5, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("single_b1", 1.0, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("double_b05", 0.5, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("double_b1", 1.0, [-0.9, -0.3, 0.5], [0.35, 0.45, 0.5], -0.3, 1.7, 1.0),
    BenchConfig("single_T0", Inf, [0.4], [0.55], 0.2, 1.7, 1.0),
]

function run_config(cfg::BenchConfig)
    bath = isinf(cfg.β) ?
           discretefermionicvacuum(cfg.ws, cfg.fs; μ = 0.0) :
           discretefermionicbath(cfg.ws, cfg.fs; β = cfg.β, μ = 0.0)
    lat = GrassmannLattice(N = NT, δt = δt, contour = :real, bands = 2)
    kernel = correlationfunction(bath, lat)
    model = AndersonIM(U = cfg.U, ϵ_d = cfg.εd)
    cn = coupling_norm(bath)
    Is, _ = hybriddynamics(lat, kernel, PartialIF())
    Is isa GTEMPO.GMPS && (Is = (Is,))
    initial = Matrix{ComplexF64}(fock_thermalstate(model, cfg.β).data)
    K = sysdynamics(lat, model; order = 2, coupling_norm = cn, initial_state = initial)
    Z = integrate(K, Is...)
    gt = [greater(lat, model, kernel, Is...; i = i, j = 0, Z = Z, order = 2,
                  coupling_norm = cn, initial_state = initial) for i in 1:NT]
    return gt
end

function main(outfile::String = joinpath(@__DIR__, "benchmark_gtempo.csv"))
    open(outfile, "w") do io
        println(io, "config,beta,t,gt")
        for cfg in CONFIGS
            gt = run_config(cfg)
            for i in 1:NT
                println(io, cfg.id, ",", isinf(cfg.β) ? "Inf" : cfg.β, ",", i * δt, ",", gt[i])
            end
            println(stderr, "done: ", cfg.id)
        end
    end
    println("written: ", outfile)
end

if abspath(PROGRAM_FILE) == @__FILE__; main(); end
