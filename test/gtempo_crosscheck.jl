# ==============================================================================
# G4 crosscheck gate.
#
# The full comparison runs as two separate processes (see scripts/): the GTEMPO
# stack pins TensorKit 0.17 and this package pins 0.13, so they cannot share an
# environment. This test consumes the GTEMPO-side CSV produced by
#   julia --project=<GTEMPO> scripts/gtempo_benchmark.jl
# and compares it against the TCMPS-side values recomputed live (seconds).
# Skips when the GTEMPO-side CSV is absent.
# ==============================================================================

const GTEMPO_CSV = joinpath(@__DIR__, "..", "scripts", "benchmark_gtempo.csv")

if isfile(GTEMPO_CSV) && filesize(GTEMPO_CSV) > 60
    @testset "G4 crosscheck vs GTEMPO (greater, thermal impurity)" begin
        rows = filter(!isempty, readlines(GTEMPO_CSV))
        gt_of = Dict{Tuple{String,Float64},ComplexF64}()
        times_of = Dict{String,Vector{Float64}}()
        for r in rows[2:end]
            f = split(r, ",")
            cfg, t, v = f[1], parse(Float64, f[3]), parse(ComplexF64, f[end])
            gt_of[(cfg, t)] = v
            push!(get!(times_of, cfg, Float64[]), t)
        end

        using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum, fermidirac
        δt = 0.02
        trunc = truncdimcutoff(D = 256, ϵ = 1e-10)
        for (cfg, bath, εd) in (
            ("single_b05", discretefermionicbath([-0.9, -0.3, 0.5], [0.35, 0.45, 0.5]; β = 0.5, μ = 0.0), -0.3),
            ("single_T0", discretefermionicvacuum([0.4], [0.55]; μ = 0.0), 0.2),
        )
            haskey(times_of, cfg) || continue
            x = SISBD(bath; U = 1.7, εd = εd, J = 1.1)
            times = sort(unique(times_of[cfg]))
            gt, = gf_greater_lesser_t(x, times; δt, trunc)
            β = FermionicTCMPS._inverse_temperature(bath)
            nd = isinf(β) ? Float64(εd < 0.0) : fermidirac(β, 0.0, εd)
            ours = (1.0 - nd) .* gt
            worst = maximum(abs(ours[i] - gt_of[(cfg, t)]) for (i, t) in enumerate(times) if t > 0)
            @test worst < 1e-2
        end
    end
else
    println("G4 crosscheck skipped: run scripts/gtempo_benchmark.jl in the GTEMPO env first")
end
