# ==============================================================================
# Compare the two benchmark CSVs (run after both benchmark scripts):
#
#   julia --project=. scripts/compare_gtempo.jl
#
# Tolerances per plan G4: ≤ 1e-2 pointwise on the greater function.
# Exits non-zero on failure.
# ==============================================================================

using DelimitedFiles

const TOL = 1.0e-2

function read_csv(path)
    isfile(path) || error("missing $path — run the corresponding benchmark script first")
    rows = readlines(path)
    data = Dict{String,Dict{Float64,ComplexF64}}()
    for r in rows[2:end]
        isempty(r) && continue
        f = split(r, ",")
        cfg, β, t, v = f[1], f[2], parse(Float64, f[3]), parse(ComplexF64, f[end])
        get!(get!(data, cfg, Dict{Float64,ComplexF64}()), t, v)
    end
    return data
end

function main()
    ours = read_csv(joinpath(@__DIR__, "benchmark_tcmps.csv"))
    theirs = read_csv(joinpath(@__DIR__, "benchmark_gtempo.csv"))
    worst = 0.0
    failures = Tuple{String,Float64,Float64,Float64}[]
    for cfg in sort!(collect(keys(ours)))
        haskey(theirs, cfg) || (println("SKIP (no GTEMPO data): ", cfg); continue)
        for (t, v) in sort(collect(ours[cfg]))
            t == 0 && continue                          # GTEMPO grid starts at δt
            haskey(theirs[cfg], t) || continue
            Δ = abs(v - theirs[cfg][t])
            worst = max(worst, Δ)
            Δ > TOL && push!(failures, (cfg, t, abs(v), abs(theirs[cfg][t])))
        end
    end
    println("worst |Δgt| = ", worst)
    if isempty(failures)
        println("G4 PASS (tolerance ", TOL, ")")
    else
        println("G4 FAIL: ", length(failures), " points above tolerance:")
        for f in failures[1:min(end, 20)]
            println("   ", f[1], "  t = ", f[2], "  |ours| = ", f[3], "  |gtempo| = ", f[4])
        end
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__; main(); end
