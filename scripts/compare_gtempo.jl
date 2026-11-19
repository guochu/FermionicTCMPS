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
        # TCMPS csv: config,nbaths,beta,t,gt,lt — the GTEMPO csv: config,beta,t,gt
        if basename(path) == "benchmark_tcmps.csv"
            cfg, t, v = f[1], parse(Float64, f[4]), parse(ComplexF64, f[5])
        else
            cfg, t, v = f[1], parse(Float64, f[3]), parse(ComplexF64, f[end])
        end
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
        for (t, v) in sort(collect(ours[cfg]), by = first)
            t == 0 && continue                          # GTEMPO grid starts at δt
            haskey(theirs[cfg], t) || continue
            Δ = abs(v - theirs[cfg][t])
            worst = max(worst, Δ)
            Δ > TOL && push!(failures, (cfg, t, abs(v), abs(theirs[cfg][t])))
        end
    end
    println("worst |Δgt| = ", worst)
    # accuracy envelope: the largest t such that ALL configs stay within the
    # tolerance for every grid point in (0, t]
    all_times = sort(collect(union([Set(keys(theirs[c])) for c in keys(theirs)]...)))
    env = 0.0
    for t in all_times
        t == 0 && continue
        ok = true
        for cfg in keys(ours)
            (haskey(ours[cfg], t) && haskey(theirs[cfg], t)) || continue
            abs(ours[cfg][t] - theirs[cfg][t]) > TOL && (ok = false; break)
        end
        ok || break
        env = t
    end
    println("accuracy envelope (all configs within ", TOL, "): t ≤ ", env)
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
