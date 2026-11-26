# Measure the PLAIN-route errors of the two G5 testsets that were expand-based,
# to set their post-removal gates from data (not guesses):
#   (a) L1 SIDBD double vacuum-sea GF vs dense ED (δt=0.01, t ≤ 0.3)
#   (b) L2 U=0 vacuum/sea configs, hot impurity: occupation + 2·Im⟨J⟩ vs the
#       strict solution (δt=0.01, t ≤ 0.4), plain fiducial + default trunc
using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicvacuum
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))

Random.seed!(20261126)
const PLAIN_TRUNC = truncdimcutoff(D = 256, ϵ = 1e-10, add_back = 10)

# ---------- (a) L1 SIDBD GF ----------
let δt = 0.01, times = collect(0.05:0.05:0.3)
    xl = discretefermionicvacuum([-0.8], [0.6]; μ = 0.0)
    xr = discretefermionicvacuum([0.5], [0.45]; μ = 0.0)
    x = SIDBD(xl, xr; U = 1.0, εd = -0.3, Jl = 1.2, Jr = 0.9)
    spec = FermionicTCMPS.chain_spec(x); imp = spec.imp_site
    wsl, vl1, vl2, nsl = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2, nsr = FermionicTCMPS.thermofield_couplings(x.bathr)
    couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k]) for k in eachindex(wsl)],
                   [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k]) for k in eachindex(wsr)])
    H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
    modes = [(c[1], c[2], c[6]) for c in couples]
    ψ0 = dense_tfd_vector(spec, modes, (0, 0))
    gt_d, lt_d = dense_gf_sandwich(H, ψ0, imp, times)
    gt, lt = gf_greater_lesser_t(x, times; δt, trunc = PLAIN_TRUNC)
    println("(a) L1 SIDBD plain GF: max|gt−gt_d| = ", maximum(abs.(gt .- gt_d)),
            "   max|lt−lt_d| = ", maximum(abs.(lt .- lt_d)))
    flush(stdout)
end

# ---------- (b) L2 occupation + current ----------
let δt = 0.01, times = collect(0.1:0.1:0.4)
    configs = Tuple{String,Any,Vector{Tuple{Int,Int,Float64,Float64,Float64,Float64}}}[]
    let x = SISBD(discretefermionicvacuum([0.4], [0.55]; μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
        _, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
        push!(configs, ("vacuum K=1 (above μ)", x, [(2k, 2k + 1, [0.4][k], x.J * v1[k], x.J * v2[k], ns[k])
                                                    for k in eachindex(v1)]))
    end
    let x = SISBD(discretefermionicvacuum([-0.9, 0.5], [0.35, 0.55]; μ = 0.0);
                  U = 0.0, εd = -0.3, J = 1.1)
        _, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
        push!(configs, ("vacuum K=2 (straddling sea)", x, [(2k, 2k + 1, [-0.9, 0.5][k], x.J * v1[k], x.J * v2[k], ns[k])
                                                           for k in eachindex(v1)]))
    end
    let x = SIDBD(discretefermionicvacuum([-0.8], [0.6]; μ = 0.0),
                  discretefermionicvacuum([0.5], [0.45]; μ = 0.0);
                  U = 0.0, εd = -0.3, Jl = 1.0, Jr = 1.0)
        spec = FermionicTCMPS.chain_spec(x); imp = spec.imp_site
        wsl, vl1, vl2, nsl = FermionicTCMPS.thermofield_couplings(x.bathl)
        wsr, vr1, vr2, nsr = FermionicTCMPS.thermofield_couplings(x.bathr)
        couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k]) for k in eachindex(wsl)],
                       [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k]) for k in eachindex(wsr)])
        push!(configs, ("double vacuum 1+1 (straddling sea)", x, couples))
    end

    for (name, x, couples) in configs
        spec = FermionicTCMPS.chain_spec(x)
        L = length(spec.kinds); imp = spec.imp_site
        H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
        F = eigen(Hermitian(H))
        modes = [(c[1], c[2], c[6]) for c in couples]
        ψ0 = dense_tfd_vector(spec, modes, (1, 0))
        g0 = F.vectors' * ψ0
        nddiag = [(i - 1) ÷ 4^(L - imp) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]
        n_strict(t) = begin
            ψt = F.vectors * (exp.(-im .* F.values .* t) .* g0)
            real(ψt' * Diagonal(nddiag) * ψt)
        end
        hfd = 0.005
        jref = [(n_strict(t - hfd) - n_strict(t + hfd)) / (2 * hfd) for t in times]

        h = mpo(x)
        jop = current_op(x; spin = :both)
        occ = occupation_op(x; spin = :up)
        ψ = fiducial_state(x; imp_up = 1, imp_down = 0)
        cache = environments(h, ψ)
        alg = FermionicTCMPS._tdvp2(δt, PLAIN_TRUNC)
        steps = round.(Int, times ./ δt)
        prev = 0
        dn = 0.0; dj = 0.0
        for (i, t) in enumerate(times)
            for _ in (prev + 1):steps[i]; sweep!(cache, alg); end
            prev = steps[i]
            nv = real(expectationvalue(occ, ψ))
            jd = 2 * imag(expectationvalue(jop, ψ))
            dn = max(dn, abs(nv - n_strict(t)))
            dj = max(dj, abs(jd - jref[i]))
        end
        println("(b) ", rpad(name, 34), "  max|Δn| = ", dn, "   max|ΔJ| = ", dj)
        flush(stdout)
    end
end
