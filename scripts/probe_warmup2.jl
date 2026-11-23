# warm-up across the config classes: vacuum-sea (D), hot imp (F), SIDBD finite-β (E/G class)
using Random, TensorKit, DMRG, FermionicTCMPS
import ImpurityModelBase
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum, freefermions_greater_lesser
using LinearAlgebra
const F = FermionicTCMPS
include(joinpath(pwd(), "test", "anchor_util.jl"))

function star_helpers(x)
    function _sh(x::F.SISBD)
        ws = ImpurityModelBase.frequencies(x.bath.f)
        V = ImpurityModelBase.spectrumcouplings(x.bath.f)
        K = length(ws)
        hu = zeros(1 + K, 1 + K)
        hu[1, 1] = x.εd
        for k in 1:K
            hu[1, 1+k] = x.J * V[k]; hu[1+k, 1] = hu[1, 1+k]; hu[1+k, 1+k] = ws[k]
        end
        hu
    end
    _sh, nothing
end

function run_case(name, x, times, δt; imp_up = 0, imp_down = 0)
    # dense reference
    spec = F.chain_spec(x)
    L = length(spec.kinds)
    wsall = Float64[]; couples = Tuple{Int,Int,Float64,Float64,Float64,Float64}[]
    if x isa F.SISBD
        wsb, v1b, v2b, nsb = F.thermofield_couplings(x.bath)
        for k in eachindex(wsb)
            push!(couples, (2k, 2k + 1, wsb[k], x.J * v1b[k], x.J * v2b[k], nsb[k]))
        end
    else
        wsl, vl1, vl2, nsl = F.thermofield_couplings(x.bathl)
        wsr, vr1, vr2, nsr = F.thermofield_couplings(x.bathr)
        imp = spec.imp_site
        for k in eachindex(wsl); push!(couples, (2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k])); end
        for k in eachindex(wsr); push!(couples, (imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k])); end
    end
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    Fh = eigen(Hermitian(H))
    modes = [(c[1], c[2], c[6]) for c in couples]
    ψ0 = dense_tfd_vector(spec, modes, (imp_up, imp_down))
    g0 = Fh.vectors' * ψ0
    gt_exact, lt_exact = dense_gf_sandwich(H, ψ0, spec.imp_site, times; imp_up, imp_down)

    gt0, lt0 = gf_greater_lesser_t(x, times; δt, imp_up, imp_down)
    gt1, lt1 = gf_greater_lesser_t(x, times; δt, imp_up, imp_down, warmup = :WI)
    gt2, lt2 = gf_greater_lesser_t(x, times; δt, imp_up, imp_down, warmup = :WII)
    gt3, lt3 = gf_greater_lesser_t(x, times; δt, imp_up, imp_down, expand = true)
    println(name, ":")
    println("   plain  : Δgt = ", round(maximum(abs.(gt0 .- gt_exact)), digits = 6),
            "  Δlt = ", round(maximum(abs.(lt0 .- lt_exact)), digits = 6))
    println("   WI     : Δgt = ", round(maximum(abs.(gt1 .- gt_exact)), digits = 6),
            "  Δlt = ", round(maximum(abs.(lt1 .- lt_exact)), digits = 6))
    println("   WII    : Δgt = ", round(maximum(abs.(gt2 .- gt_exact)), digits = 6),
            "  Δlt = ", round(maximum(abs.(lt2 .- lt_exact)), digits = 6))
    println("   expand : Δgt = ", round(maximum(abs.(gt3 .- gt_exact)), digits = 8),
            "  Δlt = ", round(maximum(abs.(lt3 .- lt_exact)), digits = 8))
end

Random.seed!(31)
δt = 0.02
times = collect(0.0:0.1:0.5)

# D class: straddling vacuum (sharp sea), U=1.7
xD = SISBD(discretefermionicvacuum([-0.9, 0.5], [0.35, 0.55]; μ = 0.0); U = 1.7, εd = -0.3)
run_case("D straddling vacuum U=1.7", xD, times, δt)

# F class: vacuum above μ, hot impurity
xF = SISBD(discretefermionicvacuum([0.4, 0.9], [0.55, 0.35]; μ = 0.0); U = 1.7, εd = 0.2)
run_case("F vacuum-above hot imp", xF, times, δt; imp_up = 1)

# E class: double bath finite β
xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
xE = SIDBD(xl, xr; U = 1.7, εd = -0.3)
run_case("E double 1+1 finiteβ", xE, times, δt)
