# Smoke: the public gf_greater_lesser_t(expand = true) path against the dense
# exact leg for the degraded classes (E/G) and one exact-class sanity (C).

using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))
include(joinpath(@__DIR__, "exact_validation.jl"))   # dense helpers, δt, TIMES, TRUNC

function smoke(name, x, couples, L; εd, U, imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd, U)
    occs = ntuple(i -> i == spec.imp_site ? (imp_up, imp_down) :
                  (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
    ψ0 = zeros(ComplexF64, 4^L)
    ψ0[_flat(L, occs)] = 1.0
    E0 = real(H[_flat(L, occs), _flat(L, occs)])
    χ_up = _insert_vector(ψ0, :up, 1, L, spec.imp_site)
    χ_dn = _insert_vector(ψ0, :down, 1, L, spec.imp_site)
    η_up = _insert_vector(ψ0, :up, -1, L, spec.imp_site)
    η_dn = _insert_vector(ψ0, :down, -1, L, spec.imp_site)
    overlaps, = _exact_overlaps(H, Any[χ_up, χ_dn, η_up, η_dn], TIMES)
    gt_exact = [-im * exp(im * E0 * t) * (overlaps[1][i] + overlaps[2][i])
                for (i, t) in enumerate(TIMES)]
    lt_exact = [im * exp(-im * E0 * t) * conj(overlaps[3][i] + overlaps[4][i])
                for (i, t) in enumerate(TIMES)]

    for expand in (false, true)
        gt, lt = gf_greater_lesser_t(x, TIMES; imp_up, imp_down, δt, expand)
        println(name, " expand=$expand: Δgt = ", round(maximum(abs.(gt - gt_exact)); digits = 6),
                "   Δlt = ", round(maximum(abs.(lt - lt_exact)); digits = 6))
    end
    return nothing
end

xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
spec = FermionicTCMPS.chain_spec(x)
wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
imp = spec.imp_site
couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in 1:length(wsl)],
               [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in 1:length(wsr)])

xc = SISBD(discretefermionicvacuum([0.4, 0.9], [0.55, 0.35]; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
_, vc1, vc2 = FermionicTCMPS.thermofield_couplings(xc.bath)
smoke("C single 2mode vacuum U=1.7", xc,
      [(2k, 2k + 1, [0.4, 0.9][k], xc.J * vc1[k], xc.J * vc2[k]) for k in 1:2], 5;
      εd = 0.2, U = 1.7)
smoke("E double 1+1mode finiteβ U=1.7", x, couples, 5; εd = -0.3, U = 1.7)
smoke("G double 1+1mode finiteβ U=1.7 hot imp", x, couples, 5;
      εd = -0.3, U = 1.7, imp_up = 1)
