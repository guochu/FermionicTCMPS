# Smoke: the public gf_greater_lesser_t(expand = true) path against the dense
# exact leg for the degraded vacuum classes (C) and the finite-β classes (E/G,
# plain pipeline only — the expansion is vacuum-only under the exact-image
# construction).

using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))
include(joinpath(@__DIR__, "exact_validation.jl"))   # δt, TIMES, TRUNC

function smoke(name, x, couples, L; εd, U, imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd, U)
    modes = [(c[1], c[2], c[6]) for c in couples]
    ψ0 = dense_tfd_vector(spec, modes, (imp_up, imp_down))
    gt_exact, lt_exact = dense_gf_sandwich(H, ψ0, spec.imp_site, TIMES; imp_up, imp_down)

    hot = any(c -> !iszero(c[6]), couples)
    for expand in (hot ? (false,) : (false, true))
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
wsl, vl1, vl2, nsl = FermionicTCMPS.thermofield_couplings(x.bathl)
wsr, vr1, vr2, nsr = FermionicTCMPS.thermofield_couplings(x.bathr)
imp = spec.imp_site
couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k]) for k in 1:length(wsl)],
               [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k]) for k in 1:length(wsr)])

xc = SISBD(discretefermionicvacuum([0.4, 0.9], [0.55, 0.35]; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
_, vc1, vc2, nsc = FermionicTCMPS.thermofield_couplings(xc.bath)
smoke("C single 2mode vacuum U=1.7", xc,
      [(2k, 2k + 1, [0.4, 0.9][k], xc.J * vc1[k], xc.J * vc2[k], nsc[k]) for k in 1:2], 5;
      εd = 0.2, U = 1.7)
smoke("E double 1+1mode finiteβ U=1.7", x, couples, 5; εd = -0.3, U = 1.7)
smoke("G double 1+1mode finiteβ U=1.7 hot imp", x, couples, 5;
      εd = -0.3, U = 1.7, imp_up = 1)
