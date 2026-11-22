# ==============================================================================
# Models: discrete single-impurity Anderson models and their thermo-field chains.
#
# Physics convention: every bath mode k (frequency ω_k, spectral coupling
# V_k = √f_k, thermal occupation n_k = fermidirac(β, μ, ω_k)) is thermofield
# doubled into an interleaved (a₁, a₂) pair of PARTICLE-labelled sites
# representing (f̂₁ₖ, f̂₂ₖ) with c̃ₖ = √(1−n_k) f̂₁ₖ + √n_k f̂₂ₖ. The chain carries
# the EXACT Bogoliubov image of the star Hamiltonian,
#
#   H = εd·n_d + U·n↑n↓
#     + Σ_k [ ω_k(1−n_k)·n_a₁ + ω_k n_k·n_a₂ + ω_k√(n_k(1−n_k))(a₁†a₂ + a₂†a₁) ]
#     + Σ_k [ J·V_k√(1−n_k)(d†a₁ + h.c.) + J·V_k√n_k(d†a₂ + h.c.) ],
#
# and the fiducial state (impurity occupations ⊗ ⊗_k σ_σ(√(1−n_k)|00⟩ + √n_k|11⟩))
# is the exact pure-state image of the thermal star density matrix: every impurity
# observable of the sudden protocol — equal-time occupation/current AND (0,t)
# Green's functions alike — is reproduced exactly (the span of the impurity
# orbital and each mode's thermofield combination w_k = √(1−n)f̂₁ + √n f̂₂† is an
# invariant subspace carrying the star single-particle problem, and the fiducial
# one-particle density matrix restricts to the thermal star data). The earlier
# construction (KS 2012.01424 Eq. (9): unweighted on-sites, particle-hole-labelled
# a₂ branch, product fiducial) only mimics the bath occupation factors (Eq. (7));
# it reproduces the (0,t) correlators of a diagonal star ρ₀ but fails the
# equal-time dynamics at O(t³). Sites are interleaved (a₁, a₂) per mode; the
# impurity sits at the left end (single bath) or in the middle (double bath).
# ==============================================================================

"""Discrete fermionic bath: `discretefermionicbath` (finite β) or `discretefermionicvacuum`."""
const DiscreteFermionicBath = Union{Bath{Fermion,DiscreteSpectrum}, Vacuum{Fermion,DiscreteSpectrum}}

_inverse_temperature(b::Bath) = b.β
_inverse_temperature(::Vacuum) = Inf

"""
    thermofield_couplings(bath) -> (ws, v1, v2, ns)

Thermo-field splitting of a discrete fermionic bath: for each mode frequency ω_k and
spectral coupling V_k = √f_k with thermal occupation n_k, returns the per-branch
couplings `v1_k = V_k·√(1−n_k)` (a₁ branch), `v2_k = V_k·√n_k` (a₂ branch) and the
occupations `ns`.
"""
function thermofield_couplings(bath::DiscreteFermionicBath)
    ws = frequencies(bath.f)
    vs = spectrumcouplings(bath.f)
    β = _inverse_temperature(bath)
    ns = fermidirac.(β, bath.μ, ws)
    v1 = [v * sqrt(1 - n) for (v, n) in zip(vs, ns)]
    v2 = [v * sqrt(n) for (v, n) in zip(vs, ns)]
    return ws, v1, v2, ns
end

"""
    SISBD(bath; U=0, εd=0, J=1)

Single-Impurity Single-Bath discrete Anderson model,
H = εd·n_d + U·n↑n↓ + Σ_j ω_j n_j + Σ_j J·V_j (d†c_j + h.c.),
with the bath encoded by a `discretefermionicbath` / `discretefermionicvacuum`.
Finite temperature is handled internally by the thermo-field doubling described
in the module header; the chain is [d, a₁₁, a₂₁, a₁₂, a₂₂, …].
"""
struct SISBD
    bath::DiscreteFermionicBath
    U::Float64
    εd::Float64
    J::Float64
end
SISBD(bath::DiscreteFermionicBath; U::Real = 0.0, εd::Real = 0.0, J::Real = 1.0) =
    SISBD(bath, Float64(U), Float64(εd), Float64(J))

"""
    SIDBD(bathl, bathr; U=0, εd=0, Jl=1, Jr=1)

Single-Impurity Double-Bath discrete Anderson model with the impurity in the middle:
chain [left bath (a₁/a₂ interleaved), d, right bath (a₁/a₂ interleaved)].
The left bath follows the coefficient-at-the-opening-end convention and the right
bath the coefficient-at-the-closing-end one (both ends are non-shared, so the
placement rule "coefficients at the non-shared end" holds on both sides).
"""
struct SIDBD
    bathl::DiscreteFermionicBath
    bathr::DiscreteFermionicBath
    U::Float64
    εd::Float64
    Jl::Float64
    Jr::Float64
end
SIDBD(bathl::DiscreteFermionicBath, bathr::DiscreteFermionicBath;
      U::Real = 0.0, εd::Real = 0.0, Jl::Real = 1.0, Jr::Real = 1.0) =
    SIDBD(bathl, bathr, Float64(U), Float64(εd), Float64(Jl), Float64(Jr))

# ---------------- finite-state-machine channel description ----------------

struct OpenEnd
    site::Int
    coef::Float64
    amps::OpTable
end
struct CloseEnd
    site::Int
    coef::Float64
    amps::OpTable
end
struct FSMChannel
    name::Symbol
    side::Symbol      # :right — opens left of the impurity or in a right bath;
                      # :left — opens in a left bath
    charge::Charge
    opens::Vector{OpenEnd}
    closes::Vector{CloseEnd}
end

const PhysSpace = typeof(particle_space())

struct ChainSpec
    imp_site::Int
    kinds::Vector{Symbol}             # :particle per chain site (kept for callers)
    spaces::Vector{PhysSpace}         # physical space per chain site
    onsite::Vector{OpTable}           # on-site amplitude table per site
    channels::Vector{FSMChannel}
end

# on-site amplitude tables
_impurity_onsite(εd::Real, U::Real) = table(
    (P_UP, P_UP) => εd,
    (P_DOWN, P_DOWN) => εd,
    (P_FULL, P_FULL) => 2 * εd + U,
)
# branch on-site energy ω·w per spin: a₁ branch w = 1−n, a₂ branch w = n
_branch_onsite(ω::Real, w::Real) =
    iszero(ω * w) ? OpTable([]) :
    table((P_UP, P_UP) => ω * w, (P_DOWN, P_DOWN) => ω * w, (P_FULL, P_FULL) => 2 * ω * w)

"interleaved (a₁ site, a₂ site, n_k) pairs of a model, in chain order"
_mode_pairs(x::SISBD) = begin
    ws, _, _, ns = thermofield_couplings(x.bath)
    [(2k, 2k + 1, ns[k]) for k in eachindex(ws)]
end
_mode_pairs(x::SIDBD) = begin
    wsl, _, _, nsl = thermofield_couplings(x.bathl)
    wsr, _, _, nsr = thermofield_couplings(x.bathr)
    imp = 2length(wsl) + 1
    vcat([(2k - 1, 2k, nsl[k]) for k in eachindex(wsl)],
         [(imp + 2k - 1, imp + 2k, nsr[k]) for k in eachindex(nsr)])
end

function chain_spec(x::SISBD)
    ws, v1, v2, ns = thermofield_couplings(x.bath)
    K = length(ws)
    L = 2K + 1
    kinds = fill(:particle, L)
    spaces = fill(particle_space(), L)
    onsite = Vector{OpTable}(undef, L)
    kinds[1], spaces[1], onsite[1] = :particle, particle_space(), _impurity_onsite(x.εd, x.U)
    for k in 1:K
        s1, s2 = 2k, 2k + 1
        kinds[s1], spaces[s1], onsite[s1] = :particle, particle_space(), _branch_onsite(ws[k], 1 - ns[k])
        kinds[s2], spaces[s2], onsite[s2] = :particle, particle_space(), _branch_onsite(ws[k], ns[k])
    end

    # impurity → branch hoppings (both branches number-conserving now)
    hop_up = FSMChannel(:up_hop, :right, (1, 0), [OpenEnd(1, 1.0, compose(_adag_up, _parity))], CloseEnd[])
    hop_down = FSMChannel(:down_hop, :right, (0, 1), [OpenEnd(1, 1.0, compose(_adag_down, _parity))], CloseEnd[])
    hc_up = FSMChannel(:up_hc, :right, (-1, 0), [OpenEnd(1, 1.0, compose(_parity, _a_up))], CloseEnd[])
    hc_down = FSMChannel(:down_hc, :right, (0, -1), [OpenEnd(1, 1.0, compose(_parity, _a_down))], CloseEnd[])
    pairs = FSMChannel[]
    for k in 1:K
        s1, s2 = 2k, 2k + 1
        push!(hop_up.closes, CloseEnd(s1, x.J * v1[k], _a_up))
        push!(hop_up.closes, CloseEnd(s2, x.J * v2[k], _a_up))
        push!(hop_down.closes, CloseEnd(s1, x.J * v1[k], _a_down))
        push!(hop_down.closes, CloseEnd(s2, x.J * v2[k], _a_down))
        push!(hc_up.closes, CloseEnd(s1, x.J * v1[k], _adag_up))
        push!(hc_up.closes, CloseEnd(s2, x.J * v2[k], _adag_up))
        push!(hc_down.closes, CloseEnd(s1, x.J * v1[k], _adag_down))
        push!(hc_down.closes, CloseEnd(s2, x.J * v2[k], _adag_down))
        # intra-mode pairing ω√(n(1−n))(a₁†a₂ + a₂†a₁), one channel pair per spin
        gk = ws[k] * sqrt(ns[k] * (1 - ns[k]))
        iszero(gk) && continue
        push!(pairs, FSMChannel(:up_pair, :right, (1, 0),
                                [OpenEnd(s1, gk, compose(_adag_up, _parity))],
                                [CloseEnd(s2, 1.0, _a_up)]))
        push!(pairs, FSMChannel(:up_pair_ad, :right, (-1, 0),
                                [OpenEnd(s1, gk, compose(_parity, _a_up))],
                                [CloseEnd(s2, 1.0, _adag_up)]))
        push!(pairs, FSMChannel(:down_pair, :right, (0, 1),
                                [OpenEnd(s1, gk, compose(_adag_down, _parity))],
                                [CloseEnd(s2, 1.0, _a_down)]))
        push!(pairs, FSMChannel(:down_pair_ad, :right, (0, -1),
                                [OpenEnd(s1, gk, compose(_parity, _a_down))],
                                [CloseEnd(s2, 1.0, _adag_down)]))
    end
    channels = vcat(FSMChannel[hop_up, hop_down, hc_up, hc_down], pairs)
    _drop_zero_ends!.(channels)
    return ChainSpec(1, kinds, spaces, onsite, channels)
end

function chain_spec(x::SIDBD)
    wsl, vl1, vl2, nsl = thermofield_couplings(x.bathl)
    wsr, vr1, vr2, nsr = thermofield_couplings(x.bathr)
    Kl, Kr = length(wsl), length(wsr)
    L = 2(Kl + Kr) + 1
    imp = 2Kl + 1
    kinds = fill(:particle, L)
    spaces = fill(particle_space(), L)
    onsite = Vector{OpTable}(undef, L)
    for k in 1:Kl
        s1, s2 = 2k - 1, 2k
        kinds[s1], spaces[s1], onsite[s1] = :particle, particle_space(), _branch_onsite(wsl[k], 1 - nsl[k])
        kinds[s2], spaces[s2], onsite[s2] = :particle, particle_space(), _branch_onsite(wsl[k], nsl[k])
    end
    kinds[imp], spaces[imp], onsite[imp] = :particle, particle_space(), _impurity_onsite(x.εd, x.U)
    for k in 1:Kr
        s1, s2 = imp + 2k - 1, imp + 2k
        kinds[s1], spaces[s1], onsite[s1] = :particle, particle_space(), _branch_onsite(wsr[k], 1 - nsr[k])
        kinds[s2], spaces[s2], onsite[s2] = :particle, particle_space(), _branch_onsite(wsr[k], nsr[k])
    end

    Jl, Jr = x.Jl, x.Jr
    # left bath: coefficient at the opening (bath) end, closes at the impurity
    l_hop_up = FSMChannel(:up_hop, :left, (-1, 0), OpenEnd[], [CloseEnd(imp, 1.0, _adag_up)])
    l_hop_down = FSMChannel(:down_hop, :left, (0, -1), OpenEnd[], [CloseEnd(imp, 1.0, _adag_down)])
    l_hc_up = FSMChannel(:up_hc, :left, (1, 0), OpenEnd[], [CloseEnd(imp, 1.0, _a_up)])
    l_hc_down = FSMChannel(:down_hc, :left, (0, 1), OpenEnd[], [CloseEnd(imp, 1.0, _a_down)])
    # right bath: coefficient at the closing (bath) end, opens at the impurity
    r_hop_up = FSMChannel(:up_hop, :right, (1, 0), [OpenEnd(imp, 1.0, compose(_adag_up, _parity))], CloseEnd[])
    r_hop_down = FSMChannel(:down_hop, :right, (0, 1), [OpenEnd(imp, 1.0, compose(_adag_down, _parity))], CloseEnd[])
    r_hc_up = FSMChannel(:up_hc, :right, (-1, 0), [OpenEnd(imp, 1.0, compose(_parity, _a_up))], CloseEnd[])
    r_hc_down = FSMChannel(:down_hc, :right, (0, -1), [OpenEnd(imp, 1.0, compose(_parity, _a_down))], CloseEnd[])

    pairs = FSMChannel[]
    for k in 1:Kl
        s1, s2 = 2k - 1, 2k
        push!(l_hop_up.opens, OpenEnd(s1, Jl * vl1[k], compose(_parity, _a_up)))
        push!(l_hop_up.opens, OpenEnd(s2, Jl * vl2[k], compose(_parity, _a_up)))
        push!(l_hop_down.opens, OpenEnd(s1, Jl * vl1[k], compose(_parity, _a_down)))
        push!(l_hop_down.opens, OpenEnd(s2, Jl * vl2[k], compose(_parity, _a_down)))
        push!(l_hc_up.opens, OpenEnd(s1, Jl * vl1[k], compose(_adag_up, _parity)))
        push!(l_hc_up.opens, OpenEnd(s2, Jl * vl2[k], compose(_adag_up, _parity)))
        push!(l_hc_down.opens, OpenEnd(s1, Jl * vl1[k], compose(_adag_down, _parity)))
        push!(l_hc_down.opens, OpenEnd(s2, Jl * vl2[k], compose(_adag_down, _parity)))
        gk = wsl[k] * sqrt(nsl[k] * (1 - nsl[k]))
        iszero(gk) && continue
        push!(pairs, FSMChannel(:up_pair, :left, (1, 0),
                                [OpenEnd(s1, gk, compose(_adag_up, _parity))],
                                [CloseEnd(s2, 1.0, _a_up)]))
        push!(pairs, FSMChannel(:up_pair_ad, :left, (-1, 0),
                                [OpenEnd(s1, gk, compose(_parity, _a_up))],
                                [CloseEnd(s2, 1.0, _adag_up)]))
        push!(pairs, FSMChannel(:down_pair, :left, (0, 1),
                                [OpenEnd(s1, gk, compose(_adag_down, _parity))],
                                [CloseEnd(s2, 1.0, _a_down)]))
        push!(pairs, FSMChannel(:down_pair_ad, :left, (0, -1),
                                [OpenEnd(s1, gk, compose(_parity, _a_down))],
                                [CloseEnd(s2, 1.0, _adag_down)]))
    end
    for k in 1:Kr
        s1, s2 = imp + 2k - 1, imp + 2k
        push!(r_hop_up.closes, CloseEnd(s1, Jr * vr1[k], _a_up))
        push!(r_hop_up.closes, CloseEnd(s2, Jr * vr2[k], _a_up))
        push!(r_hop_down.closes, CloseEnd(s1, Jr * vr1[k], _a_down))
        push!(r_hop_down.closes, CloseEnd(s2, Jr * vr2[k], _a_down))
        push!(r_hc_up.closes, CloseEnd(s1, Jr * vr1[k], _adag_up))
        push!(r_hc_up.closes, CloseEnd(s2, Jr * vr2[k], _adag_up))
        push!(r_hc_down.closes, CloseEnd(s1, Jr * vr1[k], _adag_down))
        push!(r_hc_down.closes, CloseEnd(s2, Jr * vr2[k], _adag_down))
        gk = wsr[k] * sqrt(nsr[k] * (1 - nsr[k]))
        iszero(gk) && continue
        push!(pairs, FSMChannel(:up_pair, :right, (1, 0),
                                [OpenEnd(s1, gk, compose(_adag_up, _parity))],
                                [CloseEnd(s2, 1.0, _a_up)]))
        push!(pairs, FSMChannel(:up_pair_ad, :right, (-1, 0),
                                [OpenEnd(s1, gk, compose(_parity, _a_up))],
                                [CloseEnd(s2, 1.0, _adag_up)]))
        push!(pairs, FSMChannel(:down_pair, :right, (0, 1),
                                [OpenEnd(s1, gk, compose(_adag_down, _parity))],
                                [CloseEnd(s2, 1.0, _a_down)]))
        push!(pairs, FSMChannel(:down_pair_ad, :right, (0, -1),
                                [OpenEnd(s1, gk, compose(_parity, _a_down))],
                                [CloseEnd(s2, 1.0, _adag_down)]))
    end
    channels = vcat(FSMChannel[l_hop_up, l_hop_down, l_hc_up, l_hc_down,
                               r_hop_up, r_hop_down, r_hc_up, r_hc_down], pairs)
    _drop_zero_ends!.(channels)
    return ChainSpec(imp, kinds, spaces, onsite, channels)
end

# zero-coefficient ends would produce zero-valued tensor entries, which the
# SchurMPO sparse iteration drops — breaking the environment propagation of the
# channel through those sites. Drop them at construction instead; the MPO is
# unchanged (they contributed nothing) and the channel carriers stay connected.
function _drop_zero_ends!(c::FSMChannel)
    filter!(o -> !iszero(o.coef), c.opens)
    filter!(cl -> !iszero(cl.coef), c.closes)
    return c
end
