# ==============================================================================
# Models: discrete single-impurity Anderson models and their thermo-field chains.
#
# Physics convention (Kohn–Santoro, arXiv:2012.01424): every bath mode k with
# frequency ω_k, spectral coupling V_k = √f_k and thermal occupation
# n_k = fermidirac(β, μ, ω_k) is doubled into
#   * an a₁ branch site (particle labels), coupled with J·V_k·√(1−n_k),
#   * an a₂ branch site (particle–hole transformed: HOLE sector labels), coupled
#     with J·V_k·√(n_k) through ordinary hopping terms (the PHT absorbs the
#     pair-creation nature of the thermofield coupling).
# The supplemented chain Hamiltonian is
#   H = εd·n_d + U·n↑n↓ + Σ_k ω_k(n_k1 + n_k2^hole)
#       + Σ_k [branch couplings (d† c + h.c.)],
# and the separable initial state |∅₁⟩⊗|F₂⟩ (a₁ vacuum, a₂ filled) reproduces the
# thermal impurity dynamics of the star Anderson model. Sites are interleaved
# (a₁, a₂) per mode; the impurity sits at the left end (single bath) or in the
# middle (double bath).
# ==============================================================================

"""Discrete fermionic bath: `discretefermionicbath` (finite β) or `discretefermionicvacuum`."""
const DiscreteFermionicBath = Union{Bath{Fermion,DiscreteSpectrum}, Vacuum{Fermion,DiscreteSpectrum}}

_inverse_temperature(b::Bath) = b.β
_inverse_temperature(::Vacuum) = Inf

"""
    thermofield_couplings(bath) -> (ws, v1, v2)

Thermo-field splitting of a discrete fermionic bath: for each mode frequency ω_k and
spectral coupling V_k = √f_k, returns the per-branch couplings
`v1_k = V_k·√(1−n_k)` (a₁ branch) and `v2_k = V_k·√(n_k)` (a₂ branch).
"""
function thermofield_couplings(bath::DiscreteFermionicBath)
    ws = frequencies(bath.f)
    vs = spectrumcouplings(bath.f)
    β = _inverse_temperature(bath)
    nk = fermidirac.(β, bath.μ, ws)
    v1 = [v * sqrt(1 - n) for (v, n) in zip(vs, nk)]
    v2 = [v * sqrt(n) for (v, n) in zip(vs, nk)]
    return ws, v1, v2
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
    side::Symbol      # :right — opens at the impurity; :left — closes at the impurity
    charge::Charge
    opens::Vector{OpenEnd}
    closes::Vector{CloseEnd}
end

const PhysSpace = typeof(particle_space())

struct ChainSpec
    imp_site::Int
    kinds::Vector{Symbol}             # :particle / :hole per chain site
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
_particle_onsite(ω::Real) = table((P_UP, P_UP) => ω, (P_DOWN, P_DOWN) => ω, (P_FULL, P_FULL) => 2ω)
_hole_onsite(ω::Real) = table((H_EMPTY, H_EMPTY) => 2ω, (H_UP, H_UP) => ω, (H_DOWN, H_DOWN) => ω)

function chain_spec(x::SISBD)
    ws, v1, v2 = thermofield_couplings(x.bath)
    K = length(ws)
    L = 2K + 1
    kinds = Vector{Symbol}(undef, L)
    spaces = Vector{PhysSpace}(undef, L)
    onsite = Vector{OpTable}(undef, L)
    kinds[1], spaces[1], onsite[1] = :particle, particle_space(), _impurity_onsite(x.εd, x.U)
    for k in 1:K
        s1, s2 = 2k, 2k + 1
        kinds[s1], spaces[s1], onsite[s1] = :particle, particle_space(), _particle_onsite(ws[k])
        kinds[s2], spaces[s2], onsite[s2] = :hole, hole_space(), _hole_onsite(ws[k])
    end

    hop_up = FSMChannel(:up_hop, :right, (1, 0), [OpenEnd(1, 1.0, compose(_adag_up, _parity))], CloseEnd[])
    hop_down = FSMChannel(:down_hop, :right, (0, 1), [OpenEnd(1, 1.0, compose(_adag_down, _parity))], CloseEnd[])
    hc_up = FSMChannel(:up_hc, :right, (-1, 0), [OpenEnd(1, 1.0, compose(_parity, _a_up))], CloseEnd[])
    hc_down = FSMChannel(:down_hc, :right, (0, -1), [OpenEnd(1, 1.0, compose(_parity, _a_down))], CloseEnd[])
    for k in 1:K
        s1, s2 = 2k, 2k + 1
        push!(hop_up.closes, CloseEnd(s1, x.J * v1[k], _a_up))
        push!(hop_up.closes, CloseEnd(s2, x.J * v2[k], _b_up))
        push!(hop_down.closes, CloseEnd(s1, x.J * v1[k], _a_down))
        push!(hop_down.closes, CloseEnd(s2, x.J * v2[k], _b_down))
        push!(hc_up.closes, CloseEnd(s1, x.J * v1[k], _adag_up))
        push!(hc_up.closes, CloseEnd(s2, x.J * v2[k], _bdag_up))
        push!(hc_down.closes, CloseEnd(s1, x.J * v1[k], _adag_down))
        push!(hc_down.closes, CloseEnd(s2, x.J * v2[k], _bdag_down))
    end
    return ChainSpec(1, kinds, spaces, onsite, FSMChannel[hop_up, hop_down, hc_up, hc_down])
end

function chain_spec(x::SIDBD)
    wsl, vl1, vl2 = thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = thermofield_couplings(x.bathr)
    Kl, Kr = length(wsl), length(wsr)
    L = 2(Kl + Kr) + 1
    imp = 2Kl + 1
    kinds = Vector{Symbol}(undef, L)
    spaces = Vector{PhysSpace}(undef, L)
    onsite = Vector{OpTable}(undef, L)
    for k in 1:Kl
        s1, s2 = 2k - 1, 2k
        kinds[s1], spaces[s1], onsite[s1] = :particle, particle_space(), _particle_onsite(wsl[k])
        kinds[s2], spaces[s2], onsite[s2] = :hole, hole_space(), _hole_onsite(wsl[k])
    end
    kinds[imp], spaces[imp], onsite[imp] = :particle, particle_space(), _impurity_onsite(x.εd, x.U)
    for k in 1:Kr
        s1, s2 = imp + 2k - 1, imp + 2k
        kinds[s1], spaces[s1], onsite[s1] = :particle, particle_space(), _particle_onsite(wsr[k])
        kinds[s2], spaces[s2], onsite[s2] = :hole, hole_space(), _hole_onsite(wsr[k])
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

    for k in 1:Kl
        s1, s2 = 2k - 1, 2k
        push!(l_hop_up.opens, OpenEnd(s1, Jl * vl1[k], compose(_parity, _a_up)))
        push!(l_hop_up.opens, OpenEnd(s2, Jl * vl2[k], compose(_hole_parity, _b_up)))
        push!(l_hop_down.opens, OpenEnd(s1, Jl * vl1[k], compose(_parity, _a_down)))
        push!(l_hop_down.opens, OpenEnd(s2, Jl * vl2[k], compose(_hole_parity, _b_down)))
        push!(l_hc_up.opens, OpenEnd(s1, Jl * vl1[k], compose(_adag_up, _parity)))
        push!(l_hc_up.opens, OpenEnd(s2, Jl * vl2[k], compose(_bdag_up, _hole_parity)))
        push!(l_hc_down.opens, OpenEnd(s1, Jl * vl1[k], compose(_adag_down, _parity)))
        push!(l_hc_down.opens, OpenEnd(s2, Jl * vl2[k], compose(_bdag_down, _hole_parity)))
    end
    for k in 1:Kr
        s1, s2 = imp + 2k - 1, imp + 2k
        push!(r_hop_up.closes, CloseEnd(s1, Jr * vr1[k], _a_up))
        push!(r_hop_up.closes, CloseEnd(s2, Jr * vr2[k], _b_up))
        push!(r_hop_down.closes, CloseEnd(s1, Jr * vr1[k], _a_down))
        push!(r_hop_down.closes, CloseEnd(s2, Jr * vr2[k], _b_down))
        push!(r_hc_up.closes, CloseEnd(s1, Jr * vr1[k], _adag_up))
        push!(r_hc_up.closes, CloseEnd(s2, Jr * vr2[k], _bdag_up))
        push!(r_hc_down.closes, CloseEnd(s1, Jr * vr1[k], _adag_down))
        push!(r_hc_down.closes, CloseEnd(s2, Jr * vr2[k], _bdag_down))
    end
    channels = FSMChannel[l_hop_up, l_hop_down, l_hc_up, l_hc_down,
                          r_hop_up, r_hop_down, r_hc_up, r_hc_down]
    return ChainSpec(imp, kinds, spaces, onsite, channels)
end
