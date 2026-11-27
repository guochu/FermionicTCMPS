# ==============================================================================
# U₁×U₁ local operator machinery for the spinful fermionic chain.
#
# The whole package works in a bosonic U₁×U₁ symmetric representation with
# explicit Jordan–Wigner bookkeeping (no fermionic grading):
#
#   * local basis |0⟩, |↑⟩, |↓⟩, |↑↓⟩ with |↑↓⟩ = c̃↑† c̃↓†|0⟩ (up-first ordering);
#     the local operators c̃ carry the intra-site JW phases;
#   * P = (−1)^n is the local parity (charge (0,0));
#   * a two-body term c_a† c_b between chain positions p_a ≠ p_b reads
#         c̃_a† (∏_{min(p) ≤ t < max(p)} P_t) c̃_b ,
#     i.e. the leftmost endpoint carries the adjacent parity — trailing for a
#     creation (c̃†P) and leading for an annihilation (Pc̃) — while the rightmost
#     endpoint is bare.
#
# These conventions were inherited from GeneralHamiltonians `sp_u1u1` and are
# pinned down by the dense G2 anchor (test/g2_dense_anchor.jl).
# ==============================================================================

# charge labels (spin-up number, spin-down number)
const Charge = Tuple{Int,Int}

_sector(q::Charge) = Irrep[U₁](q[1]) ⊠ Irrep[U₁](q[2])
_add(a::Charge, b::Charge) = (a[1] + b[1], a[2] + b[2])

"""
    particle_space()

Physical space of a chain site: sector labels are the occupations (n↑, n↓).
"""
particle_space() = Rep[U₁×U₁]((0, 0) => 1, (1, 0) => 1, (0, 1) => 1, (1, 1) => 1)

# FSM bond slots: every slot is a one-dimensional single-charge space. Closed and
# final slots carry the trivial charge, channel slots carry the charge of the
# fermionic operator that opened the channel.
slot_space(q::Charge) = Rep[U₁×U₁](_sector(q) => 1)
const SLOT_CLOSED = (0, 0)
const SLOT_FINAL = (0, 0)

# operator amplitude tables, keyed by ((q_out, q_in) => ⟨q_out|op|q_in⟩)
const OpTable = Vector{Pair{Tuple{Charge,Charge},Float64}}

table(pairs::Pair{Tuple{Charge,Charge},<:Real}...) =
    OpTable([p.first => Float64(p.second) for p in pairs])

"matrix product of two amplitude tables (t2 acts first)"
function compose(t1::OpTable, t2::OpTable)
    acc = Dict{Tuple{Charge,Charge},Float64}()
    for ((qo, qm), a1) in t1, ((qm2, qi), a2) in t2
        (qm == qm2) || continue
        key = (qo, qi)
        acc[key] = get(acc, key, 0.0) + a1 * a2
    end
    return OpTable([Pair(key, v) for (key, v) in acc if !iszero(v)])
end

scale(t::OpTable, c::Real) = OpTable([Pair(k, c * v) for (k, v) in t])

# ---------------- particle-label tables ----------------
const P_VAC = (0, 0); const P_UP = (1, 0); const P_DOWN = (0, 1); const P_FULL = (1, 1)

const _a_up = table(
    (P_VAC, P_UP) => 1.0,
    (P_DOWN, P_FULL) => 1.0,
)
const _a_down = table(
    (P_VAC, P_DOWN) => 1.0,
    (P_UP, P_FULL) => -1.0,
)
const _adag_up = table(
    (P_UP, P_VAC) => 1.0,
    (P_FULL, P_DOWN) => 1.0,
)
const _adag_down = table(
    (P_DOWN, P_VAC) => 1.0,
    (P_FULL, P_UP) => -1.0,
)
"local parity P = (−1)^n, charge (0,0)"
const _parity = table(
    (P_VAC, P_VAC) => 1.0,
    (P_UP, P_UP) => -1.0,
    (P_DOWN, P_DOWN) => -1.0,
    (P_FULL, P_FULL) => 1.0,
)
const _n_up = table((P_UP, P_UP) => 1.0, (P_FULL, P_FULL) => 1.0)
const _n_down = table((P_DOWN, P_DOWN) => 1.0, (P_FULL, P_FULL) => 1.0)

# All chain sites (impurity, a₁ and a₂ branch alike) use the particle-label
# tables below: the a₂ branch represents the f̂₂ fermions directly (the
# particle–hole labelling of the previous construction could not host the
# intra-mode pairing term of the exact Bogoliubov image under charge-pure MPO
# slots).

"""
    entry_tensor([T=Float64], ql, qr, pspace, amps)

Build an MPO-tensor entry `W[ql ⊗ pspace ← qr ⊗ pspace]` (DMRG convention: legs are
`[left-bond, phys-out | right-bond, phys-in]`) filled from amplitude pairs
`((q_out, q_in) => value)`, with charges expressed in the labels of `pspace`.
Every FSM bond slot carries a single charge sector, hence each block is 1×1.
"""
function entry_tensor(::Type{T}, ql::Charge, qr::Charge, pspace::S, amps::OpTable) where {T<:Number, S<:ElementarySpace}
    W = zeros(T, slot_space(ql) ⊗ pspace ← slot_space(qr) ⊗ pspace)
    for ((q_out, q_in), v) in amps
        g = _add(ql, q_out)
        (g == _add(q_in, qr)) || throw(ArgumentError("charge mismatch: $ql ⊗ $q_out vs $q_in ⊗ $qr"))
        copy!(block(W, _sector(g)), fill(T(v), 1, 1))
    end
    return W
end
entry_tensor(ql::Charge, qr::Charge, pspace::ElementarySpace, amps::OpTable) =
    entry_tensor(Float64, ql, qr, pspace, amps)
