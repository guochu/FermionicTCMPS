# ==============================================================================
# FSM MPO construction.
#
# Per chain site the MPO tensor is a Schur (upper triangular) operator matrix over
# the bond slots [closed, channel 1 …, channel n, final]:
#
#       1   openings   on-site
#       0   P-lifts    closings
#       0   0          1
#
# Coefficients are always placed at the NON-SHARED end (the impurity is shared by
# all hopping terms): right-bath channels carry them in the closing entries,
# left-bath channels in the opening entries. Every bond slot is a one-dimensional
# single-charge space (channels carry the charge of their opening operator), so
# the bond dimension is constant along the chain (vacuum + channels + final = 6
# for single and double bath alike). The Schur tensors are assembled with fully
# explicit slot spaces — no inference and no `isid` scalarization, so charged
# slots keep their charges uniformly along the chain.
# ==============================================================================

"""
    mpo(x::Union{SISBD,SIDBD}) -> DMRG.MPOHamiltonian

Finite-state-machine MPO of the thermo-field chain Hamiltonian, constructed
channel-wise and stored as `DMRG.SchurMPOTensor`s.
"""
mpo(x::Union{SISBD,SIDBD}) = _fsm_mpo(chain_spec(x), chain_spec(x).channels; onsite = true)

_split(channels::Vector{FSMChannel}) =
    (filter(c -> c.side === :left, channels), filter(c -> c.side === :right, channels))

# a site's rows carry channels reaching leftward of it, its columns channels
# reaching rightward; single-bath models have only :right channels, so their
# impurity rows are padded with the :right family (the extra slots stay dead).
function _row_channels(spec::ChainSpec, channels::Vector{FSMChannel}, k::Int)
    left, right = _split(channels)
    (k > spec.imp_site || isempty(left)) && return right
    return left
end
function _col_channels(spec::ChainSpec, channels::Vector{FSMChannel}, k::Int)
    left, right = _split(channels)
    (k < spec.imp_site && !isempty(left)) && return left
    return right
end

"""
    _schur_tensor([T=Float64], mat, row_charges, col_charges, pspace)

Assemble a `DMRG.SchurMPOTensor` from a `Matrix{Any}` of scalars/tensors plus the
explicit per-slot charge lists.
"""
function _schur_tensor(::Type{T}, mat::Matrix{Any},
                       row_charges::Vector{Charge}, col_charges::Vector{Charge},
                       pspace::S) where {T<:Number, S<:ElementarySpace}
    M = mpotensortype(S, T)
    nr, nc = size(mat)
    (length(row_charges) == nr && length(col_charges) == nc) ||
        throw(DimensionMismatch("slot charge lists must match the matrix size"))
    Os = Matrix{Union{M,T}}(undef, nr, nc)
    for i in 1:nr, j in 1:nc
        e = mat[i, j]
        if e isa AbstractTensorMap
            Os[i, j] = convert(M, e)
        elseif e isa Number
            Os[i, j] = convert(T, e)
        else
            throw(ArgumentError("entries must be tensors or numbers"))
        end
    end
    leftspaces = S[slot_space(q) for q in row_charges]
    rightspaces = S[slot_space(q) for q in col_charges]
    return SchurMPOTensor{S,M,T}(Os, leftspaces, rightspaces, pspace)
end

function _fsm_mpo(spec::ChainSpec, channels::Vector{FSMChannel}; onsite::Bool, T::Type = Float64)
    L = length(spec.kinds)
    boundary = _boundary_channels(channels)
    ws = Vector{Any}(undef, L)
    for k in 1:L
        row_ch = (k == 1) ? boundary : _row_channels(spec, channels, k)
        col_ch = (k == L) ? boundary : _col_channels(spec, channels, k)
        ws[k] = _site_tensor(spec, k, row_ch, col_ch, onsite, T)
    end
    return SchurMPO([w for w in ws])
end

# MPOHamiltonian requires the first tensor's left slots to equal the last tensor's
# right slots; double-bath models get the union of both families there (the foreign
# slots stay dead).
function _boundary_channels(channels::Vector{FSMChannel})
    left, right = _split(channels)
    (isempty(left) || isempty(right)) && return isempty(left) ? right : left
    return vcat(left, right)
end

"propagate a channel through site k iff it opened strictly left of k and can still close to its right"
_should_propagate(c::FSMChannel, k::Int) =
    (c.opens[1].site < k) && any(cl.site > k for cl in c.closes)

function _site_tensor(spec::ChainSpec, k::Int, row_ch::Vector{FSMChannel}, col_ch::Vector{FSMChannel},
                      onsite::Bool, ::Type{T}) where {T<:Number}
    nr, nc = length(row_ch) + 2, length(col_ch) + 2
    mat = Matrix{Any}(undef, nr, nc)
    fill!(mat, zero(T))
    mat[1, 1] = one(T)
    mat[nr, nc] = one(T)
    ps = spec.spaces[k]
    par = _parity
    if onsite && !isempty(spec.onsite[k])
        mat[1, nc] = entry_tensor(T, SLOT_CLOSED, SLOT_FINAL, ps, spec.onsite[k])
    end
    for (j, c) in enumerate(col_ch)
        cj = j + 1
        for o in c.opens
            (o.site == k) && (mat[1, cj] = entry_tensor(T, SLOT_CLOSED, c.charge, ps, scale(o.amps, o.coef)))
        end
        _should_propagate(c, k) && (mat[cj, cj] = entry_tensor(T, c.charge, c.charge, ps, par))
    end
    for (j, c) in enumerate(row_ch)
        rj = j + 1
        for cl in c.closes
            (cl.site == k) && (mat[rj, nc] = entry_tensor(T, c.charge, SLOT_FINAL, ps, scale(cl.amps, cl.coef)))
        end
    end
    row_charges = vcat([SLOT_CLOSED], [c.charge for c in row_ch], [SLOT_FINAL])
    col_charges = vcat([SLOT_CLOSED], [c.charge for c in col_ch], [SLOT_FINAL])
    return _schur_tensor(T, mat, row_charges, col_charges, ps)
end

const _id_table = table(
    (P_VAC, P_VAC) => 1.0,
    (P_UP, P_UP) => 1.0,
    (P_DOWN, P_DOWN) => 1.0,
    (P_FULL, P_FULL) => 1.0,
)

"""
    site_operator_mpo([T=Float64], spaces, position, amps)

Identity MPO with the single-site operator `amps` inserted at `position`. The two
slots are [closed (charge 0), inserted (charge = q_out − q_in of the operator)], so
charge-changing insertions (e.g. Green's function operators) contract cleanly
between states whose total charges differ by exactly that operator charge.
"""
function site_operator_mpo(::Type{T}, spaces::Vector{<:ElementarySpace}, position::Int, amps::OpTable) where {T<:Number}
    deltas = Charge[_add(p.first[1], (-p.first[2][1], -p.first[2][2])) for p in amps]
    isempty(deltas) && throw(ArgumentError("empty operator table"))
    all(==(deltas[1]), deltas) || throw(ArgumentError("operator table is not charge-pure"))
    q_op = deltas[1]

    ws = Vector{Any}(undef, length(spaces))
    for k in eachindex(spaces)
        mat = Matrix{Any}(undef, 2, 2)
        fill!(mat, zero(T))
        mat[1, 1] = one(T)
        mat[2, 2] = entry_tensor(T, q_op, q_op, spaces[k], _id_table)
        (k == position) && (mat[1, 2] = entry_tensor(T, SLOT_CLOSED, q_op, spaces[k], amps))
        ws[k] = _schur_tensor(T, mat, [SLOT_CLOSED, q_op], [SLOT_CLOSED, q_op], spaces[k])
    end
    return SchurMPO([w for w in ws])
end
site_operator_mpo(spaces::Vector{<:ElementarySpace}, position::Int, amps::OpTable) =
    site_operator_mpo(Float64, spaces, position, amps)
