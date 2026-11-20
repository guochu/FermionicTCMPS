# ==============================================================================
# Initial states.
# ==============================================================================

"""
    separable_state(x::Union{SISBD,SIDBD}; imp_up=0, imp_down=0, T=ComplexF64)

Separable product state of the thermo-field chain: the impurity orbital carries the
given up/down occupations and every bath branch is in its fiducial state
(a₁ sites empty, a₂ sites fully occupied = zero holes), i.e. |∅₁⟩⊗|F₂⟩.
"""
function separable_state(x::Union{SISBD,SIDBD}; imp_up::Int = 0, imp_down::Int = 0,
                         T::Type{<:Number} = ComplexF64)
    (imp_up in (0, 1) && imp_down in (0, 1)) || throw(ArgumentError("impurity occupations must be 0 or 1"))
    spec = chain_spec(x)
    sectors = [_sector(k == spec.imp_site ? (imp_up, imp_down) : (0, 0)) for k in 1:length(spec.kinds)]
    return _widen_svectors(prodmps(T, spec.spaces, sectors))
end

# DMRG's TDVP stores plain (1,1) tensors into `mps.s`, while `prodmps` leaves the
# container typed to DiagonalTensorMap; widen it so the assignment type-checks.
function _widen_svectors(ψ::DMRG.MPS)
    B = DMRG.MPSBondTensor{spacetype(ψ[1])}
    raw = getfield(ψ, :s)
    sv = Vector{Union{Missing, B}}(missing, length(raw))
    for i in eachindex(raw)
        raw[i] isa B && (sv[i] = raw[i])
    end
    return DMRG.MPS(ψ.data, sv)
end

# ==============================================================================
# Equivalent initial-state expansion (the OptimalExpansion route).
#
# The two-site TDVP sweep cannot seed the far-pair amplitudes (impurity ↔ bath
# mode beyond the current pair) that the first H application creates on the
# sharp insertion chains (Kohn–Santoro 2107.02807, App. C): those directions lie
# outside the tangent space and are projected away. The cure here is an
# EQUIVALENT expansion of the chain — |χ⟩ → |χ⟩ ⊕ |φ_j⟩ ⊕ (−|φ_j⟩) over direct
# sums of bond sectors — which leaves the state exactly unchanged but makes
# every far-pair direction a tangent-space generator, so the amplitudes seed at
# O(δt·v) from the first step. With the expansion the truncation needs neither
# the KS MinD crutch nor a large bond dimension (see TRUNC_EXPANDED in
# timeevo.jl); directions without amplitude are simply truncated away.
#
# Implementation notes (all pinned by scripts/diag_expansion.jl dot checks):
#   * TensorKit fused blocks are keyed by the single fused charge and their rows
#     follow the canonical fusion-tree order (phys sector outer sorted × bond
#     sector inner sorted × slot), NOT the channel construction order;
#   * the ± coefficients fold into the first/last site tensors because the MPS
#     boundaries must carry one-dimensional spaces;
#   * TDVP never divides by Schmidt values, so the exactly-zero weights of the
#     redundant directions are safe.
# ==============================================================================

"per-site charge label of a particle-occupation tuple (hole sites relabeled)"
_site_label(spec::ChainSpec, i::Int, occ::Charge) =
    spec.kinds[i] === :hole ? (1 - occ[1], 1 - occ[2]) : (occ[1], occ[2])

"""
    expanded_state(spec, occs0, phis; base_vals = nothing, T = ComplexF64) -> DMRG.MPS

Represent the product basis state `occs0` (particle occupations) exactly as the
channel sum |occs0⟩ ⊕ Σ_j (|φ_j⟩ ⊕ −|φ_j⟩): bond spaces are direct sums over
channels (ψ first, then per φ a ± pair) and site tensors carry one block entry
per channel in the canonical fusion-tree row order. Boundary coefficients fold
into the first/last site tensor. `base_vals[i]` overrides the ψ entry at site i
(the insertion amplitude on the impurity).
"""
function expanded_state(spec::ChainSpec, occs0, phis; base_vals = nothing, T::Type{<:Number} = ComplexF64)
    L = length(spec.kinds)
    base_vals = base_vals === nothing ? ones(L) : base_vals
    channels = [(collect(occs0), 1.0)]
    for φ in phis
        push!(channels, (collect(φ), 1.0))
        push!(channels, (collect(φ), -1.0))
    end
    nch = length(channels)

    labels = [[_site_label(spec, i, ch[1][i]) for i in 1:L] for ch in channels]
    bondq = [Vector{Charge}(undef, L + 1) for _ in 1:nch]
    for c in 1:nch
        bondq[c][1] = (0, 0)
        for i in 1:L
            bondq[c][i + 1] = _add(bondq[c][i], labels[c][i])
        end
    end
    total = bondq[1][L + 1]
    all(==(total), (bondq[c][L + 1] for c in 1:nch)) ||
        throw(ArgumentError("expansion directions must conserve the total charge"))

    S = eltype(spec.spaces)
    bonds = Vector{S}(undef, L + 1)
    bonds[1] = oneunit(S)
    for i in 2:L
        bonds[i] = foldl(TensorKit.:⊕, (S(_sector(bondq[c][i]) => 1) for c in 1:nch))
    end
    bonds[L + 1] = S(_sector(total) => 1)

    Ws = Vector{Any}(undef, L)
    lookup = Dict(_sector((q1, q2)) => (q1, q2) for q1 in 0:L for q2 in 0:L)
    for i in 1:L
        W = zeros(T, bonds[i] ⊗ spec.spaces[i] ← bonds[i + 1])
        Psectors = sort(collect(sectors(spec.spaces[i])))
        Vsectors = sort(collect(sectors(bonds[i])))
        for c in 1:nch
            g = _add(bondq[c][i], labels[c][i])
            g == bondq[c][i + 1] || error("charge mismatch at site $i channel $c")
            key = _sector(g)
            psec = _sector(labels[c][i])
            qlsec = _sector(bondq[c][i])
            row = 1
            for p2 in Psectors
                p2 < psec || continue
                for a2 in Vsectors
                    _sector(_add(lookup[a2], lookup[p2])) == key && (row += dim(bonds[i], a2))
                end
            end
            for a2 in Vsectors
                (a2 < qlsec && _sector(_add(lookup[a2], lookup[psec])) == key) &&
                    (row += dim(bonds[i], a2))
            end
            row += i == 1 ? 0 : count(cc -> _sector(bondq[cc][i]) == qlsec, 1:(c - 1))
            col = i == L ? 1 : 1 + count(cc -> _sector(bondq[cc][i + 1]) == key, 1:(c - 1))
            val = channels[c][2] * (c == 1 ? base_vals[i] : 1.0)
            block(W, key)[row, col] = T(val)
        end
        Ws[i] = W
    end
    return MPS(Vector{typeof(Ws[1])}(Ws))
end

"""
    far_pair_directions(spec::ChainSpec) -> occs0 -> Vector{Vector{Charge}}

Structural far-pair directions of an insertion chain: for every
impurity-annihilating channel (`*_hc`, the only half that acts on the fiducial
bath) and every non-impurity channel end with a nonzero coefficient NON-adjacent
to the impurity, the pair direction the hopping populates in one application —
impurity spin emptied and the bath site changed per branch kind (particle site
filled, hole site particle-removed), subject to the occupation availability.
Adjacent pairs are omitted (the sweep window seeds them by itself). Verified
element-for-element against the dense residual support in
scripts/diag_expansion.jl.
"""
function far_pair_directions(spec::ChainSpec)
    imp = spec.imp_site
    return occs0 -> begin
        out = Vector{Vector{Charge}}()
        for ch in spec.channels
            str = String(ch.name)
            (endswith(str, "_hc") && (startswith(str, "up") || startswith(str, "down"))) || continue
            σ = startswith(str, "up") ? 1 : 2
            for end_ in (ch.opens..., ch.closes...)
                s = end_.site
                (iszero(end_.coef) || s == imp || abs(s - imp) == 1) && continue
                occs0[imp][σ] == 1 || continue
                φ = collect(occs0)
                φ[imp] = σ == 1 ? (0, occs0[imp][2]) : (occs0[imp][1], 0)
                if spec.kinds[s] === :hole
                    occs0[s][σ] == 1 || continue
                    φ[s] = σ == 1 ? (0, occs0[s][2]) : (occs0[s][1], 0)
                else
                    occs0[s][σ] == 0 || continue
                    φ[s] = σ == 1 ? (1, occs0[s][2]) : (occs0[s][1], 1)
                end
                φ in out || push!(out, φ)
            end
        end
        out
    end
end

"""
    expanded_shifted_state(x, imp_up, imp_down, spin, delta, phis; T = ComplexF64) -> MPS

The insertion chain d^delta_σ|imp⟩ as an equivalently expanded product MPS
(see `expanded_state`); `nothing` when the insertion kills the state. `phis`
may be `far_pair_directions(chain_spec(x))(occs0)` or anything charge-conserving.
"""
function expanded_shifted_state(x::Union{SISBD,SIDBD}, imp_up::Int, imp_down::Int,
                                spin::Symbol, delta::Int, phis; T::Type{<:Number} = ComplexF64)
    spec = chain_spec(x)
    q_in = (imp_up, imp_down)
    q_out = spin === :up ? (q_in[1] + delta, q_in[2]) : (q_in[1], q_in[2] + delta)
    ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || return nothing
    tbl = if spin === :up
        delta == 1 ? _adag_up : _a_up
    else
        delta == 1 ? _adag_down : _a_down
    end
    amp = _amplitude(tbl, q_out, q_in)
    amp === nothing && return nothing
    occs0 = ntuple(i -> i == spec.imp_site ? q_out :
                   (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), length(spec.kinds))
    base_vals = ones(length(spec.kinds))
    base_vals[spec.imp_site] = amp
    return expanded_state(spec, occs0, phis; base_vals, T)
end
