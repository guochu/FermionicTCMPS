# ==============================================================================
# Initial states.
# ==============================================================================

"""
    SectorSumState

Direct sum of definite-total-charge MPS components with normalized weights:
|ψ⟩ = Σ_j phase·√w_j |ψ_j⟩ (phase shared, ±1). The thermofield-double fiducial
superposes conserved-charge sectors (each TFD pair mixes the |00⟩ and |11⟩
branches), so it cannot live in a single U₁×U₁ sector; the chain Hamiltonian
conserves the charge, hence the components evolve independently and every
observable is a weighted sum. `charges[j]` is the total charge of component j.
"""
struct SectorSumState
    comps::Vector{DMRG.MPS}
    weights::Vector{Float64}
    charges::Vector{Charge}
    phase::Float64
end
SectorSumState(comps, weights, charges) = SectorSumState(comps, weights, charges, 1.0)

"normalized sector-sum expectation value: Σ_j w_j ⟨op⟩_j"
function expectationvalue(op, s::SectorSumState)
    return sum(w * expectationvalue(op, c) for (w, c) in zip(s.weights, s.comps))
end

# different sectors are orthogonal; shared phase aside, the norm² is Σ weights
function dot(s::SectorSumState, t::SectorSumState)
    s.charges == t.charges || return zero(sum(s.weights))
    return s.phase * t.phase * sum(min.(s.weights, t.weights))
end

function evolve!(s::SectorSumState, h, tspan::Tuple{<:Real,<:Real}; δt::Real, kwargs...)
    foreach(c -> evolve!(c, h, tspan; δt, kwargs...), s.comps)
    return s
end

function evolve!(s::SectorSumState, h; δt::Real, nsteps::Integer, kwargs...)
    foreach(c -> evolve!(c, h; δt = δt, nsteps = nsteps, kwargs...), s.comps)
    return s
end

"per-component TDVP caches; `sweep!` advances every component in lockstep"
struct SectorSumCaches
    caches::Vector{Any}
end
function environments(h, s::SectorSumState)
    return SectorSumCaches([environments(h, c) for c in s.comps])
end
function sweep!(c::SectorSumCaches, alg)
    foreach(sub -> sweep!(sub, alg), c.caches)
    return c
end

"""
    fiducial_state(x::Union{SISBD,SIDBD}; imp_up=0, imp_down=0, T=ComplexF64) -> SectorSumState

Fiducial state of the thermo-field chain: the impurity orbital carries the given
up/down occupations and every bath mode k is in its thermofield-double state
⊗_σ (√(1−n_k)|0ₐ₁ 0ₐ₂⟩ + √n_k|1ₐ₁ 1ₐ₂⟩) on the interleaved (a₁, a₂) pair — the
exact pure-state image of the thermal bath (see the models.jl header). For
n_k = 0 this reduces to the (single-component) product vacuum.
"""
function fiducial_state(x::Union{SISBD,SIDBD}; imp_up::Int = 0, imp_down::Int = 0,
                        T::Type{<:Number} = ComplexF64)
    L = length(chain_spec(x).kinds)
    lookup = Dict(_sector((q1, q2)) => (q1, q2) for q1 in 0:(2L) for q2 in 0:(2L))
    ψ = _fiducial_mps(x; imp_up, imp_down, T)
    return _split_sectors(ψ, lookup, T)
end

"raw thermofield-double fiducial: ONE MPS with a multi-sector final bond"
function _fiducial_mps(x::Union{SISBD,SIDBD}; imp_up::Int, imp_down::Int,
                       T::Type{<:Number} = ComplexF64)
    (imp_up in (0, 1) && imp_down in (0, 1)) || throw(ArgumentError("impurity occupations must be 0 or 1"))
    spec = chain_spec(x)
    L = length(spec.kinds)
    imp = spec.imp_site
    a1n = Dict{Int,Float64}(i => n for (i, j, n) in _mode_pairs(x))   # a₁ site => n_k

    # running cut: branch charges (each a dim-1 sector) and their amplitudes
    lookup = Dict(_sector((q1, q2)) => (q1, q2) for q1 in 0:(2L) for q2 in 0:(2L))
    qs::Vector{Charge} = Charge[(0, 0)]
    amps = Dict{Charge,Float64}((0, 0) => 1.0)
    # entries seeded by the current pair's a₁ site: (pre-pair charge, τ, intermediate, amp)
    pairentries = Tuple{Charge,Charge,Charge,Float64}[]
    Ws = Vector{Any}(undef, L)
    for i in 1:L
        ps = spec.spaces[i]
        sortqs = sort(qs)
        if i == imp
            p = (imp_up, imp_down)
            shifted = [_add(q, p) for q in sortqs]     # mid-chain: all branches shift
            W = zeros(T, _bond_space(sortqs) ⊗ ps ← _bond_space(shifted))
            rowmap = _block_rowmap(W, lookup)
            for q in sortqs
                _set_mps_block!(W, q, p, amps[q], rowmap)
            end
            Ws[i] = W
            qs = shifted
            amps = Dict{Charge,Float64}(q => amps[qq] for (q, qq) in zip(shifted, sortqs))
        elseif haskey(a1n, i)
            n = a1n[i]
            inter = Dict{Charge,Float64}()     # intermediate branch => summed amplitude
            pairentries = Tuple{Charge,Charge,Charge,Float64}[]
            for qin in sortqs, τ in ((0, 0), (1, 0), (0, 1), (1, 1))
                v = amps[qin] * sqrt((1 - n)^(2 - τ[1] - τ[2]) * n^(τ[1] + τ[2]))
                iszero(v) && continue
                m = _add(qin, τ)
                inter[m] = get(inter, m, 0.0) + v
                push!(pairentries, (qin, τ, m, v))
            end
            outqs = sort(collect(keys(inter)))
            W = zeros(T, _bond_space(sortqs) ⊗ ps ← _bond_space(outqs))
            rowmap = _block_rowmap(W, lookup)
            for (qin, τ, m, v) in pairentries
                _set_mps_block!(W, qin, τ, v, rowmap)
            end
            Ws[i] = W
            qs = outqs
            amps = inter
        else                                     # a₂ site: close the pair with unit blocks
            newamps = Dict{Charge,Float64}()
            for (qin, τ, m, _) in pairentries
                o = _add(m, τ)
                newamps[o] = get(newamps, o, 0.0) + amps[m]
            end
            outqs = sort(collect(keys(newamps)))
            W = zeros(T, _bond_space(sortqs) ⊗ ps ← _bond_space(outqs))
            rowmap = _block_rowmap(W, lookup)
            for (qin, τ, m, _) in pairentries
                _set_mps_block!(W, m, τ, 1.0, rowmap)
            end
            Ws[i] = W
            qs = outqs
            amps = newamps
            pairentries = Tuple{Charge,Charge,Charge,Float64}[]
        end
    end
    return MPS(Vector{typeof(Ws[1])}(Ws))
end

"split a multi-sector-final-bond MPS into definite-charge normalized components"
function _split_sectors(ψ::DMRG.MPS, lookup, ::Type{T}) where {T<:Number}
    comps = DMRG.MPS[]
    weights = Float64[]
    charges = Charge[]
    for s in sort(collect(sectors(space(ψ[length(ψ)], 3))))
        # the final bond is stored dualized: convert the sector label back to the charge
        q = lookup[conj(s)]
        φ = DMRG.MPS(deepcopy(ψ.data), deepcopy(getfield(ψ, :s)))
        # project the last tensor onto the single output sector s (rebuild with a
        # 1-sector domain so the component is a genuine single-sector MPS)
        WL = φ[end]
        WLp = zeros(T, space(WL, 1) ⊗ space(WL, 2) ← Rep[U₁×U₁](conj(s) => 1))
        copyto!(block(WLp, conj(s)), block(WL, conj(s)))
        φ.data[end] = WLp
        w = real(dot(φ, φ))
        w > 1e-14 || continue
        comp = _widen_svectors(_right_canonicalize!(φ))
        comp[1] = comp[1] * (1 / sqrt(w))       # components must be normalized:
        push!(comps, comp)                       # |ψ⟩ = Σ_j phase·√w_j |ψ_j⟩
        push!(weights, w)
        push!(charges, q)
    end
    isempty(comps) && error("empty sector split")
    total = sum(weights)
    return SectorSumState(comps, weights ./ total, charges)
end

"dim-1-per-charge bond space"
_bond_space(qs::Vector{Charge}) =
    isempty(qs) ? oneunit(Rep[U₁×U₁]) :
    foldl(⊕, (Rep[U₁×U₁](_sector(q) => 1) for q in qs))

# per-tensor row lookup: fused charge => ((bond sector, phys sector) => row),
# derived from TensorKit's own `fusiontrees` enumeration (the storage order of the
# codomain decomposition) — no hand-counted canonical ordering anywhere.
function _block_rowmap(W, lookup)
    rows = Dict{Charge,Dict{Tuple{Charge,Charge},Int}}()
    for ft in fusiontrees(W)
        ctree = ft[1]
        length(ctree.uncoupled) == 2 || error("unexpected fusion tree arity")
        ql = lookup[ctree.uncoupled[1]]
        p = lookup[ctree.uncoupled[2]]
        key = lookup[ctree.coupled]
        d = get!(Dict{Tuple{Charge,Charge},Int}, rows, key)
        d[(ql, p)] = length(d) + 1
    end
    return rows
end

function _set_mps_block!(W, ql::Charge, p::Charge, val, rowmap)
    key = _add(ql, p)
    d = get(rowmap, key, nothing)
    d === nothing && error("no block for charge $key")
    row = get(d, (ql, p), nothing)
    row === nothing && error("no slot for $ql ⊗ $p in block $key")
    block(W, _sector(key))[row, 1] = val
    return W
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
# NOTE: only valid on VACUUM baths (n_k = 0), where the fiducial state is the
# sharp product vacuum. For finite-temperature baths the fiducial is the
# thermofield-double superposition and this machinery does not apply.
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

"per-site charge label of a particle-occupation tuple"
_site_label(spec::ChainSpec, i::Int, occ::Charge) = occ

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
    ψ = MPS(Vector{typeof(Ws[1])}(Ws))
    return _right_canonicalize!(ψ)
end

# Bring sites 2..L to right-canonical form (QR-based, so the exactly-zero-weight
# redundancy channels survive — an SVD truncation would drop them); site 1
# carries the leftover norms, which is the mixed form the sweep expects as
# input. The bond dimensions are untouched.
function _right_canonicalize!(ψ::DMRG.MPS)
    for i in length(ψ):-1:2
        v, Q = DMRG.rightorth(ψ[i], (1,), (2, 3); alg = LQ())
        ψ[i] = permute(Q, (1, 2), (3,))
        ψ[i - 1] = ψ[i - 1] * v
    end
    return ψ
end

"""
    far_pair_directions(spec::ChainSpec) -> occs0 -> Vector{Vector{Charge}}

Structural far-pair directions of an insertion chain: for every
impurity-annihilating channel (`*_hc`, the only half that acts on the fiducial
bath) and every non-impurity channel end with a nonzero coefficient NON-adjacent
to the impurity, the pair direction the hopping populates in one application —
impurity spin emptied and the bath site filled, subject to the occupation
availability. Adjacent pairs are omitted (the sweep window seeds them by
itself). Verified element-for-element against the dense residual support in
scripts/diag_expansion.jl.
"""
function far_pair_directions(spec::ChainSpec)
    imp = spec.imp_site
    return occs0 -> begin
        out = Vector{Vector{Charge}}()
        for ch in spec.channels
            str = String(ch.name)
            (startswith(str, "up") || startswith(str, "down")) || continue
            σ = startswith(str, "up") ? 1 : 2
            filling = endswith(str, "_hop")       # imp filled + bath emptied
            for end_ in (ch.opens..., ch.closes...)
                s = end_.site
                (iszero(end_.coef) || s == imp || abs(s - imp) == 1) && continue
                if filling
                    occs0[imp][σ] == 0 || continue
                    occs0[s][σ] == 1 || continue
                else
                    occs0[imp][σ] == 1 || continue
                    occs0[s][σ] == 0 || continue
                end
                φ = collect(occs0)
                φ[imp] = σ == 1 ? (filling ? 1 : 0, occs0[imp][2]) :
                         (occs0[imp][1], filling ? 1 : 0)
                φ[s] = σ == 1 ? (filling ? 0 : 1, occs0[s][2]) :
                       (occs0[s][1], filling ? 0 : 1)
                φ in out || push!(out, φ)
            end
        end
        out
    end
end

"sharp product occupations of the saturated (T = 0) fiducial: the a₂ branch
carries the saturated f̂₂ occupations (n_k per spin), a₁ the vacuum"
function _sharp_occs0(x::Union{SISBD,SIDBD}, imp_occs::Charge)
    spec = chain_spec(x)
    return ntuple(i -> i == spec.imp_site ? imp_occs :
                  begin
                      n = get(Dict{Int,Float64}(j => n for (i2, j, n) in _mode_pairs(x)), i, 0.0)
                      (round(Int, n), round(Int, n))
                  end, length(spec.kinds))
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
    q_out = spin === :up ? (q_in[1] + delta, imp_down) : (imp_up, q_in[2] + delta)
    ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || return nothing
    tbl = if spin === :up
        delta == 1 ? _adag_up : _a_up
    else
        delta == 1 ? _adag_down : _a_down
    end
    amp = _amplitude(tbl, q_out, q_in)
    amp === nothing && return nothing
    occs0 = _sharp_occs0(x, q_out)
    base_vals = ones(length(spec.kinds))
    base_vals[spec.imp_site] = amp
    return expanded_state(spec, occs0, phis; base_vals, T)
end
