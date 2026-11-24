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
    return _merge_channels(chain_spec(x), [_ExpansionChannel(_fiducial_entries(x, imp_up, imp_down), 1.0)]; T)
end

# ==============================================================================
# Channel-entry machinery shared by the plain fiducial and the equivalent
# expansion. A channel is a (possibly multi-branch) product-chain description as
# per-site tensor entries (ql, occupation, cr, amplitude); assembling several of
# them embeds each channel into its own direct-sum bond slots — the `changebond!`
# mechanics of the GrassmannMPS stack: the fused bond spaces are allocated
# zero-filled and every entry lands on its own slot (shared charge labels pick up
# multiplicity slots in channel order), so channels never interfere and the
# represented vector is exactly the sum of the channel vectors. The chain is
# finalized right-canonical with QR/LQ (never SVD, so the zero-net-weight
# redundancy channels survive).
# ==============================================================================

"one assembly channel: per-site entries (ql, occ, cr, amp) and a boundary
amplitude folded into the first site tensor at merge time"
struct _ExpansionChannel
    entries::Vector{Vector{Tuple{Charge,Charge,Charge,Float64}}}
    amp::Float64
end

"per-cut bond sector label lists of a channel (each a dim-1 sector)"
function _channel_bonds(entries::Vector{Vector{Tuple{Charge,Charge,Charge,Float64}}})
    L = length(entries)
    bonds = Vector{Vector{Charge}}(undef, L + 1)
    bonds[1] = [(0, 0)]
    for i in 1:L
        bonds[i + 1] = sort(unique(e[3] for e in entries[i]))
    end
    return bonds
end

"sharp product-chain entries of an occupation tuple (single bond path)"
function _sharp_entries(occs::Vector{Charge})
    L = length(occs)
    entries = Vector{Vector{Tuple{Charge,Charge,Charge,Float64}}}(undef, L)
    q = (0, 0)
    for i in 1:L
        entries[i] = [(q, occs[i], _add(q, occs[i]), 1.0)]
        q = _add(q, occs[i])
    end
    return entries
end

"""
    _fiducial_entries(x, imp_up, imp_down; shift = (0,0), impamp = _ -> 1.0)

Thermofield-double fiducial as per-site channel entries: per pair the four branch
occupations τ (identically on the (a₁, a₂) sites) with thermal branch weights
√(n^|τ|(1−n)^{2−|τ|}), the impurity carrying its sharp occupations. `shift` and
`impamp` express the mechanical impurity insertion d^δ|…⟩ in the entry language:
past the impurity every bond label is shifted by `shift`, the impurity entry
carries the occupation `(imp_up, imp_down) + shift` and the amplitude
`impamp(ql)` — the operator-table amplitude times the mechanical Jordan–Wigner
string (−1)^{ql₁+ql₂} of the insertion across the left part.
"""
function _fiducial_entries(x::Union{SISBD,SIDBD}, imp_up::Int, imp_down::Int;
                           shift::Charge = (0, 0), impamp::Function = _ -> 1.0)
    spec = chain_spec(x)
    L = length(spec.kinds)
    imp = spec.imp_site
    a1n = Dict{Int,Float64}(i => n for (i, j, n) in _mode_pairs(x))   # a₁ site => n_k
    τs = ((0, 0), (1, 0), (0, 1), (1, 1))
    dql(i, q) = i > imp ? _add(q, shift) : q
    dcr(i, q) = i ≥ imp ? _add(q, shift) : q

    qs = Charge[(0, 0)]                                # running branches (unshifted)
    pairentries = Tuple{Charge,Charge,Charge,Float64}[]   # (qin, τ, m, v) seeded by the pair's a₁ site
    entries = Vector{Vector{Tuple{Charge,Charge,Charge,Float64}}}(undef, L)
    for i in 1:L
        es = Tuple{Charge,Charge,Charge,Float64}[]
        if i == imp
            p = (imp_up, imp_down)
            for q in sort(qs)
                push!(es, (q, _add(p, shift), _add(_add(q, p), shift), impamp(q)))
            end
            qs = [_add(q, p) for q in sort(qs)]
        elseif haskey(a1n, i)
            n = a1n[i]
            pairentries = Tuple{Charge,Charge,Charge,Float64}[]
            for qin in sort(qs), τ in τs
                # per-pair weight factor ONLY: multiplying running amplitudes here
                # would double-count every earlier pair's weight
                v = sqrt((1 - n)^(2 - τ[1] - τ[2]) * n^(τ[1] + τ[2]))
                iszero(v) && continue
                push!(pairentries, (qin, τ, _add(qin, τ), v))
            end
            for (qin, τ, m, v) in pairentries
                push!(es, (dql(i, qin), τ, dcr(i, m), v))
            end
            qs = sort(unique(_add(m, τ) for (qin, τ, m, _) in pairentries))
        else                                     # a₂ site: close the pair with unit blocks
            for (qin, τ, m, v) in pairentries
                push!(es, (dql(i, m), τ, dcr(i, _add(m, τ)), 1.0))
            end
            qs = sort(unique(_add(m, τ) for (qin, τ, m, _) in pairentries))
        end
        entries[i] = es
    end
    return entries
end

"""
    _merge_channels(spec, channels; T = ComplexF64) -> DMRG.MPS

Embed the channels into one MPS by direct-sum bond embedding: the bond space at
each cut is the ⊕ over the channels' label lists (shared labels pick up
multiplicity slots in channel order) and every entry is placed zero-filled at
its own (row, col) slot — the state is exactly `Σ_c amp_c·|chain_c⟩`. The
boundary amplitudes fold into the first site tensor; the chain is finalized
right-canonical (QR/LQ) with site 1 in the mixed form the sweep expects.
"""
function _merge_channels(spec::ChainSpec, channels::Vector{_ExpansionChannel};
                         T::Type{<:Number} = ComplexF64)
    L = length(spec.kinds)
    nch = length(channels)
    nch ≥ 1 || throw(ArgumentError("no channels"))
    cbonds = [_channel_bonds(ch.entries) for ch in channels]

    S = eltype(spec.spaces)
    bonds = Vector{S}(undef, L + 1)
    bonds[1] = oneunit(S)
    for i in 2:(L + 1)
        bonds[i] = foldl(TensorKit.:⊕, (_bond_space(cbonds[c][i]) for c in 1:nch))
    end

    lookup = Dict(_sector((q1, q2)) => (q1, q2) for q1 in 0:L for q2 in 0:L)
    cof(s) = haskey(lookup, s) ? lookup[s] : lookup[conj(s)]
    Ws = Vector{Any}(undef, L)
    for i in 1:L
        Bin, Bout = bonds[i], bonds[i + 1]
        W = zeros(T, Bin ⊗ spec.spaces[i] ← Bout)
        # row base per (fused key, bond ql, phys p) and column base per bond cr,
        # derived from TensorKit's own enumeration — never hand-derive the
        # canonical block layout
        rowstarts = Dict{Tuple{Charge,Charge,Charge},Int}()
        rows = Dict{Charge,Int}()
        colbase = Dict{Charge,Int}()
        cols = Dict{Charge,Int}()
        for (f1, f2) in fusiontrees(W)
            ql, p = cof(f1.uncoupled[1]), cof(f1.uncoupled[2])
            k3 = (cof(f1.coupled), ql, p)
            if !haskey(rowstarts, k3)
                rowstarts[k3] = get(rows, k3[1], 1)
                rows[k3[1]] = get(rows, k3[1], 1) + dim(Bin, _sector(ql)) * dim(spec.spaces[i], _sector(p))
            end
            cr = cof(f2.uncoupled[1])
            if !haskey(colbase, cr)
                colbase[cr] = get(cols, cr, 1)
                cols[cr] = get(cols, cr, 1) + dim(Bout, _sector(cr))
            end
        end
        # channel order = slot order within a shared charge label; the first bond
        # is trivial — all channels share its single row, columns discriminate
        prefin = Dict{Charge,Int}()
        prefout = Dict{Charge,Int}()
        for (c, ch) in enumerate(channels)
            offin = Dict{Charge,Int}()
            offout = Dict{Charge,Int}()
            for q in cbonds[c][i]
                if i == 1
                    offin[q] = 0
                else
                    offin[q] = get(prefin, q, 0)
                    prefin[q] = get(prefin, q, 0) + 1
                end
            end
            for q in cbonds[c][i + 1]
                offout[q] = get(prefout, q, 0)
                prefout[q] = get(prefout, q, 0) + 1
            end
            for (ql, occ, cr, a) in ch.entries[i]
                g = _add(ql, occ)
                g == cr || error("channel entry charge mismatch at site $i")
                row = get(rowstarts, (g, ql, occ), nothing)
                row === nothing && error("no row slot for $ql ⊗ $occ at site $i")
                col = get(colbase, cr, nothing)
                col === nothing && error("no column slot for $cr at site $i")
                block(W, _sector(g))[row + get(offin, ql, 0), col + get(offout, cr, 0)] =
                    T((i == 1 ? ch.amp : 1.0) * a)
            end
        end
        Ws[i] = W
    end
    ψ = MPS(Vector{typeof(Ws[1])}(Ws))
    return _collapse_final(_right_canonicalize!(ψ))
end

"rebuild the last tensor of an MPS with the final bond re-grouped by the
abs-decoded charge, in the CANONICAL positive-label convention (fusion-tree
labels happen to be the constructed domain labels, but `space(t, 3)` presents
their dualized view — decoding through |charge| and rebuilding with `_sector`
is the only convention-stable route, which is what the same-charge sandwich
dots require). `keep` restricts the charges, `collapse` sums each charge's
multiplicity slots into one; the move is tree-level throughout."
function _rebuild_final(ψ::DMRG.MPS; keep::Union{Nothing,Set{Charge}} = nothing,
                        collapse::Bool = false)
    WL = ψ[end]
    T = scalartype(WL)
    groups = Dict{Charge,Vector{Tuple{Any,Any}}}()
    for (f1, f2) in fusiontrees(WL)
        length(f2.uncoupled) == 1 || error("unexpected fusion tree arity")
        q = (abs(Int(f2.uncoupled[1][1].charge)), abs(Int(f2.uncoupled[1][2].charge)))
        (keep === nothing || q ∈ keep) || continue
        push!(get!(Vector{Tuple{Any,Any}}, groups, q), (f1, f2))
    end
    isempty(groups) && error("empty final bond")
    dom2 = foldl(TensorKit.:⊕,
                 (Rep[U₁×U₁](_sector(q) => (collapse ? 1 : size(WL[groups[q][1]...], ndims(WL[groups[q][1]...]))))
                  for q in sort(collect(keys(groups)))))
    W2 = zeros(T, space(WL, 1) ⊗ space(WL, 2) ← dom2)
    dst = Dict{Tuple{Any,Charge},Tuple{Any,Any}}()
    for (f1, f2) in fusiontrees(W2)
        dst[(f1, (abs(Int(f2.uncoupled[1][1].charge)), abs(Int(f2.uncoupled[1][2].charge))))] = (f1, f2)
    end
    for q in keys(groups), t in groups[q]
        b = WL[t...]
        W2[dst[(t[1], q)]...] .= collapse ? sum(b; dims = ndims(b)) : b
    end
    data = copy(getfield(ψ, :data))
    data[end] = W2
    return W2, DMRG.MPS(data, copy(getfield(ψ, :s)))
end

"collapse the final bond's per-sector multiplicity slots into one slot per sector
(summed) — the represented vector is unchanged, but the final bond then carries
the canonical positive charge labels only, which is what the same-charge dots of
the GF sandwich (and per-sector splits) require. Returns a copy; a no-op when
every sector is already one-dimensional."
function _collapse_final(ψ::DMRG.MPS)
    WL = ψ[end]
    trees = fusiontrees(WL)
    all(t -> size(WL[t...], ndims(WL[t...])) == 1, trees) && return ψ
    return _rebuild_final(ψ; collapse = true)[2]
end

"final-bond charges of an MPS, abs-decoded from the fusion trees (convention-free)"
function _final_charges(ψ::DMRG.MPS)
    return Set{Charge}((abs(Int(f2.uncoupled[1][1].charge)), abs(Int(f2.uncoupled[1][2].charge)))
                       for (f1, f2) in fusiontrees(ψ[end]))
end

"restrict the final bond to the given charges — the sectors missing from the
other chain carry truncation-level noise only, so the GF sandwich dots restrict
both chains to the common set before contracting. Always rebuilds, so both
sides of the dot end up in the same canonical label convention."
function _restrict_final(ψ::DMRG.MPS, keep::Set{Charge})
    return _rebuild_final(ψ; keep = keep)[2]
end

"split a multi-sector-final-bond MPS into definite-charge normalized components.
The last tensor is grouped at the FUSION-TREE level by its abs-decoded domain
charge (tree labels are the constructed domain labels, so the grouping is
convention-free)."
function _split_sectors(ψ::DMRG.MPS, lookup, ::Type{T}) where {T<:Number}
    comps = DMRG.MPS[]
    weights = Float64[]
    charges = Charge[]
    WL = ψ[end]
    groups = Dict{Charge,Vector{Tuple{Any,Any}}}()
    for (f1, f2) in fusiontrees(WL)
        length(f2.uncoupled) == 1 || error("unexpected fusion tree arity")
        s = f2.uncoupled[1]
        q = (abs(Int(s[1].charge)), abs(Int(s[2].charge)))
        push!(get!(Vector{Tuple{Any,Any}}, groups, q), (f1, f2))
    end
    isempty(groups) && error("empty sector split")
    for q in sort(collect(keys(groups)))
        # project the last tensor onto the single output charge (rebuild with a
        # 1-sector domain so the component is a genuine single-sector MPS); the
        # sector's full multiplicity moves along (expanded states carry several
        # channels per charge label)
        WLp, φ = _rebuild_final(ψ; keep = Set{Charge}([q]))
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
# Equivalent initial-state expansion (the OptimalExpansion route), compact form.
#
# The expanded state is ONE multi-sector MPS: the base channel is the full
# thermofield-double fiducial in its running-branch form (the per-pair branch
# superposition — never the 4^K configuration enumeration, which degenerates at
# large K where the thermal weight spreads over exponentially many
# configurations and every per-configuration slot loses its weight), and every
# far-pair seeding direction enters as a sharp ± folded pair of channels. The
# bond embedding and the right-canonical finalization follow the `changebond!`
# mechanics of the GrassmannMPS stack: bonds grow by direct-sum embedding
# (zero-filled slots placed channel by channel) and the chain is finalized
# right-canonical with QR/LQ (never SVD, so the zero-weight redundancy channels
# survive).
#
# The two-site TDVP sweep cannot seed the far-pair amplitudes (impurity ↔ bath
# mode beyond the current pair) that the first H application creates on the
# insertion chains (Kohn–Santoro 2107.02807, App. C): those directions lie
# outside the tangent space and are projected away. The cure here is an
# EQUIVALENT expansion of the chain — |χ⟩ → |χ⟩ ⊕ |φ_j⟩ ⊕ (−|φ_j⟩) over direct
# sums of bond sectors — which leaves the state exactly unchanged but makes
# every far-pair direction a tangent-space generator, so the amplitudes seed at
# O(δt·v) from the first step. With the expansion the truncation needs neither
# the KS MinD crutch nor a large bond dimension (see the dynamic truncation in
# gf_greater_lesser_t); directions without amplitude are simply truncated away.
# The base channel's charge branches are intrinsic to the exact state (its bond
# carries the reachable partial-charge labels, ~K² for K modes); the budget D
# bounds the ADDED seeding slots.
#
# Implementation notes (all pinned by dot checks):
#   * TensorKit fused blocks are keyed by the single fused charge and their rows
#     follow the canonical fusion-tree order (phys sector outer sorted × bond
#     sector inner sorted × slot), NOT the channel construction order;
#   * the ± coefficients fold into the first site tensors because the MPS
#     boundaries must carry one-dimensional spaces; they cancel pairwise, so
#     WHICH directions are admitted never changes the represented state — only
#     the seeding quality;
#   * TDVP never divides by Schmidt values, so the exactly-zero weights of the
#     redundant directions are safe.
# ==============================================================================

"per-site charge label of a particle-occupation tuple"
_site_label(spec::ChainSpec, i::Int, occ::Charge) = occ

"""
    expanded_state(spec, channels; base_vals = nothing, T = ComplexF64) -> DMRG.MPS

Channel-sum MPS over the direct-sum bond sectors: `channels` is a vector of
(per-site occupation tuple, amplitude) pairs and the state is Σ_c amp_c|occs_c⟩
exactly. Boundary coefficients fold into the first/last site tensor.
`base_vals[i]` multiplies every channel's entry at site i (insertion amplitudes).

    expanded_state(spec, occs0, phis; base_vals = nothing, T = ComplexF64)

Sharp wrapper: represent the product basis state `occs0` as the channel sum
|occs0⟩ ⊕ Σ_j (|φ_j⟩ ⊕ −|φ_j⟩). All φ charge must conserve the total.
"""
function expanded_state(spec::ChainSpec, channels::Vector{<:Tuple{Vector{Charge},<:Real}};
                        base_vals = nothing, T::Type{<:Number} = ComplexF64)
    L = length(spec.kinds)
    base_vals = base_vals === nothing ? ones(L) : base_vals
    nch = length(channels)

    labels = [[ch[1][i] for i in 1:L] for ch in channels]
    bondq = [Vector{Charge}(undef, L + 1) for _ in 1:nch]
    for c in 1:nch
        bondq[c][1] = (0, 0)
        for i in 1:L
            bondq[c][i + 1] = _add(bondq[c][i], labels[c][i])
        end
    end
    total = bondq[1][L + 1]
    all(==(total), (bondq[c][L + 1] for c in 1:nch)) ||
        throw(ArgumentError("expansion channels must conserve the total charge"))

    S = eltype(spec.spaces)
    bonds = Vector{S}(undef, L + 1)
    bonds[1] = oneunit(S)
    for i in 2:L
        bonds[i] = foldl(TensorKit.:⊕, (S(_sector(bondq[c][i]) => 1) for c in 1:nch))
    end
    bonds[L + 1] = S(_sector(total) => 1)

    Ws = Vector{Any}(undef, L)
    lookup = Dict(_sector((q1, q2)) => (q1, q2) for q1 in 0:L for q2 in 0:L)
    cof(s) = haskey(lookup, s) ? lookup[s] : lookup[conj(s)]
    for i in 1:L
        W = zeros(T, bonds[i] ⊗ spec.spaces[i] ← bonds[i + 1])
        # fusion-tree row base per (fused key, bond ql, phys p), derived from
        # TensorKit's own enumeration — never hand-derive the canonical block
        # row order (phys-outer vs bond-outer depends on the codomain split)
        rowstarts = Dict{Tuple{Charge,Charge,Charge},Int}()
        r = Dict{Charge,Int}()                # running row per fused-charge block
        for (f1, f2) in fusiontrees(W)
            k3 = (cof(f1.coupled), cof(f1.uncoupled[1]), cof(f1.uncoupled[2]))
            haskey(rowstarts, k3) || begin
                rowstarts[k3] = get(r, k3[1], 1)
                r[k3[1]] = get(r, k3[1], 1) + dim(bonds[i], _sector(k3[2])) * dim(spec.spaces[i], _sector(k3[3]))
            end
        end
        for c in 1:nch
            g = _add(bondq[c][i], labels[c][i])
            g == bondq[c][i + 1] || error("charge mismatch at site $i channel $c")
            key = _sector(g)
            ql = bondq[c][i]
            slot = i == 1 ? 0 : count(cc -> _sector(bondq[cc][i]) == _sector(ql), 1:(c - 1))
            row = rowstarts[(g, ql, labels[c][i])] + slot
            col = i == L ? 1 : 1 + count(cc -> _sector(bondq[cc][i + 1]) == key, 1:(c - 1))
            # the channel amplitude folds into the FIRST site only (per-site
            # factors would multiply to amp^L along the path)
            val = (i == 1 ? channels[c][2] : 1.0) * base_vals[i]
            block(W, key)[row, col] = T(val)
        end
        Ws[i] = W
    end
    ψ = MPS(Vector{typeof(Ws[1])}(Ws))
    return _right_canonicalize!(ψ)
end

function expanded_state(spec::ChainSpec, occs0, phis; base_vals = nothing,
                        T::Type{<:Number} = ComplexF64)
    channels = [(collect(occs0), 1.0)]
    for φ in phis
        push!(channels, (collect(φ), 1.0))
        push!(channels, (collect(φ), -1.0))
    end
    return expanded_state(spec, channels; base_vals, T)
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

# ---------------- thermofield-double initial-state expansion ----------------

"""
    _top_configs(x, imp_occs; maxn, minweight) -> Vector{Tuple{Vector{Charge},Float64}}

Best-first (descending-weight) enumeration of the thermofield-double fiducial
configurations WITHOUT the full 4^K sweep: a k-best-product walk over the
per-pair branch choices (each pair carries τ ∈ {(0,0),(1,0),(0,1),(1,1)}
identically on its (a₁, a₂) branches with weight n^|τ|(1−n)^{2−|τ|}, the
impurity carrying `imp_occs`). Stops after `maxn` configurations or when the
weight falls below `minweight` RELATIVE to the best configuration — at large K
every absolute weight is exponentially small and an absolute cutoff would
reject everything.
"""
function _top_configs(x::Union{SISBD,SIDBD}, imp_occs::Charge; maxn::Int, minweight::Float64)
    spec = chain_spec(x)
    L = length(spec.kinds)
    pairs = _mode_pairs(x)
    K = length(pairs)
    isempty(pairs) && return Tuple{Vector{Charge},Float64}[]
    choices = [[(t, nk^(t[1] + t[2]) * (1 - nk)^(2 - t[1] - t[2]))
                for t in ((0, 0), (1, 0), (0, 1), (1, 1))] for (_, _, nk) in pairs]
    foreach(v -> sort!(v; by = p -> -p[2]), choices)
    occsof(code) = begin
        o = Vector{Charge}(undef, L)
        for (k, (i1, i2, _)) in enumerate(pairs)
            t = choices[k][code[k]][1]
            o[i1] = t
            o[i2] = t
        end
        o[spec.imp_site] = imp_occs
        o
    end
    weightof(code) = prod(choices[k][code[k]][2] for k in 1:K)
    bestcode = ones(Int, K)
    out = Tuple{Vector{Charge},Float64}[]
    pending = Tuple{Float64,Vector{Int}}[(-weightof(bestcode), bestcode)]
    seen = Set{Vector{Int}}([bestcode])
    wbest = 0.0
    while !isempty(pending) && length(out) < maxn
        sort!(pending; by = e -> e[1])
        w, code = popfirst!(pending)
        if isempty(out)
            wbest = -w
            iszero(wbest) && break
        end
        (-w < minweight * wbest) && break
        push!(out, (occsof(code), -w))
        for k in 1:K, j in (code[k] + 1):length(choices[k])
            code2 = copy(code)
            code2[k] = j
            (code2 in seen) && continue
            push!(seen, code2)
            push!(pending, (-weightof(code2), code2))
        end
    end
    return out
end

"""
    _seed_directions(x, imp_occs; shifted_imp = nothing, npairs = typemax(Int))
        -> Vector{Tuple{Vector{Charge},Float64}}

Far-pair seeding directions of the thermo-field chain, in TWO tiers under the
±-pair budget `npairs`:

  * breadth — one representative per direction descriptor (σ, bath site, kind)
    over the impurity-hopping channels (`*_hop`: the impurity fills while the
    bath site empties; `*_hc`: the reverse; the intra-mode pairing channels are
    nearest-neighbor and seeded by the sweep window itself), admitted by
    coupling strength `|coef|`. Every descriptor the dynamics can reach thus
    becomes a tangent generator;
  * depth — the per-configuration φ family of the previous construction,
    admitted best-first over the thermal config weights × coupling (the
    first-order far residual is Σ_j V_j Σ_c √w_c|φ_j^c⟩, so covering the
    family restores its exact projection; the configs come from a
    k-best-product walk with a RELATIVE weight cutoff, which stays valid at
    large K where every absolute config weight is exponentially small).

Each direction carries the occupation tuple it points to: the dominant
fiducial configuration conditioned on the availability (pair-level: both the
a₁ and the a₂ site carry the branch) with the impurity and the bath site
flipped to the target; the depth tier flips inside the given config. The
availability at the sharp (possibly inserted) impurity occupation and at the
pair branch decides existence, exactly as the per-configuration filter of
`far_pair_directions` but aggregated over configurations.
"""
function _seed_directions(x::Union{SISBD,SIDBD}, imp_occs::Charge;
                          shifted_imp::Union{Nothing,Charge} = nothing,
                          npairs::Int = typemax(Int))
    spec = chain_spec(x)
    imp = spec.imp_site
    L = length(spec.kinds)
    imp_target = shifted_imp === nothing ? imp_occs : shifted_imp

    cand = Tuple{Int,Int,Bool,Float64}[]
    for ch in spec.channels
        str = String(ch.name)
        (startswith(str, "up") || startswith(str, "down")) || continue
        occursin("pair", str) && continue
        σ = startswith(str, "up") ? 1 : 2
        filling = endswith(str, "_hop")        # impurity fills + bath site empties
        for end_ in (ch.opens..., ch.closes...)
            s = end_.site
            (iszero(end_.coef) || s == imp || abs(s - imp) == 1) && continue
            push!(cand, (σ, s, filling, end_.coef))
        end
    end

    sitepair = Dict{Int,Float64}()
    for (i1, i2, n) in _mode_pairs(x)
        sitepair[i1] = n
        sitepair[i2] = n
    end
    weights(n) = [(t, n^(t[1] + t[2]) * (1 - n)^(2 - t[1] - t[2]))
                  for t in ((0, 0), (1, 0), (0, 1), (1, 1))]
    dominant(n) = first(sort!(weights(n); by = p -> -p[2]))[1]

    out = Tuple{Vector{Charge},Float64}[]
    seen = Set{Vector{Charge}}()
    # breadth tier: one representative per descriptor, coupling strength first
    for (σ, s, filling, coef) in sort!(cand; by = p -> -p[4])
        length(out) ≥ npairs && return out
        (filling ? imp_target[σ] == 0 : imp_target[σ] == 1) || continue
        n = get(sitepair, s, 0.0)
        compat = filter(p -> p[1][σ] == (filling ? 1 : 0), weights(n))
        isempty(compat) && continue
        iszero(maximum(p[2] for p in compat)) && continue
        τs = first(sort!(compat; by = p -> -p[2]))[1]
        φ = Vector{Charge}(undef, L)
        for (i1, i2, nk) in _mode_pairs(x)
            τd = dominant(nk)
            φ[i1] = τd
            φ[i2] = τd
        end
        for (i1, i2, _) in _mode_pairs(x)
            (s == i1 || s == i2) || continue
            φ[i1] = τs
            φ[i2] = τs
        end
        φ[imp] = σ == 1 ? (filling ? 1 : 0, imp_target[2]) : (imp_target[1], filling ? 1 : 0)
        φ[s] = σ == 1 ? (filling ? 0 : 1, τs[2]) : (τs[1], filling ? 0 : 1)
        push!(seen, φ)
        push!(out, (φ, coef))
    end
    # depth tier: the per-configuration φ family, config weight × coupling first
    2 * length(out) ≥ 2 * npairs && return out
    configs = _top_configs(x, imp_occs; maxn = 65536, minweight = 1.0e-14)
    for (occs, w) in configs, (σ, s, filling, coef) in cand
        length(out) ≥ npairs && return out
        (filling ? imp_target[σ] == 0 : imp_target[σ] == 1) || continue
        (filling ? occs[s][σ] == 1 : occs[s][σ] == 0) || continue
        φ = copy(occs)
        φ[imp] = σ == 1 ? (filling ? 1 : 0, imp_target[2]) : (imp_target[1], filling ? 1 : 0)
        φ[s] = σ == 1 ? (filling ? 0 : 1, occs[s][2]) : (occs[s][1], filling ? 0 : 1)
        φ in seen && continue
        push!(seen, φ)
        push!(out, (φ, coef * w))
    end
    return out
end

"""
    expanded_fiducial(x; imp_up=0, imp_down=0, D=EXPANDED_DBUDGET, T=ComplexF64) -> SectorSumState

The thermofield-double fiducial as the EQUIVALENT expansion: ONE multi-sector
MPS whose base channel is the full fiducial (the running-branch pair
superposition — never the 4^K configuration enumeration) and whose remaining
channels are the far-pair ± seeding pairs admitted under the bond budget `D`
(see `_seed_directions`, best coupling first) — the state is exactly
`fiducial_state(x)` but every far-pair direction is a tangent-space generator,
so two-site TDVP seeds the impurity ↔ bath amplitudes beyond its current pair
from the first step (KS 2107.02807 App. C).

The budget `D` bounds the ADDED seeding slots (`(D − basebond) ÷ 2` ± pairs);
the base channel's charge branches are intrinsic to the exact state. Dropped
directions only delay their amplitudes — the state stays exact regardless of
which directions are admitted (the ± pairs cancel). The component spans every
total-charge sector: `weights = [1.0]`, `charges = []`.
"""
function expanded_fiducial(x::Union{SISBD,SIDBD}; imp_up::Int = 0, imp_down::Int = 0,
                           D::Int = EXPANDED_DBUDGET, T::Type{<:Number} = ComplexF64)
    spec = chain_spec(x)
    base = _ExpansionChannel(_fiducial_entries(x, imp_up, imp_down), 1.0)
    basebond = maximum(length(_channel_bonds(base.entries)[i]) for i in 2:length(spec.kinds))
    directions = _seed_directions(x, (imp_up, imp_down);
                                  npairs = max(EXPANDED_MIN_DIR, D - basebond) ÷ 2)
    channels = _ExpansionChannel[base]
    for (φ, _) in directions
        push!(channels, _ExpansionChannel(_sharp_entries(φ), 1.0))
        push!(channels, _ExpansionChannel(_sharp_entries(φ), -1.0))
    end
    ψ = _merge_channels(spec, channels; T)
    nrm = sqrt(real(dot(ψ, ψ)))
    ψ[1] = ψ[1] * (1 / nrm)     # normalizes (a no-op up to roundoff: the ± pairs cancel)
    return SectorSumState([ψ], [1.0], Charge[])
end

"""
    _expanded_inserted(x, ψ, spin, delta; imp_up=0, imp_down=0, D=EXPANDED_DBUDGET, T) -> SectorSumState

The inserted chain d^delta_σ|ψ₀⟩ as the equivalent expansion of the (finite-
temperature) fiducial: ONE multi-sector MPS whose base channel is the
mechanically inserted fiducial entry stream (operator-table amplitude and the
per-branch Jordan–Wigner string at the impurity, every downstream bond shifted)
plus the far-pair ± seeding pairs of the inserted chain, admitted under the
bond budget `D` exactly as in `expanded_fiducial`. `nothing` when the insertion
annihilates the state.
"""
function _expanded_inserted(x::Union{SISBD,SIDBD}, ψ::SectorSumState, spin::Symbol,
                            delta::Int; imp_up::Int = 0, imp_down::Int = 0,
                            D::Int = EXPANDED_DBUDGET, T::Type{<:Number} = ComplexF64)
    spec = chain_spec(x)
    imp = spec.imp_site
    q_in = (imp_up, imp_down)
    q_out = spin === :up ? (q_in[1] + delta, q_in[2]) : (q_in[1], q_in[2] + delta)
    ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || return nothing
    tbl = spin === :up ? (delta == 1 ? _adag_up : _a_up) : (delta == 1 ? _adag_down : _a_down)
    amp = _amplitude(tbl, q_out, q_in)
    amp === nothing && return nothing
    q_d = q_out .- q_in
    jw(ql) = amp * (isodd(ql[1] + ql[2]) ? -1.0 : 1.0)
    base = _ExpansionChannel(_fiducial_entries(x, imp_up, imp_down; shift = q_d, impamp = jw), 1.0)
    basebond = maximum(length(_channel_bonds(base.entries)[i]) for i in 2:length(spec.kinds))
    directions = _seed_directions(x, (imp_up, imp_down); shifted_imp = q_out,
                                  npairs = max(EXPANDED_MIN_DIR, D - basebond) ÷ 2)
    nleft(φ) = sum(φ[t][1] + φ[t][2] for t in 1:(imp - 1); init = 0)
    channels = _ExpansionChannel[base]
    for (φ, _) in directions
        a = amp * (isodd(nleft(φ)) ? -1.0 : 1.0)
        push!(channels, _ExpansionChannel(_sharp_entries(φ), a))
        push!(channels, _ExpansionChannel(_sharp_entries(φ), -a))
    end
    χ = _merge_channels(spec, channels; T)
    nrm = sqrt(real(dot(χ, χ)))
    χ[1] = χ[1] * (1 / nrm)
    return SectorSumState([χ], copy(ψ.weights), Charge[], ψ.phase)
end
