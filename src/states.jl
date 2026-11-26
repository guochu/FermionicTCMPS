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
# Channel-entry machinery of the fiducial construction. A channel is a (possibly
# multi-branch) product-chain description as per-site tensor entries
# (ql, occupation, cr, amplitude); assembling several of them embeds each channel
# into its own direct-sum bond slots — the `changebond!` mechanics of the
# GrassmannMPS stack: the fused bond spaces are allocated zero-filled and every
# entry lands on its own slot (shared charge labels pick up multiplicity slots
# in channel order), so channels never interfere and the represented vector is
# exactly the sum of the channel vectors. The chain is finalized
# right-canonical with QR/LQ (a pure factorization — nothing is truncated).
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

"""
    _fiducial_entries(x, imp_up, imp_down)

Thermofield-double fiducial as per-site channel entries: per pair the four branch
occupations τ (identically on the (a₁, a₂) sites) with thermal branch weights
√(n^|τ|(1−n)^{2−|τ|}), the impurity carrying its sharp occupations.
"""
function _fiducial_entries(x::Union{SISBD,SIDBD}, imp_up::Int, imp_down::Int)
    spec = chain_spec(x)
    L = length(spec.kinds)
    imp = spec.imp_site
    a1n = Dict{Int,Float64}(i => n for (i, j, n) in _mode_pairs(x))   # a₁ site => n_k
    τs = ((0, 0), (1, 0), (0, 1), (1, 1))

    qs = Charge[(0, 0)]                                # running branches (unshifted)
    pairentries = Tuple{Charge,Charge,Charge,Float64}[]   # (qin, τ, m, v) seeded by the pair's a₁ site
    entries = Vector{Vector{Tuple{Charge,Charge,Charge,Float64}}}(undef, L)
    for i in 1:L
        es = Tuple{Charge,Charge,Charge,Float64}[]
        if i == imp
            p = (imp_up, imp_down)
            for q in sort(qs)
                push!(es, (q, p, _add(q, p), 1.0))
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
                push!(es, (qin, τ, m, v))
            end
            qs = sort(unique(_add(m, τ) for (qin, τ, m, _) in pairentries))
        else                                     # a₂ site: close the pair with unit blocks
            for (qin, τ, m, v) in pairentries
                push!(es, (m, τ, _add(m, τ), 1.0))
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
        # sector's full multiplicity moves along
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


# Bring sites 2..L to right-canonical form (QR-based — a pure factorization,
# nothing is truncated); site 1 carries the leftover norms, which is the mixed
# form the sweep expects as input. The bond dimensions are untouched.
function _right_canonicalize!(ψ::DMRG.MPS)
    for i in length(ψ):-1:2
        v, Q = DMRG.rightorth(ψ[i], (1,), (2, 3); alg = LQ())
        ψ[i] = permute(Q, (1, 2), (3,))
        ψ[i - 1] = ψ[i - 1] * v
    end
    return ψ
end


