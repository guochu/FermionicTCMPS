using LinearAlgebra: diagm, kron, I
using FermionicTCMPS: _sector

# ==============================================================================
# G2 dense anchor: the FSM MPO must reproduce, element by element, the dense
# Hamiltonian assembled independently by mechanical Jordan–Wigner transformation
# (plain Julia matrices, no shared code with the MPO construction path).
#
# This pins down: local operator tables (intra-site JW phases), the string rule
# c̃_a† (∏_{min ≤ t < max} P_t) c̃_b, the channel layout with coefficients at the
# non-shared end (方案 b), and the a₂-branch particle–hole convention.
# ==============================================================================

# ---------- canonical local operators (independent construction) ----------
# basis (empty, up, down, updown) with |↑↓⟩ = c↑†c↓†|0⟩ (up-first)
function _canonical_local_ops()
    cup = zeros(4, 4);   cup[1, 2] = 1.0;   cup[3, 4] = 1.0
    cdn = zeros(4, 4);   cdn[1, 3] = 1.0;   cdn[2, 4] = -1.0
    cupd = cup'; cdnd = cdn'
    par = diagm([1.0, -1.0, -1.0, 1.0])
    nup = diagm([0.0, 1.0, 0.0, 1.0])
    ndn = diagm([0.0, 0.0, 1.0, 1.0])
    nn = diagm([0.0, 0.0, 0.0, 1.0])
    # hole modes: b ≡ particle creation; hole number n_h = 1 − n per spin
    bup = copy(cupd); bdn = copy(cdnd)
    nhup = diagm([1.0, 0.0, 1.0, 0.0])
    nhdn = diagm([1.0, 1.0, 0.0, 0.0])
    return (cup = cup, cdn = cdn, par = par, nup = nup, ndn = ndn, nn = nn,
            bup = bup, bdn = bdn, nhup = nhup, nhdn = nhdn)
end

const _LOC = _canonical_local_ops()

_suffix_id(n::Int) = n == 0 ? ones(1, 1) : Matrix(1.0 * I, 4^n, 4^n)
_kron3(a, b, c) = kron(kron(a, b), c)

"""
    dense_hamiltonian(kinds, couples; εd, U) -> Matrix

Mechanical JW assembly of the thermo-field chain Hamiltonian in the 4ᴸ-dimensional
particle basis. `kinds[i]` ∈ (:imp, :particle, :hole) and `couples` is a
`Vector{Tuple{Int,Float64,Float64}}` of per bath mode `(site, v1, v2)`.
Site ordering is the documented one (site-major, up-first within a site).
"""
function dense_hamiltonian(kinds::Vector{Symbol}, imp::Int,
                           couples::Vector{Tuple{Int,Int,Float64,Float64,Float64}};
                           εd::Float64 = 0.0, U::Float64 = 0.0)
    L = length(kinds)
    prefixes = Vector{Any}(undef, L + 1)
    prefixes[1] = ones(1, 1)
    for s in 1:L
        prefixes[s+1] = kron(prefixes[s], _LOC.par)
    end
    sitename(M, s) = _kron3(prefixes[s], M, _suffix_id(L - s))

    dim = 4^L
    H = zeros(dim, dim)
    # local impurity terms — even operators carry no JW string, so embed them
    # WITHOUT the parity prefix: sitename would attach ∏_{t<imp} P_t, which flips
    # the diagonal sign on bath states with odd parity left of a mid-chain
    # impurity (every double-bath model; single-bath models have imp = 1 and an
    # empty prefix, which is how this slipped past the off-diagonal-only sampling)
    sitename_bare(M, s) = _kron3(_suffix_id(s - 1), M, _suffix_id(L - s))
    H .+= εd .* (sitename_bare(_LOC.nup, imp) .+ sitename_bare(_LOC.ndn, imp))
    (U != 0.0) && (H .+= U .* sitename_bare(_LOC.nn, imp))
    # chain annihilation operators for the hopping terms
    imp_ops = (_kron3(prefixes[imp], _LOC.cup, _suffix_id(L - imp)),
               _kron3(prefixes[imp], _LOC.cdn, _suffix_id(L - imp)))
    for (s1, s2, ω, v1, v2) in couples
        # number operators must be assembled as C†C so the JW parity cancels
        ops1 = (_kron3(prefixes[s1], _LOC.cup, _suffix_id(L - s1)),
                _kron3(prefixes[s1], _LOC.cdn, _suffix_id(L - s1)))
        ops2 = (_kron3(prefixes[s2], _LOC.bup, _suffix_id(L - s2)),
                _kron3(prefixes[s2], _LOC.bdn, _suffix_id(L - s2)))
        H .+= ω .* (ops1[1]' * ops1[1] .+ ops1[2]' * ops1[2])   # a₁: ω(n↑+n↓)
        H .+= ω .* (ops2[1]' * ops2[1] .+ ops2[2]' * ops2[2])   # a₂: ω(n_h↑+n_h↓)
        for σ in 1:2
            H .+= v1 .* (imp_ops[σ]' * ops1[σ] .+ ops1[σ]' * imp_ops[σ])
            H .+= v2 .* (imp_ops[σ]' * ops2[σ] .+ ops2[σ]' * imp_ops[σ])
        end
    end
    return H
end

# ---------- MPO-side element extraction ----------
# particle occupation tuple (n↑, n↓) per site -> prodmps sector (hole sites relabeled)
function _basis_state(spec, occs; T = Float64)
    occs = collect(occs)
    sectors = [_sector(k == spec.imp_site ? occs[k] :
                       (spec.kinds[k] === :hole ? (1 - occs[k][1], 1 - occs[k][2]) : occs[k]))
               for k in 1:length(spec.kinds)]
    return prodmps(T, spec.spaces, sectors)
end

function _mpo_elements(x, H_dense, states)
    spec = FermionicTCMPS.chain_spec(x)
    h = mpo(x)
    cache = Dict{Any,Any}()
    state_mps(v) = get!(cache, v) do; _basis_state(spec, v); end
    # H conserves (particles at particle sites + holes at hole sites) per spin:
    # only pairs with equal conserved charges contract
    err = 0.0
    for (v, w) in states
        _qtotal(spec, v) == _qtotal(spec, w) || continue
        mpo_val = expectation(state_mps(v), h, state_mps(w))
        dense_val = H_dense[_flat(spec, v), _flat(spec, w)]
        err = max(err, abs(mpo_val - dense_val))
        err = max(err, abs(conj(expectation(state_mps(w), h, state_mps(v))) - mpo_val))
    end
    return err
end

function _qtotal(spec, occs)
    q1 = q2 = 0
    for k in 1:length(occs)
        if spec.kinds[k] === :hole
            q1 += 1 - occs[k][1]
            q2 += 1 - occs[k][2]
        else
            q1 += occs[k][1]
            q2 += occs[k][2]
        end
    end
    return (q1, q2)
end

# flat index of a particle-occupation basis vector — kron convention (site 1 is the
# most significant factor), site-major with up-first in-site ordering (↑ + 2↓ per site)
_flat(spec, occs) = 1 + sum((occs[s][1] + 2 * occs[s][2]) * 4^(length(occs) - s) for s in 1:length(occs))

# enumerate all particle-occupation tuples for L sites
_all_occs(L) = [Tuple(ntuple(i -> ((v ÷ 4^(i-1)) % 4 == 0 ? (0, 0) :
                                    (v ÷ 4^(i-1)) % 4 == 1 ? (1, 0) :
                                    (v ÷ 4^(i-1)) % 4 == 2 ? (0, 1) : (1, 1)), L))
                for v in 0:(4^L - 1)]

function _anchor_model(x, H_dense; nsample::Int = 250)
    L = length(FermionicTCMPS.chain_spec(x).kinds)
    occs = _all_occs(L)
    # structured references: vacuum, filled, up-only, down-only — crossed with a sample
    refs = [occs[1], occs[end], occs[2], occs[3]]  # vacuum, full, first-excited pair
    sample = [occs[rand(1:length(occs))] for _ in 1:64]
    pairs = Any[]
    for r in refs, o in sample
        push!(pairs, (r, o))
        push!(pairs, (o, r))
    end
    # diagonal elements: random pairs essentially never repeat a state, so the
    # on-site (diagonal) part of H is only exercised here
    for o in sample, r in refs
        push!(pairs, (o, o))
    end
    for _ in 1:nsample
        push!(pairs, (rand(occs), rand(occs)))
    end
    return _mpo_elements(x, H_dense, pairs)
end

