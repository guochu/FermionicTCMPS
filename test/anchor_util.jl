using LinearAlgebra: diagm, kron, I, eigen, Hermitian, dot
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
particle basis. `kinds[i]` ∈ (:imp, :particle) and `couples` is a per-mode
`(s1, s2, ω, v1, v2, n)` tuple: weighted on-sites ω(1−n) at s1 and ω·n at s2,
impurity hoppings v1/v2, and the intra-mode pairing g = ω√(n(1−n)) between the
adjacent pair. Site ordering is the documented one (site-major, up-first within
a site).
"""
function dense_hamiltonian(kinds::Vector{Symbol}, imp::Int,
                           couples::Vector{Tuple{Int,Int,Float64,Float64,Float64,Float64}};
                           εd::Float64 = 0.0, U::Float64 = 0.0)
    L = length(kinds)
    prefixes = Vector{Any}(undef, L + 1)
    prefixes[1] = ones(1, 1)
    for s in 1:L
        prefixes[s+1] = kron(prefixes[s], _LOC.par)
    end

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
    # chain annihilation operators for the hopping/pairing terms
    imp_ops = (_kron3(prefixes[imp], _LOC.cup, _suffix_id(L - imp)),
               _kron3(prefixes[imp], _LOC.cdn, _suffix_id(L - imp)))
    for (s1, s2, ω, v1, v2, n) in couples
        # number operators must be assembled as C†C so the JW parity cancels;
        # the pairing c₁†c₂ as a raw product of JW-dressed operators carries the
        # string P_{s1} automatically (Prefix²  = 1 algebra)
        ops1 = (_kron3(prefixes[s1], _LOC.cup, _suffix_id(L - s1)),
                _kron3(prefixes[s1], _LOC.cdn, _suffix_id(L - s1)))
        ops2 = (_kron3(prefixes[s2], _LOC.cup, _suffix_id(L - s2)),
                _kron3(prefixes[s2], _LOC.cdn, _suffix_id(L - s2)))
        H .+= ω * (1 - n) .* (ops1[1]' * ops1[1] .+ ops1[2]' * ops1[2])   # a₁: ω(1−n)(n↑+n↓)
        H .+= ω * n .* (ops2[1]' * ops2[1] .+ ops2[2]' * ops2[2])         # a₂: ω·n(n↑+n↓)
        g = ω * sqrt(n * (1 - n))
        for σ in 1:2
            H .+= v1 .* (imp_ops[σ]' * ops1[σ] .+ ops1[σ]' * imp_ops[σ])
            H .+= v2 .* (imp_ops[σ]' * ops2[σ] .+ ops2[σ]' * imp_ops[σ])
            iszero(g) || (t1 = ops1[σ]' * ops2[σ]; H .+= g .* (t1 .+ t1'))
        end
    end
    return H
end

# ---------- MPO-side element extraction ----------
# particle occupation tuple (n↑, n↓) per site -> prodmps sector
function _basis_state(spec, occs; T = Float64)
    occs = collect(occs)
    sectors = [_sector(occs[k]) for k in 1:length(spec.kinds)]
    return prodmps(T, spec.spaces, sectors)
end

function _mpo_elements(x, H_dense, states)
    spec = FermionicTCMPS.chain_spec(x)
    h = mpo(x)
    cache = Dict{Any,Any}()
    state_mps(v) = get!(cache, v) do; _basis_state(spec, v); end
    # H conserves the total particle number per spin: only pairs with equal
    # conserved charges contract
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
        q1 += occs[k][1]
        q2 += occs[k][2]
    end
    return (q1, q2)
end

# ---------- dense fiducial + exact Green's functions ----------

"occupation tuple of a flat basis index (site 1 most significant, up-first)"
_occ_of(L, idx) = ntuple(s -> begin
    d = (idx - 1) ÷ 4^(L - s) % 4
    (d == 1 || d == 3 ? 1 : 0, d >= 2 ? 1 : 0)
end, L)

"""
    dense_tfd_vector(spec, modes, imp_occs; T=ComplexF64)

Dense 4ᴸ vector of the thermofield-double fiducial state: impurity occupations
`imp_occs` and, per mode `(s1, s2, n)` in chain order, the doubles
⊗_σ(√(1−n)|0 0⟩ + √n|1 1⟩) on the interleaved pair.
"""
function dense_tfd_vector(spec, modes, imp_occs; T = ComplexF64)
    L = length(spec.kinds)
    imp = spec.imp_site
    ψ = zeros(T, 4^L)
    for idx in 1:4^L
        occs = _occ_of(L, idx)
        occs[imp] == imp_occs || continue
        amp = 1.0
        for (s1, s2, n) in modes
            for σ in 1:2
                occs[s1][σ] == occs[s2][σ] || (amp = 0.0; break)
                amp *= (occs[s1][σ] == 0 ? sqrt(1 - n) : sqrt(n))
            end
            iszero(amp) && break
        end
        ψ[idx] = amp
    end
    return ψ
end

"""
    dense_gf_sandwich(H, ψ0, imp, times; imp_up=0, imp_down=0) -> (gt, lt)

Exact spin-summed impurity greater/lesser GF of the chain-frame state `ψ0` under
the dense H (the sandwich estimator of `gf_greater_lesser_t`, evaluated exactly):
G^>(t) = −i⟨ψ(t)|d_σ|χ(t)⟩ with χ = d†|ψ₀⟩, G^<(t) = i⟨η(t)|d_σ|ψ(t)⟩ with
η = d|ψ₀⟩.
"""
function dense_gf_sandwich(H, ψ0, imp, times; imp_up = 0, imp_down = 0)
    L = round(Int, log(4, size(H, 1)))
    prefixes = Vector{Any}(undef, L + 1)
    prefixes[1] = ones(1, 1)
    for s in 1:L
        prefixes[s+1] = kron(prefixes[s], _LOC.par)
    end
    F = eigen(Hermitian(H))
    g0 = F.vectors' * ψ0
    gt = zeros(ComplexF64, length(times))
    lt = zeros(ComplexF64, length(times))
    for (spin, occ) in ((:up, imp_up), (:down, imp_down))
        cup = spin === :up ? _LOC.cup : _LOC.cdn
        ann = _kron3(prefixes[imp], cup, _suffix_id(L - imp))
        χ = ann' * ψ0
        η = ann * ψ0
        χg = iszero(χ) ? nothing : F.vectors' * χ
        ηg = iszero(η) ? nothing : F.vectors' * η
        for (j, t) in enumerate(times)
            ph = exp.(-im .* F.values .* t)
            ψt = F.vectors * (ph .* g0)
            if χg !== nothing
                χt = F.vectors * (ph .* χg)
                gt[j] += -im * dot(ψt, ann * χt)
            end
            if ηg !== nothing
                ηt = F.vectors * (ph .* ηg)
                lt[j] += im * dot(ηt, ann * ψt)
            end
        end
    end
    return gt, lt
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

