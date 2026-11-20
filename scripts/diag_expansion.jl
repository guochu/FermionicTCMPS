# ==============================================================================
# Diagnostic: equivalent (channel-doubled) expansion of the sharp initial state
# against the vacuum-class TDVP2 seeding defect.
#
# Mechanism under test (Kohn–Santoro 2107.02807 App. C): the pipeline evolves the
# insertion chains χ = d†|ψ₀⟩ / η = d|ψ₀⟩; for vacuum baths the first H action
# creates pair excitations (impurity ↔ far bath mode) that lie OUTSIDE the
# two-site TDVP tangent space of the sharp product state, so their amplitudes
# never seed (sweep errors 0.3–1.0 at t = 1 for configs C/D/F/E/G).
#
# Proposed cure: an EQUIVALENT expansion of the insertion chain —
#     |χ⟩ → |χ⟩ ⊕ |φ⟩ ⊕ (−|φ⟩)
# with φ = the residual (far-pair) directions of Hχ. The channels carry φ with
# opposite signs so the total state is exactly unchanged, but each far-pair
# direction becomes a tangent-space generator of the expanded state, restoring
# the O(δt) seeding. The φ set is computed model-agnostically from the dense
# Hamiltonian: the support of Hχ outside the product-state tangent span
# (basis vectors at Hamming distance ≤ 1, or exactly 2 on ADJACENT sites).
#
# Metrics: baseline vs expanded sweep error Δgt/Δlt against the dense exact leg,
# plus a δt scan on config C to expose the convergence order.
#
#   julia --project=. scripts/diag_expansion.jl
# ==============================================================================

using Random, LinearAlgebra
using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum
using FermionicTCMPS: _sector, _add, _amplitude, _adag_up, _adag_down, _a_up, _a_down, Charge

include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))
include(joinpath(@__DIR__, "exact_validation.jl"))   # dense helpers, δt, TIMES, TRUNC

# ---------------- equivalent expansion construction ----------------

"per-site charge label of a particle-occupation tuple (hole sites relabeled)"
_site_label(spec, i, occ::Charge) =
    spec.kinds[i] === :hole ? (1 - occ[1], 1 - occ[2]) : (occ[1], occ[2])

"""
    expanded_state(spec, occs0, phis; base_vals, T) -> MPS

Represent the product basis state `occs0` exactly as the channel sum
|occs0⟩ ⊕ Σ_j (|φ_j⟩ ⊕ −|φ_j⟩): bond spaces are direct sums over channels (ψ
first, then per φ a ± pair), site tensors carry one 1×1 block entry per channel
(rows/cols ordered by channel construction order). Boundary coefficients fold
into the first/last site tensor. `base_vals[i]` is the ψ entry at site i (the
insertion amplitude on the impurity).
"""
function expanded_state(spec, occs0, phis; base_vals = nothing, T = ComplexF64)
    L = length(spec.kinds)
    base_vals = base_vals === nothing ? ones(L) : base_vals
    occ(ch) = ch[1]
    coeff(ch) = ch[2]
    channels = [(collect(occs0), 1.0)]
    for φ in phis
        push!(channels, (collect(φ), 1.0))
        push!(channels, (collect(φ), -1.0))
    end
    nch = length(channels)

    labels = [[_site_label(spec, i, occ(ch)[i]) for i in 1:L] for ch in channels]
    bondq = [Vector{Charge}(undef, L + 1) for _ in 1:nch]
    for c in 1:nch
        bondq[c][1] = (0, 0)
        for i in 1:L
            bondq[c][i + 1] = _add(bondq[c][i], labels[c][i])
        end
    end
    total = bondq[1][L + 1]
    all(==(total), (bondq[c][L + 1] for c in 1:nch)) ||
        throw(ArgumentError("φ directions must conserve the total charge"))

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
            # canonical codomain row: fusion-tree order = phys sector (outer, sorted) ×
            # bond sector (inner, sorted) × slot within the bond sector; at the left
            # boundary the bond is one-dimensional and the row is always 1
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
            val = coeff(channels[c]) * (c == 1 ? base_vals[i] : 1.0)
            block(W, key)[row, col] = T(val)
        end
        Ws[i] = W
    end
    return MPS(Vector{typeof(Ws[1])}(Ws))
end

"""
    residual_directions(spec, occs0, H; tol=1e-12) -> Vector{Vector{Charge}}

Support of H|occs0⟩ (flat) outside the product-state tangent span: basis vectors
differing from `occs0` at one site anywhere, or at exactly two ADJACENT sites,
are projected out (those amplitudes the two-site sweep can already seed). Each
surviving basis vector becomes one expansion direction.
"""
function residual_directions(spec, occs0, H::Matrix{<:Real}; tol = 1.0e-12)
    L = length(spec.kinds)
    occs_all = [_occ_of(L, idx) for idx in 1:size(H, 1)]   # site-major flat order
    v0 = zeros(ComplexF64, size(H, 1))
    v0[_flat(L, occs0)] = 1.0
    χ = H * v0
    out = Vector{Vector{Charge}}()
    for (idx, b) in enumerate(occs_all)
        abs(χ[idx]) <= tol && continue
        diff = [s for s in 1:L if b[s] != occs0[s]]
        seeded = isempty(diff) || length(diff) == 1 ||
                 (length(diff) == 2 && abs(diff[1] - diff[2]) == 1)
        seeded && continue
        push!(out, collect(b))
    end
    return out
end

"""
    structural_directions(spec, couples, imp) -> occs0 -> Vector{Vector{Charge}}

Model-level φ construction without a dense Hamiltonian: the pair directions the
hopping terms can populate in one application from `occs0` —
  v1-type: impurity spin σ emptied, a₁ site `s1` filled with σ;
  v2-type: impurity spin σ emptied, a₂ site `s2` particle σ removed —
kept only when the two sites are NON-adjacent (adjacent pairs are already
inside the two-site sweep window) and the occupations allow the move.
Must agree with `residual_directions` wherever the dense Hamiltonian is available.
"""
function structural_directions(spec, couples, imp)
    L = length(spec.kinds)
    return occs0 -> begin
        out = Vector{Vector{Charge}}()
        for (s1, s2, ω, v1, v2) in couples, σ in 1:2
            occs0[imp][σ] == 1 || continue
            if !iszero(v1) && abs(s1 - imp) != 1 && occs0[s1][σ] == 0
                φ = collect(occs0)
                φ[imp] = σ == 1 ? (0, occs0[imp][2]) : (occs0[imp][1], 0)
                φ[s1] = σ == 1 ? (1, occs0[s1][2]) : (occs0[s1][1], 1)
                φ in out || push!(out, φ)
            end
            if !iszero(v2) && abs(s2 - imp) != 1 && occs0[s2][σ] == 1
                φ = collect(occs0)
                φ[imp] = σ == 1 ? (0, occs0[imp][2]) : (occs0[imp][1], 0)
                φ[s2] = σ == 1 ? (0, occs0[s2][2]) : (occs0[s2][1], 0)
                φ in out || push!(out, φ)
            end
        end
        out
    end
end

"expanded insertion chain χ = d†_σ|ψ₀⟩ or η = d_σ|ψ₀⟩ (nothing if it kills the state)"
function expanded_shifted_state(spec, imp_up, imp_down, spin, delta, phis, T = ComplexF64)
    L = length(spec.kinds)
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

    base = Charge[charge for charge in
                  (i == spec.imp_site ? q_out :
                   (spec.kinds[i] === :hole ? (1, 1) : (0, 0)) for i in 1:L)]
    base_vals = ones(L)
    base_vals[spec.imp_site] = amp
    return expanded_state(spec, Tuple(base), phis; base_vals, T)
end

# ---------------- sweep leg with expanded chains ----------------

"Green's function via the production sweep, with optionally expanded chains.
Direction source: dense residual (H) or a structural callback phis_fn(occs0) -> Vector{Vector{Charge}}."
function gf_sweep_expanded(x, times; imp_up, imp_down, δt, trunc, expand, H = nothing, phis_fn = nothing)
    h = mpo(x)
    spec = FermionicTCMPS.chain_spec(x)
    steps = round.(Int, times ./ δt)

    chains = NTuple{3,Any}[]
    greater = NTuple{2,Any}[]
    lesser = NTuple{2,Any}[]
    for (spin, nσ) in ((:up, imp_up), (:down, imp_down))
        for (kind, delta) in ((:greater, 1), (:lesser, -1))
            q_out = spin === :up ? (imp_up + delta, imp_down) : (imp_up, imp_down + delta)
            ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || continue
            occs0 = ntuple(i -> i == spec.imp_site ? q_out :
                           (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), length(spec.kinds))
            phis = Vector{Vector{Charge}}()
            if expand
                phis = phis_fn === nothing ? residual_directions(spec, occs0, H) : phis_fn(occs0)
            end
            χ = expand ?
                expanded_shifted_state(spec, imp_up, imp_down, spin, delta, phis) :
                FermionicTCMPS._shifted_state(x, imp_up, imp_down, spin, delta)
            χ === nothing && continue
            cache = environments(h, χ)
            push!(chains, (χ, cache, kind))
            snapshot = expand ?
                expanded_shifted_state(spec, imp_up, imp_down, spin, delta, phis) :
                FermionicTCMPS._shifted_state(x, imp_up, imp_down, spin, delta)
            (kind === :greater ? push!(greater, (length(chains), snapshot)) :
             push!(lesser, (length(chains), snapshot)))
        end
    end

    ψ0 = separable_state(x; imp_up, imp_down)
    E = real(expectationvalue(h, ψ0))
    alg = FermionicTCMPS._tdvp2(δt, trunc)
    gt = zeros(ComplexF64, length(times))
    lt = zeros(ComplexF64, length(times))
    prev = 0
    for (i, n) in enumerate(steps)
        for _ in (prev + 1):n
            foreach(c -> sweep!(c[2], alg), chains)
        end
        prev = n
        gtphase = -im * exp(im * E * times[i])
        ltphase = im * exp(-im * E * times[i])
        for (cid, χ0) in greater
            gt[i] += gtphase * dot(χ0, chains[cid][1])
        end
        for (cid, η0) in lesser
            lt[i] += ltphase * conj(dot(η0, chains[cid][1]))
        end
    end
    return gt, lt
end

# dense Hamiltonian for the chain frame (kinds → anchor_util convention)
function _chain_dense_H(x)
    spec = FermionicTCMPS.chain_spec(x)
    ws, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    return dense_hamiltonian(spec.kinds, spec.imp_site,
                             [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)];
                             εd = x.εd, U = x.U)
end

"construction self-check: expanded state must equal the base state exactly"
function expansion_selfcheck(x, H; imp_up, imp_down)
    spec = FermionicTCMPS.chain_spec(x)
    worst = 0.0
    for spin in (:up, :down), delta in (1, -1)
        χ = FermionicTCMPS._shifted_state(x, imp_up, imp_down, spin, delta)
        χ === nothing && continue
        q_out = spin === :up ? (imp_up + delta, imp_down) : (imp_up, imp_down + delta)
        occs0 = ntuple(i -> i == spec.imp_site ? q_out :
                       (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), length(spec.kinds))
        phis = residual_directions(spec, occs0, H)
        isempty(phis) && continue
        println("      directions (spin=$spin delta=$delta): ",
                [φ[spec.imp_site] => "sites " * string([s for s in 1:length(spec.kinds) if φ[s] != occs0[s]]) for φ in phis])
        exp = expanded_shifted_state(spec, imp_up, imp_down, spin, delta, phis)
        d0 = abs(dot(χ, exp) - 1.0)
        dn = abs(dot(exp, exp) - 1.0)
        dφ = maximum(abs(dot(_basis_state(spec, Tuple(φ); T = ComplexF64), exp)) for φ in phis)
        println("   selfcheck spin=$spin delta=$delta: |⟨χ|exp⟩−1| = $d0, |⟨exp|exp⟩−1| = $dn, " *
                "max|⟨φ|exp⟩| = $dφ  ($(length(phis)) directions)")
        worst = max(worst, d0, dn, dφ)
    end
    return worst
end

# ---------------- config runner ----------------

"large-config runner: reference = the lossless exact_timeevolution leg (dense H infeasible)"
function run_big_config(name, x, couples, times; δt, imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    phis_fn = structural_directions(spec, couples, spec.imp_site)
    gt_ref = _et_gt(x, times; imp_up, imp_down)
    gt_base, lt_base = gf_sweep_expanded(x, times; imp_up, imp_down, δt,
                                         trunc = TRUNC, expand = false)
    gt_exp, lt_exp = gf_sweep_expanded(x, times; imp_up, imp_down, δt,
                                       trunc = TRUNC, expand = true, phis_fn)
    fmt(v) = round(v, digits = 6)
    println("== ", name, "  (imp = ", (imp_up, imp_down), ") ==")
    println("   |gt| ET at (0.1, 1.0): ", fmt(abs(gt_ref[2])), " / ", fmt(abs(gt_ref[end])))
    println("   baseline: Δgt(vs ET) = ", fmt(maximum(abs.(gt_base - gt_ref))),
            "   Δlt = ", fmt(maximum(abs.(lt_base))))
    println("   expanded: Δgt(vs ET) = ", fmt(maximum(abs.(gt_exp - gt_ref))),
            "   Δlt = ", fmt(maximum(abs.(lt_exp))))
    return nothing
end

function run_exp_config(name, x, couples, L; εd, U, imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd, U)
    size(H, 1) == 4^L || error("dim mismatch")

    worst = expansion_selfcheck(x, H; imp_up, imp_down)
    worst < 1e-8 || error("expansion equivalence FAILED for $name: $worst")

    occs = ntuple(i -> i == spec.imp_site ? (imp_up, imp_down) :
                  (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
    ψ0 = zeros(ComplexF64, 4^L)
    ψ0[_flat(L, occs)] = 1.0
    E0 = real(H[_flat(L, occs), _flat(L, occs)])

    χ_up = _insert_vector(ψ0, :up, 1, L, spec.imp_site)
    χ_dn = _insert_vector(ψ0, :down, 1, L, spec.imp_site)
    η_up = _insert_vector(ψ0, :up, -1, L, spec.imp_site)
    η_dn = _insert_vector(ψ0, :down, -1, L, spec.imp_site)
    overlaps, = _exact_overlaps(H, Any[χ_up, χ_dn, η_up, η_dn], TIMES)
    gt_exact = [-im * exp(im * E0 * t) * (overlaps[1][i] + overlaps[2][i])
                for (i, t) in enumerate(TIMES)]
    lt_exact = [im * exp(-im * E0 * t) * conj(overlaps[3][i] + overlaps[4][i])
                for (i, t) in enumerate(TIMES)]

    gt_base, lt_base = gf_sweep_expanded(x, TIMES; imp_up, imp_down, δt, trunc = TRUNC, expand = false, H)
    gt_exp, lt_exp = gf_sweep_expanded(x, TIMES; imp_up, imp_down, δt, trunc = TRUNC, expand = true, H)

    fmt(v) = round(v, digits = 6)
    Δgt_b = maximum(abs.(gt_base - gt_exact)); Δgt_e = maximum(abs.(gt_exp - gt_exact))
    Δlt_b = maximum(abs.(lt_base - lt_exact)); Δlt_e = maximum(abs.(lt_exp - lt_exact))
    println("== ", name, "  (imp = ", (imp_up, imp_down), ") ==")
    println("   baseline: Δgt = ", fmt(Δgt_b), "   Δlt = ", fmt(Δlt_b))
    println("   expanded: Δgt = ", fmt(Δgt_e), "   Δlt = ", fmt(Δlt_e))
    return (Δgt_base = Δgt_b, Δgt_exp = Δgt_e, Δlt_base = Δlt_b, Δlt_exp = Δlt_e)
end

function main()
    wsc = [0.4, 0.9]; fsc = [0.55, 0.35]

    println("# part 1: L=5 validation matrix (δt = $δt, dense-exact reference;" *
            " each config first runs the construction self-check)")
    # A. finite β, free
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 0.0, εd = -0.3, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_exp_config("A single 2mode finiteβ U=0", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 0.0)

    # B. finite β, interacting
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 1.7, εd = -0.3, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_exp_config("B single 2mode finiteβ U=1.7", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 1.7)

    # C. vacuum above μ, interacting — the stuck class
    x = SISBD(discretefermionicvacuum(wsc, fsc; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    couples_c = [(2k, 2k + 1, wsc[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(wsc)]
    run_exp_config("C single 2mode vacuum-above U=1.7", x, couples_c, 5;
        εd = 0.2, U = 1.7)
    # structural directions must coincide with the dense residual (C, χ_up)
    specc = FermionicTCMPS.chain_spec(x)
    occs_c = ntuple(i -> i == specc.imp_site ? (1, 0) :
                    (specc.kinds[i] === :hole ? (1, 1) : (0, 0)), length(specc.kinds))
    Hd_c = dense_hamiltonian(copy(specc.kinds), specc.imp_site, couples_c; εd = x.εd, U = x.U)
    sort(structural_directions(specc, couples_c, specc.imp_site)(occs_c)) ==
        sort(residual_directions(specc, occs_c, Hd_c)) ||
        error("structural directions disagree with the dense residual on C")
    println("   structural-direction check vs dense residual (C, χ_up): ✓")

    # D. straddling vacuum, interacting
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.7, εd = -0.3, J = 1.0)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    run_exp_config("D single 2mode straddling vacuum U=1.7", x,
        [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)], 5;
        εd = -0.3, U = 1.7)

    # E. double bath finite β, interacting
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
    imp = spec.imp_site
    couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in 1:length(wsl)],
                   [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in 1:length(wsr)])
    run_exp_config("E double 1+1mode finiteβ U=1.7", x, couples, 5; εd = -0.3, U = 1.7)

    # F. vacuum above μ, hot impurity
    x = SISBD(discretefermionicvacuum(wsc, fsc; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
    _, v1f, v2f = FermionicTCMPS.thermofield_couplings(x.bath)
    run_exp_config("F single 2mode vacuum-above U=1.7 hot imp", x,
        [(2k, 2k + 1, wsc[k], x.J * v1f[k], x.J * v2f[k]) for k in eachindex(wsc)], 5;
        εd = 0.2, U = 1.7, imp_up = 1)

    # G. double bath finite β, hot impurity
    x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    run_exp_config("G double 1+1mode finiteβ U=1.7 hot imp", x, couples, 5;
        εd = -0.3, U = 1.7, imp_up = 1)

    # δt scan on G: convergence order of the degraded class, baseline vs expanded
    println("\n# δt scan on G (double 1+1mode finiteβ U=1.7 hot imp)")
    xg = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(xg)
    Lc = length(spec.kinds)
    impg = spec.imp_site
    Hd = dense_hamiltonian(copy(spec.kinds), impg, couples; εd = xg.εd, U = xg.U)
    occs = ntuple(i -> i == impg ? (1, 0) :
                  (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), Lc)
    ψ0 = zeros(ComplexF64, 4^Lc); ψ0[_flat(Lc, occs)] = 1.0
    χ_up = _insert_vector(ψ0, :up, 1, Lc, impg)
    χ_dn = _insert_vector(ψ0, :down, 1, Lc, impg)
    η_up = _insert_vector(ψ0, :up, -1, Lc, impg)
    η_dn = _insert_vector(ψ0, :down, -1, Lc, impg)
    overlaps, = _exact_overlaps(Hd, Any[χ_up, χ_dn, η_up, η_dn], TIMES)
    E0 = real(Hd[_flat(Lc, occs), _flat(Lc, occs)])
    gt_exact = [-im * exp(im * E0 * t) * (overlaps[1][i] + overlaps[2][i])
                for (i, t) in enumerate(TIMES)]
    lt_exact = [im * exp(-im * E0 * t) * conj(overlaps[3][i] + overlaps[4][i])
                for (i, t) in enumerate(TIMES)]
    for dt in (0.04, 0.02, 0.01)
        gt_base, lt_base = gf_sweep_expanded(xg, TIMES; imp_up = 1, imp_down = 0,
                                             δt = dt, trunc = TRUNC, expand = false, H = Hd)
        gt_exp, lt_exp = gf_sweep_expanded(xg, TIMES; imp_up = 1, imp_down = 0,
                                           δt = dt, trunc = TRUNC, expand = true, H = Hd)
        println("   δt = $dt:  baseline Δgt = ", round(maximum(abs.(gt_base - gt_exact)); digits = 6),
                "  Δlt = ", round(maximum(abs.(lt_base - lt_exact)); digits = 6),
                "   |   expanded Δgt = ", round(maximum(abs.(gt_exp - gt_exact)); digits = 6),
                "  Δlt = ", round(maximum(abs.(lt_exp - lt_exact)); digits = 6))
    end

    # part 2: larger baths — the benchmark/GTEMPO-crosscheck regime. The dense leg
    # is infeasible here; the lossless exact_timeevolution leg is the reference.
    println("\n# part 2: larger baths (reference = lossless ET leg)")
    for K in (3, 4)
        wsb = [-0.9, -0.4, 0.5, 0.9][1:K]
        fsb = [0.55, 0.35, 0.45, 0.3][1:K]
        xb = SISBD(discretefermionicvacuum(wsb, fsb; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
        _, vb1, vb2 = FermionicTCMPS.thermofield_couplings(xb.bath)
        couples_b = [(2k, 2k + 1, wsb[k], xb.J * vb1[k], xb.J * vb2[k]) for k in eachindex(wsb)]
        run_big_config("vacuum $(K)mode U=1.7 (L=$(2K + 1))", xb, couples_b, TIMES;
                       δt = δt, imp_up = 0, imp_down = 0)
    end
    # the GTEMPO-crosscheck benchmark class: double bath, 2+1 modes, finite β
    wsb = [-0.9, -0.3, 0.5]
    fsb = [0.35, 0.45, 0.5]
    for β in (0.5, 1.0)
        xl = discretefermionicbath(wsb[1:2], fsb[1:2]; β = β, μ = 0.0)
        xr = discretefermionicbath(wsb[3:3], fsb[3:3]; β = β, μ = 0.0)
        x = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
        spec = FermionicTCMPS.chain_spec(x)
        wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
        wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
        imp = spec.imp_site
        couples_b = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in 1:length(wsl)],
                         [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in 1:length(wsr)])
        run_big_config("benchmark double bath K=3 β=$β U=1.7 (L=$(length(spec.kinds)))",
                       x, couples_b, TIMES; δt = δt, imp_up = 0, imp_down = 0)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
