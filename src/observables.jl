# ==============================================================================
# Observables: impurity occupation, current, real-time Green's function.
# ==============================================================================

"""
    occupation_op(x::Union{SISBD,SIDBD}; spin=:up)

MPO of the impurity number operator n_dσ (`spin = :up` or `:down`).
Measure with `DMRG.expectationvalue(op, ψ)`.
"""
function occupation_op(x::Union{SISBD,SIDBD}; spin::Symbol = :up)
    (spin === :up || spin === :down) || throw(ArgumentError("spin must be :up or :down"))
    spec = chain_spec(x)
    amps = spin === :up ? _n_up : _n_down
    return site_operator_mpo(spec.spaces, spec.imp_site, amps)
end

"""
    current_op(x::Union{SISBD,SIDBD}; spin=:both)

MPO of the impurity→bath current operator J = Σ_j V'_j c_j† d, summed over the
thermo-field branches (both branches carry number-conserving hoppings with the
thermo-field weights). `spin` selects one spin contribution.
The physical current is J(t) = 2·Im⟨J⟩(t).
"""
function current_op(x::Union{SISBD,SIDBD}; spin::Symbol = :both)
    spec = chain_spec(x)
    selected = filter(spec.channels) do c
        spinok = spin === :both || (spin === :up && c.name === :up_hc) || (spin === :down && c.name === :down_hc)
        spinok && endswith(String(c.name), "_hc") &&
            any(e.site == spec.imp_site for e in (c.opens..., c.closes...))
    end
    isempty(selected) && throw(ArgumentError("no current channels selected for spin = $spin"))
    return _fsm_mpo(spec, selected; onsite = false)
end

# amplitude lookup in an operator table; `nothing` marks a zero matrix element
function _amplitude(t::OpTable, q_out::Charge, q_in::Charge)
    idx = findfirst(p -> p.first == (q_out, q_in), t)
    return idx === nothing ? nothing : t[idx].second
end

"""
    _shifted_state(x, imp_up, imp_down, spin, delta)

Impurity insertion d^delta_σ |imp⟩ as a product MPS of the vacuum-fiducial chain
(bond dimensions unchanged, total charge shifted), with the intra-site phase
taken from the operator table. Returns `nothing` when the insertion annihilates
the state.
"""
function _shifted_state(x::Union{SISBD,SIDBD}, imp_up::Int, imp_down::Int, spin::Symbol, delta::Int)
    q_in = (imp_up, imp_down)
    q_out = spin === :up ? (q_in[1] + delta, q_in[2]) : (q_in[1], q_in[2] + delta)
    tbl = if spin === :up
        delta == 1 ? _adag_up : _a_up
    else
        delta == 1 ? _adag_down : _a_down
    end
    amp = _amplitude(tbl, q_out, q_in)
    amp === nothing && return nothing
    spec = chain_spec(x)
    sectors = [_sector(k == spec.imp_site ? q_out : (0, 0)) for k in 1:length(spec.kinds)]
    ψ = prodmps(ComplexF64, spec.spaces, sectors)
    ψ[spec.imp_site] = amp * ψ[spec.imp_site]
    return _widen_svectors(ψ)
end

# ⟨ψ|d_σ|χ⟩ evaluated as ⟨d†ψ|χ⟩: apply the impurity creation mechanically to the
# ket (impurity tensor rebuilt, every downstream bond sector shifted by the operator
# charge), then a same-charge dot. This sidesteps the bridging-charge three-way
# contraction, which neither the generic DMRG expectation boundary nor a hand-rolled
# one represents as a valid graded tensor.
#
# The data move runs at the FUSION-TREE level (not per block): on an evolved chain
# the bond sectors are multi-dimensional, the fused-charge blocks of the shifted
# tensor group a DIFFERENT set of (bond, phys) fusion trees than the original
# block, so whole-block copies are structurally wrong (they silently skip or
# misplace rows). Each tree entry moves individually: (ql, p, cr) →
# (ql ⊕ q_d, p, cr ⊕ q_d) past the impurity, and at the impurity
# (ql, p, cr) → (ql, p ⊕ q_d, cr ⊕ q_d) with the operator-table amplitude times
# the mechanical Jordan–Wigner string of c̃†_imp across the left part,
# (−1)^{n_left} = (−1)^{ql₁+ql₂} (ql = accumulated left charge of the tree).
# The tree-pair subblocks are strided views, and sector dims are preserved by the
# label shift, so within-tree slot order carries over.
function _apply_dag(x::Union{SISBD,SIDBD}, ψ::DMRG.MPS, spin::Symbol, delta::Int;
                    imp_up::Int = 0, imp_down::Int = 0)
    spec = chain_spec(x)
    imp = spec.imp_site
    L = length(spec.kinds)
    lookup = Dict(_sector((q1, q2)) => (q1, q2) for q1 in 0:(2L) for q2 in 0:(2L))
    q_in = (imp_up, imp_down)
    q_out = spin === :up ? (q_in[1] + delta, q_in[2]) : (q_in[1], q_in[2] + delta)
    tbl = spin === :up ? (delta == 1 ? _adag_up : _a_up) : (delta == 1 ? _adag_down : _a_down)
    _amplitude(tbl, q_out, q_in) === nothing && error("insertion annihilates the state")
    q_d = q_out .- q_in
    # strip the canonical-form bond tensors: the shifted data must not inherit
    # the unshifted .s chain
    B = DMRG.MPSBondTensor{spacetype(ψ[1])}
    sv = Vector{Union{Missing,B}}(missing, length(ψ))
    φ = DMRG.MPS(deepcopy(ψ.data), sv)
    for k in imp:length(ψ)
        T = φ[k]
        lft = k > imp ? _shifted_space(space(T, 1), q_d, lookup) : space(T, 1)
        rgt = _shifted_space(space(T, 3), q_d, lookup)
        T2 = zeros(ComplexF64, lft ⊗ space(T, 2) ← rgt)
        dst = Dict{Tuple{Charge,Charge,Charge},Any}()
        for (f1, f2) in fusiontrees(T2)
            dst[(_charge_of(f1.uncoupled[1], lookup), _charge_of(f1.uncoupled[2], lookup),
                 _charge_of(f2.uncoupled[1], lookup))] = (f1, f2)
        end
        for (f1, f2) in fusiontrees(T)
            ql = _charge_of(f1.uncoupled[1], lookup)
            p = _charge_of(f1.uncoupled[2], lookup)
            cr = _charge_of(f2.uncoupled[1], lookup)
            ql2 = k > imp ? _add(ql, q_d) : ql
            p2 = k == imp ? _add(p, q_d) : p
            val = T[f1, f2]
            if k == imp
                a = _amplitude(tbl, p2, p)
                (a === nothing || iszero(a)) && continue
                a *= isodd(ql[1] + ql[2]) ? -1.0 : 1.0
                val = a .* val
            end
            d = get(dst, (ql2, p2, _add(cr, q_d)), nothing)
            d === nothing && continue
            T2[d...] .= val
        end
        φ[k] = T2
    end
    return φ
end

"sector label → charge: bond charges are non-negative particle counts, and the
stored labels may dualize each U₁ component independently (the exact-timeevolution
surgery does), so decode by absolute irrep charge"
_charge_of(s, lookup) = (abs(Int(s[1].charge)), abs(Int(s[2].charge)))

"bond space with every sector shifted by q"
function _shifted_space(sp, q::Charge, lookup)
    elems = Pair{eltype(sp),Int}[]
    for s in sectors(sp)
        c = _charge_of(s, lookup) .+ q
        push!(elems, _sector(c) => dim(sp, s))
    end
    isempty(elems) && return oneunit(Rep[U₁×U₁])
    return foldl(⊕, (Rep[U₁×U₁](e) for e in elems))
end

"""
    gf_greater_lesser_t(x, times; imp_up=0, imp_down=0, δt, trunc) -> (gt, lt)

Spin-summed real-time impurity Green's function of the thermo-field state,
G^>(t) = −i⟨d(t)d†(0)⟩ and G^<(t) = i⟨d†(0)d(t)⟩.

With the thermofield-double fiducial the chain state is NOT a Hamiltonian
eigenstate (the branch hoppings act on the doubles), so the correlators are
sampled as sandwiches of forward-evolved states:
G^>(t) = −i·⟨ψ(t)| d_σ |χ(t)⟩ with χ = d†|ψ₀⟩ and
G^<(t) = +i·⟨η(t)| d_σ |ψ(t)⟩ with η = d|ψ₀⟩,
evolved in lockstep. Spin channels whose insertion annihilates the state
contribute zero. `times` must be sorted, non-negative and lie on the δt grid.

The sweep is the plain route: sharp thermofield-double fiducial, mechanical
insertions, two-site TDVP with the KS minimum-bond crutch (`add_back = MIND`).
The two-site projection seeds the far-pair amplitudes only at O(δt·v), so the
scheme is effectively FIRST order in δt (documented TDVP2 projection defect;
class-wise error budgets: vacuum ≈ exact, finite-β ~1e-3–5e-3, sharp sea ~1e-2).
"""
function gf_greater_lesser_t(x::Union{SISBD,SIDBD}, times::AbstractVector{<:Real};
                             imp_up::Int = 0, imp_down::Int = 0, δt::Real,
                             trunc::Union{Nothing,DMRG.TruncationDimCutoff} = nothing)
    h = mpo(x)

    steps = round.(Int, times ./ δt)
    all(abs(s * δt - t) ≤ 1e-8 * max(1.0, abs(t)) for (s, t) in zip(steps, times)) ||
        throw(ArgumentError("times must lie on the δt grid"))
    (issorted(times) && first(times) ≥ 0) || throw(ArgumentError("times must be sorted and non-negative"))

    spec = chain_spec(x)
    ψ = fiducial_state(x; imp_up, imp_down)
    dmpo = Dict(:up => site_operator_mpo(spec.spaces, spec.imp_site, _a_up),
                :down => site_operator_mpo(spec.spaces, spec.imp_site, _a_down))

    ψcaches = Any[environments(h, c) for c in ψ.comps]   # sweep! mutates ψ.comps in place
    inserts = NamedTuple{(:kind, :weights, :comps, :caches, :sign, :spin, :qeta),
                         Tuple{Symbol,Any,Any,Any,Float64,Symbol,Any}}[]
    for spin in (:up, :down)
        χ = _inserted_state(x, spin, 1; ψ, imp_up, imp_down)
        if χ !== nothing
            push!(inserts, (kind = :gt, weights = χ.weights, comps = χ.comps,
                            caches = Any[environments(h, c) for c in χ.comps],
                            sign = χ.phase * ψ.phase, spin = spin, qeta = (imp_up, imp_down)))
        end
        η = _inserted_state(x, spin, -1; ψ, imp_up, imp_down)
        if η !== nothing
            # η's impurity occupations are the annihilated ones
            qeta = spin === :up ? (imp_up - 1, imp_down) : (imp_up, imp_down - 1)
            push!(inserts, (kind = :lt, weights = η.weights, comps = η.comps,
                            caches = Any[environments(h, c) for c in η.comps],
                            sign = η.phase * ψ.phase, spin = spin, qeta = qeta))
        end
    end

    efftrunc = trunc === nothing ? truncdimcutoff(D = 256, ϵ = 1.0e-10, add_back = MIND) : trunc
    alg = _tdvp2(δt, efftrunc)

    gt = zeros(ComplexF64, length(times))
    lt = zeros(ComplexF64, length(times))
    prev = 0
    for (i, n) in enumerate(steps)
        for _ in (prev + 1):n
            foreach(c -> sweep!(c, alg), ψcaches)
            foreach(ins -> foreach(c -> sweep!(c, alg), ins.caches), inserts)
        end
        prev = n
        # the evolved chains' final bonds may carry per-sector multiplicity slots
        # (sweep splitting on top of the expansion channels): collapse to the
        # canonical charge labels, and restrict both sandwich sides to the common
        # charge set (marginal sectors truncated on one side only contribute noise)
        ψt = [_collapse_final(c.state) for c in ψcaches]
        for ins in inserts
            χt = [_collapse_final(c.state) for c in ins.caches]
            # ⟨ψ|d|χ⟩ = ⟨d†ψ|χ⟩ — same-charge overlaps of mechanically shifted states
            if ins.kind === :gt
                acc = sum(w * begin
                    bra = _apply_dag(x, ψt[j], ins.spin, 1; imp_up, imp_down)
                    common = intersect(_final_charges(bra), _final_charges(χt[j]))
                    dot(_restrict_final(bra, common), _restrict_final(χt[j], common))
                end for (j, w) in enumerate(ins.weights))
                gt[i] += -im * ins.sign * acc
            else
                acc = sum(w * begin
                    bra = _apply_dag(x, χt[j], ins.spin, 1;
                                     imp_up = ins.qeta[1], imp_down = ins.qeta[2])
                    common = intersect(_final_charges(bra), _final_charges(ψt[j]))
                    dot(_restrict_final(bra, common), _restrict_final(ψt[j], common))
                end for (j, w) in enumerate(ins.weights))
                lt[i] += im * ins.sign * acc
            end
        end
    end
    return gt, lt
end

"insertion d^delta_σ|ψ₀⟩ of the fiducial state; nothing when the insertion annihilates it.
Each fiducial component is shifted mechanically (exact per-configuration
Jordan–Wigner strings)."
function _inserted_state(x::Union{SISBD,SIDBD}, spin::Symbol, delta::Int;
                         ψ::SectorSumState, imp_up::Int = 0, imp_down::Int = 0)
    q_in = (imp_up, imp_down)
    q_out = spin === :up ? (q_in[1] + delta, q_in[2]) : (q_in[1], q_in[2] + delta)
    ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || return nothing
    tbl = spin === :up ? (delta == 1 ? _adag_up : _a_up) : (delta == 1 ? _adag_down : _a_down)
    amp = _amplitude(tbl, q_out, q_in)
    amp === nothing && return nothing
    q_d = q_out .- q_in
    comps = [_apply_dag(x, c, spin, delta; imp_up, imp_down) for c in ψ.comps]
    return SectorSumState(comps, copy(ψ.weights), [q .+ q_d for q in ψ.charges], ψ.phase)
end
