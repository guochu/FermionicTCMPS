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
    return site_operator_mpo(spec.spaces, spec.imp_site, amps; kinds = spec.kinds)
end

"""
    current_op(x::Union{SISBD,SIDBD}; spin=:both)

MPO of the impurity→bath current operator J = Σ_j V'_j c_j† d, summed over the
thermo-field branches (the particle–hole representation turns the a₂ branch into
ordinary h.c.-type hoppings). `spin` selects one spin contribution.
The physical current is J(t) = 2·Im⟨J⟩(t).
"""
function current_op(x::Union{SISBD,SIDBD}; spin::Symbol = :both)
    spec = chain_spec(x)
    selected = filter(spec.channels) do c
        spinok = spin === :both || (spin === :up && c.name === :up_hc) || (spin === :down && c.name === :down_hc)
        spinok && endswith(String(c.name), "_hc")
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

Impurity insertion d^delta_σ |imp⟩ as a product MPS (bond dimensions unchanged,
total charge shifted), with the intra-site phase taken from the operator table.
Returns `nothing` when the insertion annihilates the state.
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
    return ψ
end

"""
    gf_greater_lesser_t(x, times; imp_up=0, imp_down=0, δt, trunc) -> (gt, lt)

Spin-summed real-time impurity Green's function in the thermofield state,
G^>(t) = −i⟨d(t)d†(0)⟩ and G^<(t) = i⟨d†(0)d(t)⟩.

Algorithm: forward evolutions in lockstep from t = 0 — the state itself, the d†-inserted
chains (greater) and the d-inserted chains (lesser) — sampled at every `times` entry via
a bond-dimension-2 insertion MPO expectation (one MPS×MPO×MPS contraction per sample).
`times` must be sorted, non-negative and lie on the δt grid. A spin channel whose
insertion annihilates the impurity state contributes zero and is skipped.
"""
function gf_greater_lesser_t(x::Union{SISBD,SIDBD}, times::AbstractVector{<:Real};
                             imp_up::Int = 0, imp_down::Int = 0, δt::Real,
                             trunc::DMRG.TruncationDimCutoff = truncdimcutoff(D = 256, ϵ = 1.0e-10))
    spec = chain_spec(x)
    h = mpo(x)
    pos = spec.imp_site

    steps = round.(Int, times ./ δt)
    all(abs(s * δt - t) ≤ 1e-8 * max(1.0, abs(t)) for (s, t) in zip(steps, times)) ||
        throw(ArgumentError("times must lie on the δt grid"))
    (issorted(times) && first(times) ≥ 0) || throw(ArgumentError("times must be sorted and non-negative"))

    d_ops = Dict(:up => site_operator_mpo(spec.spaces, pos, _a_up; kinds = spec.kinds),
                 :down => site_operator_mpo(spec.spaces, pos, _a_down; kinds = spec.kinds))
    ddag_ops = Dict(:up => site_operator_mpo(spec.spaces, pos, _adag_up; kinds = spec.kinds),
                    :down => site_operator_mpo(spec.spaces, pos, _adag_down; kinds = spec.kinds))

    ψ0 = separable_state(x; imp_up, imp_down)
    ψ = separable_state(x; imp_up, imp_down)
    ψcache = environments(h, ψ)
    alg = _tdvp2(δt, trunc)

    greater = NTuple{4,Any}[]   # (χ, cache, d-op MPO, spin)
    lesser = NTuple{4,Any}[]    # (η, cache, d†-op MPO, spin)
    for (spin, nσ) in ((:up, imp_up), (:down, imp_down))
        if nσ == 0
            χ = _shifted_state(x, imp_up, imp_down, spin, 1)
            χ === nothing || push!(greater, (χ, environments(h, χ), d_ops[spin], spin))
        else
            η = _shifted_state(x, imp_up, imp_down, spin, -1)
            η === nothing || push!(lesser, (η, environments(h, η), ddag_ops[spin], spin))
        end
    end

    gt = zeros(ComplexF64, length(times))
    lt = zeros(ComplexF64, length(times))
    prev = 0
    for (i, n) in enumerate(steps)
        for _ in (prev + 1):n
            sweep!(ψcache, alg)
            foreach(g -> sweep!(g[2], alg), greater)
            foreach(l -> sweep!(l[2], alg), lesser)
        end
        prev = n
        for (χ, _, op, _) in greater
            gt[i] += -im * expectation(ψ, op, χ)
        end
        for (η, _, op, _) in lesser
            lt[i] += im * expectation(ψ0, op, η)
        end
    end
    return gt, lt
end
