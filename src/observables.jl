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
    return _widen_svectors(ψ)
end

"""
    gf_greater_lesser_t(x, times; imp_up=0, imp_down=0, δt, trunc, expand=false) -> (gt, lt)

Spin-summed real-time impurity Green's function of the thermo-field state,
G^>(t) = −i⟨d(t)d†(0)⟩ and G^<(t) = i⟨d†(0)d(t)⟩.

The separable TFD state is an eigenstate of the chain Hamiltonian (every hopping
term annihilates the fiducial bath), so the correlators reduce to same-charge
overlaps of freely evolved insertion chains,
G^>(t) = −i·e^{iEt}·⟨χ|χ(t)⟩ with χ = d†|ψ₀⟩ and
G^<(t) = +i·e^{−iEt}·conj⟨η|η(t)⟩ with η = d|ψ₀⟩,
where E = ⟨ψ₀|H|ψ₀⟩. Spin channels whose insertion annihilates the state contribute
zero. `times` must be sorted, non-negative and lie on the δt grid.

`expand = true` starts each insertion chain from its EQUIVALENT expansion
(|χ⟩ ⊕ |φ⟩ ⊕ (−|φ⟩), see `expanded_state`), which makes the far-pair amplitudes
the sharp state cannot seed (bath modes beyond the sweep window — the whole
double-bath class, and everything under a truncation without `add_back`)
first-order in the dynamics. With `expand = true` and no explicit `trunc`, the
verdict truncation `TRUNC_EXPANDED` (D = 10, add_back = 0) is used.
"""
function gf_greater_lesser_t(x::Union{SISBD,SIDBD}, times::AbstractVector{<:Real};
                             imp_up::Int = 0, imp_down::Int = 0, δt::Real,
                             expand::Bool = false,
                             trunc::Union{Nothing,DMRG.TruncationDimCutoff} = nothing)
    h = mpo(x)
    efftrunc = trunc === nothing ? (expand ? TRUNC_EXPANDED :
                                   truncdimcutoff(D = 256, ϵ = 1.0e-10, add_back = MIND)) : trunc

    steps = round.(Int, times ./ δt)
    all(abs(s * δt - t) ≤ 1e-8 * max(1.0, abs(t)) for (s, t) in zip(steps, times)) ||
        throw(ArgumentError("times must lie on the δt grid"))
    (issorted(times) && first(times) ≥ 0) || throw(ArgumentError("times must be sorted and non-negative"))

    ψ0 = separable_state(x; imp_up, imp_down)
    E = real(expectationvalue(h, ψ0))
    alg = _tdvp2(δt, efftrunc)
    phis_fn = expand ? far_pair_directions(chain_spec(x)) : nothing

    chains = NTuple{3,Any}[]            # (state, cache, kind::Symbol)
    greater = NTuple{2,Any}[]           # (chain index, initial state)
    lesser = NTuple{2,Any}[]
    for (spin, nσ) in ((:up, imp_up), (:down, imp_down))
        χ = _chain_state(x, imp_up, imp_down, spin, 1, expand, phis_fn)
        if χ !== nothing
            cache = environments(h, χ)
            push!(chains, (χ, cache, :greater))
            push!(greater, (length(chains), _chain_state(x, imp_up, imp_down, spin, 1, expand, phis_fn)))
        end
        η = _chain_state(x, imp_up, imp_down, spin, -1, expand, phis_fn)
        if η !== nothing
            cache = environments(h, η)
            push!(chains, (η, cache, :lesser))
            push!(lesser, (length(chains), _chain_state(x, imp_up, imp_down, spin, -1, expand, phis_fn)))
        end
    end

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
            state = chains[cid][1]
            gt[i] += gtphase * dot(χ0, state)
        end
        for (cid, η0) in lesser
            state = chains[cid][1]
            lt[i] += ltphase * conj(dot(η0, state))
        end
    end
    return gt, lt
end

"insertion chain, plain or equivalently expanded (φ directions need occs0)"
function _chain_state(x::Union{SISBD,SIDBD}, imp_up::Int, imp_down::Int,
                      spin::Symbol, delta::Int, expand::Bool, phis_fn)
    if !expand
        return _shifted_state(x, imp_up, imp_down, spin, delta)
    end
    spec = chain_spec(x)
    q_out = spin === :up ? (imp_up + delta, imp_down) : (imp_up, imp_down + delta)
    ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || return nothing
    occs0 = ntuple(i -> i == spec.imp_site ? q_out :
                   (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), length(spec.kinds))
    return expanded_shifted_state(x, imp_up, imp_down, spin, delta, phis_fn(occs0))
end
