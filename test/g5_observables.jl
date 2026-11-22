# ==============================================================================
# G5 observable-level validation, two layers per user verdict:
#
#   L1  single-mode baths (one discrete mode per bath), interacting (U ≠ 0) —
#       reference = dense mechanical-JW ED of the chain (no strict solution
#       exists for the interacting model). Covers the greater/lesser GF, the
#       impurity occupation n_dσ and the current; the single-bath case anchors
#       the current at the full-matrix level (the dense image of current_op is
#       assembled FROM the MPO, whose elements are pinned by G2).
#   L2  non-interacting (U = 0) against the strict free-fermion solution (the
#       ImpurityModelBase star-frame family). The current is verified through
#       the exact conservation identity 2·Im⟨J⟩ = d⟨n_d↑ + n_d↓⟩/dt (with the
#       hot impurity, so the decay is nontrivial). IMB's NEQ greater/lesser are
#       (t, 0) two-point functions and do NOT reduce to the occupation for a
#       non-stationary ρ₀ — see the testset comment for the numerically pinned
#       convention.
#
# current_op is the ONE-DIRECTION operator Σ V'(c†d) (non-hermitian; the G2
# product-basis anchor is blind to a missing h.c. by construction), so the
# physical current is 2·Im⟨J⟩.
# ==============================================================================

using LinearAlgebra: eigen, Diagonal, Hermitian

"occupation tuple of a flat basis index (site 1 most significant, up-first)"
_g5_occ_of(L, idx) = ntuple(s -> begin
    d = (idx - 1) ÷ 4^(L - s) % 4
    (d == 1 || d == 3 ? 1 : 0, d >= 2 ? 1 : 0)
end, L)

"flat index of a particle-occupation tuple"
_g5_flat(L, occs) = 1 + sum((occs[s][1] + 2 * occs[s][2]) * 4^(L - s) for s in 1:L)

"dense exact greater AND lesser GF of the chain-frame fiducial state occs0"
function _g5_dense_gf(x, couples, times; imp_up = 0, imp_down = 0)
    spec = FermionicTCMPS.chain_spec(x)
    L = length(spec.kinds); imp = spec.imp_site
    H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
    occs0 = ntuple(i -> i == imp ? (imp_up, imp_down) :
                   (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
    ψ0 = zeros(ComplexF64, 4^L); ψ0[_g5_flat(L, occs0)] = 1.0
    E0 = real(H[_g5_flat(L, occs0), _g5_flat(L, occs0)])
    F = eigen(Hermitian(H))
    gt = zeros(ComplexF64, length(times))
    lt = zeros(ComplexF64, length(times))
    for spin in (:up, :down), (delta, acc) in ((1, gt), (-1, lt))
        q_out = spin === :up ? (imp_up + delta, imp_down) : (imp_up, imp_down + delta)
        ((0 <= q_out[1] <= 1) && (0 <= q_out[2] <= 1)) || continue
        tbl = if spin === :up
            delta == 1 ? FermionicTCMPS._adag_up : FermionicTCMPS._a_up
        else
            delta == 1 ? FermionicTCMPS._adag_down : FermionicTCMPS._a_down
        end
        χ = zeros(ComplexF64, 4^L); found = false
        for i in eachindex(ψ0)
            iszero(ψ0[i]) && continue
            amp = FermionicTCMPS._amplitude(tbl, q_out, _g5_occ_of(L, i)[imp])
            amp === nothing && continue
            χ[i + delta * (spin === :up ? 4^(L - imp) : 2 * 4^(L - imp))] += amp * ψ0[i]
            found = true
        end
        found || continue
        g = F.vectors' * χ
        for (j, t) in enumerate(times)
            ov = sum(abs2(g[k]) * exp(-im * F.values[k] * t) for k in eachindex(g))
            acc[j] += delta == 1 ? -im * exp(im * E0 * t) * ov :
                      im * exp(-im * E0 * t) * conj(ov)
        end
    end
    return gt, lt
end

"evolve the fiducial chain with sweeps, return observable values at times"
function _g5_sweep_obs(x, times; δt, trunc, imp_up = 0, imp_down = 0)
    h = mpo(x)
    occ = occupation_op(x; spin = :up)
    jop = current_op(x; spin = :both)
    ψ = separable_state(x; imp_up, imp_down)
    cache = environments(h, ψ)
    alg = FermionicTCMPS._tdvp2(δt, trunc)
    out = Dict{Float64,Tuple{Float64,ComplexF64}}()
    prev = 0
    for t in sort(times)
        n = round(Int, t / δt)
        for _ in (prev + 1):n; sweep!(cache, alg); end
        prev = n
        out[round(t; digits = 9)] = (real(expectationvalue(occ, ψ)),
                                     expectationvalue(jop, ψ))
    end
    return out
end

# ---------- L1: single-mode interacting (dense ED reference) ----------

@testset "G5 L1 single-mode Anderson: GF + occupation vs dense ED (SISBD)" begin
    Random.seed!(51)
    ws = [0.4]; fs = [0.55]
    δt = 0.01
    times = collect(0.05:0.05:0.3)
    trunc = truncdimcutoff(D = 64, ϵ = 1e-12)
    for (imp_up, imp_down) in ((0, 0), (1, 0))
        x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.7, εd = 0.2, J = 1.1)
        spec = FermionicTCMPS.chain_spec(x)
        L = length(spec.kinds); imp = spec.imp_site
        _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
        couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)]

        gt_d, lt_d = _g5_dense_gf(x, couples, times; imp_up, imp_down)
        gt, lt = gf_greater_lesser_t(x, times; δt, trunc, imp_up, imp_down)
        @test maximum(abs.(gt .- gt_d)) < 1e-4
        @test maximum(abs.(lt .- lt_d)) < 1e-4

        H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
        Fe = eigen(Hermitian(H))
        occs0 = ntuple(i -> i == imp ? (imp_up, imp_down) :
                       (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
        ψ0 = zeros(ComplexF64, 4^L); ψ0[_g5_flat(L, occs0)] = 1.0
        g0 = Fe.vectors' * ψ0
        nddiag = [(i - 1) ÷ 4^(L - imp) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]
        obs = _g5_sweep_obs(x, times; δt, trunc, imp_up, imp_down)
        for t in times
            ψt = Fe.vectors * (exp.(-im .* Fe.values .* t) .* g0)
            nd_dense = real(ψt' * Diagonal(nddiag) * ψt)
            @test abs(first(obs[round(t; digits = 9)]) - nd_dense) < 1e-4
        end
    end
end

@testset "G5 L1 single-mode Anderson: current vs dense ED (SISBD, full matrix)" begin
    # the dense image of current_op is assembled from the MPO itself (elements
    # pinned by G2), so this anchors sweep+expectation against exact evolution
    # of the same observable definition — complex value, not just the real part
    Random.seed!(52)
    ws = [0.4]; fs = [0.55]
    δt = 0.01
    times = collect(0.05:0.05:0.3)
    trunc = truncdimcutoff(D = 64, ϵ = 1e-12)
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.7, εd = 0.2, J = 1.1)
    spec = FermionicTCMPS.chain_spec(x)
    L = length(spec.kinds)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)]
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    F = eigen(Hermitian(H))
    jop = current_op(x; spin = :both)
    Jfull = zeros(4^L, 4^L)
    for w in 1:4^L
        ψw = _basis_state(spec, _g5_occ_of(L, w))
        for v in 1:4^L
            _qtotal(spec, _g5_occ_of(L, v)) == _qtotal(spec, _g5_occ_of(L, w)) || continue
            Jfull[v, w] = expectation(_basis_state(spec, _g5_occ_of(L, v)), jop, ψw)
        end
    end

    for (imp_up, imp_down) in ((0, 0), (1, 0))
        occs0 = ntuple(i -> i == spec.imp_site ? (imp_up, imp_down) :
                       (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
        ψ0 = zeros(ComplexF64, 4^L); ψ0[_g5_flat(L, occs0)] = 1.0
        g0 = F.vectors' * ψ0
        obs = _g5_sweep_obs(x, times; δt, trunc, imp_up, imp_down)
        for t in times
            ψt = F.vectors * (exp.(-im .* F.values .* t) .* g0)
            j_exact = ψt' * Jfull * ψt
            @test abs(last(obs[round(t; digits = 9)]) - j_exact) < 2e-4
        end
    end
end

@testset "G5 L1 single-mode-per-bath Anderson: double bath GF + occupation vs dense ED (SIDBD)" begin
    # the far pair (left bath, sites 1↔3) needs the equivalent initial-state
    # expansion for first-order seeding; with it the pipeline is exact to ~1e-6
    Random.seed!(53)
    δt = 0.01
    times = collect(0.05:0.05:0.3)
    xl = discretefermionicvacuum([-0.8], [0.6]; μ = 0.0)
    xr = discretefermionicvacuum([0.5], [0.45]; μ = 0.0)
    x = SIDBD(xl, xr; U = 1.0, εd = -0.3, Jl = 1.2, Jr = 0.9)
    spec = FermionicTCMPS.chain_spec(x)
    L = length(spec.kinds); imp = spec.imp_site
    wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
    couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in eachindex(wsl)],
                   [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in eachindex(wsr)])
    gt_d, lt_d = _g5_dense_gf(x, couples, times)
    gt, lt = gf_greater_lesser_t(x, times; δt, expand = true)
    @test maximum(abs.(gt .- gt_d)) < 1e-5
    @test maximum(abs.(lt .- lt_d)) < 1e-5

    H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
    Fe = eigen(Hermitian(H))
    occs0 = ntuple(i -> i == imp ? (0, 0) : (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
    ψ0 = zeros(ComplexF64, 4^L); ψ0[_g5_flat(L, occs0)] = 1.0
    g0 = Fe.vectors' * ψ0
    nddiag = [(i - 1) ÷ 4^(L - imp) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]
    obs = _g5_sweep_obs(x, times; δt, trunc = truncdimcutoff(D = 10, ϵ = 1e-10, add_back = 0))
    for t in times
        ψt = Fe.vectors * (exp.(-im .* Fe.values .* t) .* g0)
        nd_dense = real(ψt' * Diagonal(nddiag) * ψt)
        @test abs(first(obs[round(t; digits = 9)]) - nd_dense) < 1e-5
    end
end

# ---------- L2: non-interacting vs the ImpurityModelBase strict solution ----------

@testset "G5 L2 current + occupation vs the U=0 strict solution (vacuum, hot impurity)" begin
    # Exact conservation identity 2·Im⟨J⟩ = −d⟨n_d⟩/dt (i⟨[H,n]⟩ with
    # [H_hop, n_d] = V(c†d − d†c); positive J = impurity→bath outflow), checked
    # against the strict solution of the FREE chain: dense ED of the
    # mechanical-JW Hamiltonian (the U = 0 model is exactly solvable; the GF
    # layer against ImpurityModelBase.freefermions_greater_lesser is covered in
    # G3, whose (t, 0) two-point convention is what the GF pipeline targets —
    # NOTE IMB's NEQ greater/lesser are (t, 0) functions and do NOT reduce to
    # the equal-time occupation for a non-stationary ρ₀).
    #
    # The sweep MUST run with the equivalent initial-state expansion here: for
    # vacuum baths every impurity→bath decay channel of the hot fiducial chain
    # is a non-adjacent pair (the baseline sweep is completely stuck: n ≡ 1,
    # ⟨J⟩ ≡ 0 — that is the documented TDVP2 defect, not an observable error).
    # No spin factor: with imp = (1,0) only the up channels are populated.
    Random.seed!(54)
    δt = 0.01
    times = collect(0.1:0.1:0.4)

    configs = Tuple{String,Any,Vector{Tuple{Int,Int,Float64,Float64,Float64}}}[]
    let x = SISBD(discretefermionicvacuum([0.4], [0.55]; μ = 0.0);
                  U = 0.0, εd = -0.3, J = 1.1)
        _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
        push!(configs, ("vacuum K=1", x, [(2k, 2k + 1, [0.4][k], x.J * v1[k], x.J * v2[k])
                                         for k in eachindex(v1)]))
    end
    let x = SISBD(discretefermionicvacuum([-0.9, 0.5], [0.35, 0.55]; μ = 0.0);
                  U = 0.0, εd = -0.3, J = 1.1)
        _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
        push!(configs, ("vacuum K=2", x, [(2k, 2k + 1, [-0.9, 0.5][k], x.J * v1[k], x.J * v2[k])
                                         for k in eachindex(v1)]))
    end
    let x = SIDBD(discretefermionicvacuum([-0.8], [0.6]; μ = 0.0),
                  discretefermionicvacuum([0.5], [0.45]; μ = 0.0);
                  U = 0.0, εd = -0.3, Jl = 1.0, Jr = 1.0)
        spec = FermionicTCMPS.chain_spec(x); imp = spec.imp_site
        wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
        wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
        couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in eachindex(wsl)],
                       [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in eachindex(wsr)])
        push!(configs, ("double vacuum 1+1", x, couples))
    end

    for (name, x, couples) in configs
        spec = FermionicTCMPS.chain_spec(x)
        L = length(spec.kinds); imp = spec.imp_site
        H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
        F = eigen(Hermitian(H))
        occs0 = ntuple(i -> i == imp ? (1, 0) :
                       (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
        ψ0 = zeros(ComplexF64, 4^L); ψ0[_g5_flat(L, occs0)] = 1.0
        g0 = F.vectors' * ψ0
        nddiag = [(i - 1) ÷ 4^(L - imp) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]
        n_strict(t) = begin
            ψt = F.vectors * (exp.(-im .* F.values .* t) .* g0)
            real(ψt' * Diagonal(nddiag) * ψt)
        end
        # strict current by central differences of the strict occupation
        hfd = 0.005
        jref = [(n_strict(t - hfd) - n_strict(t + hfd)) / (2 * hfd) for t in times]

        # sweep side with the equivalent expansion
        h = mpo(x)
        jop = current_op(x; spin = :both)
        occ = occupation_op(x; spin = :up)
        ψ = FermionicTCMPS.expanded_state(spec, occs0,
                                          FermionicTCMPS.far_pair_directions(spec)(occs0))
        cache = environments(h, ψ)
        alg = FermionicTCMPS._tdvp2(δt, FermionicTCMPS.TRUNC_EXPANDED)
        steps = round.(Int, times ./ δt)
        prev = 0
        for (i, t) in enumerate(times)
            for _ in (prev + 1):steps[i]; sweep!(cache, alg); end
            prev = steps[i]
            nv = real(expectationvalue(occ, ψ))
            jd = 2 * imag(expectationvalue(jop, ψ))
            @test abs(nv - n_strict(t)) < 1e-5
            @test abs(jd - jref[i]) < 5e-4   # FD truncation of the reference
        end
    end
end
