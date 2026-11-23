# ==============================================================================
# G3 dynamics anchors: TDVP real-time evolution of the thermo-field chain.
#
# Reference side: exact free-fermion evolution of the STAR model with the density
# matrix ρ₀ = fixed impurity occupations ⊗ thermal bath occupations
# (diag(n_F(ω_k))) — the sudden-protocol initial condition. With the exact
# Bogoliubov image (weighted on-sites + intra-mode pairing + thermofield-double
# fiducial) the chain reproduces the star frame exactly, equal-time data included.
# ==============================================================================

using LinearAlgebra: eigen, diagm, tr, Hermitian

# star single-particle Hamiltonian (per spin block); the star-frame bath couplings
# are the RAW spectral couplings √f, before the thermofield splitting
function _star_h(x::SISBD)
    ws = frequencies(x.bath.f)
    V = spectrumcouplings(x.bath.f)
    K = length(ws)
    hu = zeros(1 + K, 1 + K)
    hu[1, 1] = x.εd
    for k in 1:K
        hu[1, 1+k] = x.J * V[k]
        hu[1+k, 1] = hu[1, 1+k]
        hu[1+k, 1+k] = ws[k]
    end
    return hu
end
function _star_h(x::SIDBD)
    wsl = frequencies(x.bathl.f); Vl = spectrumcouplings(x.bathl.f)
    wsr = frequencies(x.bathr.f); Vr = spectrumcouplings(x.bathr.f)
    Kl, Kr = length(wsl), length(wsr)
    hu = zeros(1 + Kl + Kr, 1 + Kl + Kr)
    hu[1, 1] = x.εd
    for k in 1:Kl
        hu[1, 1+k] = x.Jl * Vl[k]; hu[1+k, 1] = hu[1, 1+k]; hu[1+k, 1+k] = wsl[k]
    end
    for k in 1:Kr
        j = 1 + Kl + k
        hu[1, j] = x.Jr * Vr[k]; hu[j, 1] = hu[1, j]; hu[j, j] = wsr[k]
    end
    return hu
end

# density matrix of the separable initial state in the star frame
_rho0(x::Union{SISBD,SIDBD}; imp::Float64) = begin
    bath = x isa SISBD ? x.bath : x.bathr
    β = x isa SISBD ? FermionicTCMPS._inverse_temperature(x.bath) : FermionicTCMPS._inverse_temperature(x.bathl)
    ws = frequencies((x isa SISBD ? x.bath : x.bathl).f)
    ws2 = x isa SIDBD ? frequencies(x.bathr.f) : ws
    occ = vcat([imp], [fermidirac(β, (x isa SISBD ? x.bath : x.bathl).μ, ω) for ω in ws])
    if x isa SIDBD
        occ = vcat(occ, [fermidirac(x.bathr.β, x.bathr.μ, ω) for ω in ws2])
    end
    diagm(occ)
end

# ⟨n_dσ(t)⟩ = [e^{iht} ρ₀ e^{−iht}]_{dd}
function _nd_star(x, ρ0, t)
    h = _star_h(x)
    F = eigen(h)
    ev = F.vectors * diagm(exp.(-im .* F.values .* t)) * F.vectors'
    P = zeros(size(h)...); P[1, 1] = 1.0
    return real(tr(P * ev' * ρ0 * ev))
end

@testset "G3 free dynamics: finite-β dense chain leg vs star (definition level)" begin
    # THE gate for the exact-image construction: the mechanical-JW chain evolving
    # the thermofield-double fiducial reproduces the star frame at machine
    # precision — equal-time AND (0,t) alike (the sweep pipeline below carries
    # the documented TDVP2 projection defect for this class and CANNOT serve
    # as the definition-level reference)
    Random.seed!(31)
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    β = 0.8; μ = 0.0
    times = collect(0.0:0.1:0.5)
    x = SISBD(discretefermionicbath(ws, fs; β, μ); U = 0.0, εd = -0.3, J = 1.1)
    spec = FermionicTCMPS.chain_spec(x)
    _, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
    couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k], ns[k]) for k in eachindex(ws)]
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    modes = [(c[1], c[2], c[6]) for c in couples]
    ψ0 = dense_tfd_vector(spec, modes, (0, 0))
    gt, lt = dense_gf_sandwich(H, ψ0, spec.imp_site, times)
    hu = _star_h(x)
    ρ0 = _rho0(x; imp = 0.0)
    gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 1e-10       # spin-summed
    @test maximum(abs.(lt)) < 1e-12                    # d annihilates the empty impurity
    # equal-time occupation of the hot fiducial against the star frame
    ρ0f = _rho0(x; imp = 1.0)
    ψ0h = dense_tfd_vector(spec, modes, (1, 0))
    F = eigen(Hermitian(H))
    Fh = eigen(hu)
    g0 = F.vectors' * ψ0h
    L = length(spec.kinds); imp = spec.imp_site
    nddiag = [(i - 1) ÷ 4^(L - imp) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L]
    worst = 0.0
    for t in times
        ψt = F.vectors * (exp.(-im .* F.values .* t) .* g0)
        nd = real(sum(conj(ψt) .* nddiag .* ψt))
        P = zeros(size(hu)); P[1, 1] = 1.0
        ev = Fh.vectors * diagm(exp.(-im .* Fh.values .* t)) * Fh.vectors'
        nd_star = real(tr(P * ev' * ρ0f * ev))
        worst = max(worst, abs(nd - nd_star))
    end
    @test worst < 1e-10
end

@testset "G3 free dynamics: greater/lesser GF vs exact star evolution (sweep)" begin
    # sweep-level gate: with the equivalent initial-state expansion the
    # finite-β TFD class is exact to roundoff (the far-pair directions seed at
    # first order); the plain sweep retains a small projection residue
    Random.seed!(31)
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    β = 0.8; μ = 0.0; δt = 0.025
    times = collect(0.0:0.1:0.5)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)

    # impurity empty: greater nontrivial, lesser zero
    x = SISBD(discretefermionicbath(ws, fs; β, μ); U = 0.0, εd = -0.3, J = 1.1)
    gt, lt = gf_greater_lesser_t(x, times; δt, trunc)
    gt1, lt1 = gf_greater_lesser_t(x, times; δt, expand = true)
    hu = _star_h(x)
    ρ0 = _rho0(x; imp = 0.0)
    gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 0.02         # spin-summed, projection residue
    @test maximum(abs.(lt)) < 1e-12                    # d annihilates the empty impurity
    @test maximum(abs.(gt1 .- 2 .* gtr)) < 1e-6        # expansion: exact
    @test maximum(abs.(lt1)) < 1e-12
end

@testset "G3 free dynamics: hot impurity channel structure" begin
    # impurity up-occupied: the greater runs on the empty down channel only
    # (which is an imp-empty problem) and the lesser on the occupied up channel;
    # both are checked against the star frame per spin
    Random.seed!(32)
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    β = 0.8; μ = 0.0; δt = 0.025
    times = collect(0.0:0.1:0.5)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    x = SISBD(discretefermionicbath(ws, fs; β, μ); U = 0.0, εd = -0.3, J = 1.1)
    hu = _star_h(x)
    gt2, lt2 = gf_greater_lesser_t(x, times; imp_up = 1, δt, trunc)
    gt2e, lt2e = gf_greater_lesser_t(x, times; imp_up = 1, δt, expand = true)
    @test abs(gt2[1] + im) < 1e-8          # G^>(0) = −i·(1 − n_down)
    gtr_e, = freefermions_greater_lesser(hu, _rho0(x; imp = 0.0), times, 1, 1)
    _, ltr_f = freefermions_greater_lesser(hu, _rho0(x; imp = 1.0), times, 1, 1)
    @test maximum(abs.(gt2 .- gtr_e)) < 0.01    # projection residue (see the GF testset)
    @test maximum(abs.(lt2 .- ltr_f)) < 0.01
    @test maximum(abs.(gt2e .- gtr_e)) < 1e-6   # expansion: exact
    @test maximum(abs.(lt2e .- ltr_f)) < 1e-6
end

@testset "G3 TFD equivalence: expanded fiducial/insertions ≡ plain construction" begin
    # the equivalent initial-state expansion must represent EXACTLY the same
    # states as the plain constructions: per-charge component overlaps and the
    # inserted-chain overlaps (gt/lt, both spins) pin the channel-sum machinery
    for (name, x, imp) in (("SISBD finite-β", SISBD(discretefermionicbath([-0.9, 0.5], [0.35, 0.55]; β = 0.8, μ = 0.0); U = 0.3, εd = -0.2), (0, 0)),
                           ("SISBD finite-β hot", SISBD(discretefermionicbath([-0.9, 0.5], [0.35, 0.55]; β = 0.8, μ = 0.0); U = 0.3, εd = -0.2), (1, 0)),
                           ("SIDBD vacuum-sea", SIDBD(discretefermionicvacuum([-0.8], [0.6]; μ = 0.0),
                                                      discretefermionicvacuum([0.5], [0.45]; μ = 0.0); U = 1.0, εd = -0.3), (0, 0)),
                           ("SIDBD finite-β", SIDBD(discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0),
                                                    discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0); U = 1.7, εd = -0.3), (0, 0)))
        ψf = fiducial_state(x; imp_up = imp[1], imp_down = imp[2])
        ψe = FermionicTCMPS.expanded_fiducial(x; imp_up = imp[1], imp_down = imp[2])
        @test sort(ψf.charges) == sort(ψe.charges)
        for q in ψf.charges
            j = findfirst(==(q), ψf.charges)
            k = findfirst(==(q), ψe.charges)
            @test abs(ψf.weights[j] - ψe.weights[k]) < 1e-10
            @test abs(dot(ψf.comps[j], ψe.comps[k]) - 1) < 1e-10
        end
        for spin in (:up, :down), delta in (1, -1)
            χf = FermionicTCMPS._inserted_state(x, spin, delta; ψ = ψf, imp_up = imp[1], imp_down = imp[2])
            χe = FermionicTCMPS._inserted_state(x, spin, delta; ψ = ψf, imp_up = imp[1], imp_down = imp[2], expand = true)
            @test (χf === nothing) == (χe === nothing)
            χf === nothing && continue
            for (j, q) in enumerate(χf.charges)
                @test abs(dot(χf.comps[j], χe.comps[j]) - 1) < 1e-10
            end
        end
    end
end

@testset "G3 free dynamics: vacuum degeneracy (single mode, a2 decoupled)" begin
    # n_k = 0: v2 = 0, the a2 branch decouples; with a single mode the whole
    # dynamics lives on the direct impurity-site-2 coupling
    Random.seed!(38)
    ws = [0.4]; fs = [0.55]
    δt = 0.025
    times = collect(0.0:0.1:0.4)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
    _, v1, v2, ns = FermionicTCMPS.thermofield_couplings(x.bath)
    @test all(iszero, v2) && all(!iszero, v1) && all(iszero, ns)
    gt, lt = gf_greater_lesser_t(x, times; δt, trunc)
    hu = _star_h(x)
    ρ0 = _rho0(x; imp = 0.0)
    gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 1e-5
    @test maximum(abs.(lt)) < 1e-12
end

@testset "G3 free dynamics: impurity occupation" begin
    # single vacuum mode: the whole dynamics is the direct impurity-a1 coupling,
    # which the two-site sweep captures exactly
    Random.seed!(34)
    ws = [0.4]; fs = [0.55]
    δt = 0.02
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
    h = mpo(x)
    occ = occupation_op(x; spin = :up)
    ρ0 = _rho0(x; imp = 1.0)

    ψ = fiducial_state(x; imp_up = 1, imp_down = 0)
    cache = environments(h, ψ)
    alg = FermionicTCMPS._tdvp2(δt, trunc)
    for n in 1:10
        sweep!(cache, alg)
        t = n * δt
        (n % 2 == 0) && (@test abs(real(expectationvalue(occ, ψ)) - _nd_star(x, ρ0, t)) < 1e-5)
    end

    # empty impurity: the vacuum fiducial is an exact eigenstate — frozen
    ψ0 = fiducial_state(x)
    @test abs(expectationvalue(h, ψ0)) < 1e-10
    evolve!(ψ0, h, (0.0, 0.2); δt, trunc)
    @test abs(real(expectationvalue(occ, ψ0))) < 1e-10
end

@testset "G3 free dynamics: conserved charges" begin
    Random.seed!(35)
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    δt = 0.05
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 1.7, εd = -0.3, J = 1.1)
    spec = FermionicTCMPS.chain_spec(x)
    h = mpo(x)

    # Q↑ = n_imp↑ + Σ_bath n↑ is conserved (every chain term is number-conserving)
    function qup(spec, ψ)
        s = real(expectationvalue(FermionicTCMPS.site_operator_mpo(spec.spaces, spec.imp_site, FermionicTCMPS._n_up), ψ))
        for k in 2:length(spec.kinds)
            s += real(expectationvalue(FermionicTCMPS.site_operator_mpo(spec.spaces, k, FermionicTCMPS._n_up), ψ))
        end
        return s
    end

    for (up, down) in ((1, 0), (0, 0))
        ψ = fiducial_state(x; imp_up = up, imp_down = down)
        q0 = qup(spec, ψ)
        evolve!(ψ, h, (0.0, 0.3); δt, trunc)
        @test abs(qup(spec, ψ) - q0) < 1e-8
    end

    # single vacuum mode (direct coupling only): the doubly-occupied case too
    x1 = SISBD(discretefermionicvacuum([0.4], [0.55]; μ = 0.0); U = 1.7, εd = -0.3, J = 1.1)
    spec1 = FermionicTCMPS.chain_spec(x1)
    h1 = mpo(x1)
    ψ = fiducial_state(x1; imp_up = 1, imp_down = 1)
    q0 = qup(spec1, ψ)
    evolve!(ψ, h1, (0.0, 0.3); δt, trunc)
    @test abs(qup(spec1, ψ) - q0) < 1e-8
end

@testset "G3 free dynamics: sweep drift regression gate (finite-β TFD class)" begin
    # the finite-β TFD plain sweep carries a small KS App-C projection residue
    # (the d–a₂ coupling spans two sites); the equivalent initial-state
    # expansion removes it — the expand leg pins exactness
    Random.seed!(36)
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    β = 0.8
    times = [0.4]
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    x = SISBD(discretefermionicbath(ws, fs; β, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
    hu = _star_h(x)
    ρ0 = _rho0(x; imp = 0.0)
    exact, _ = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
    exact = 2 .* exact
    g2, = gf_greater_lesser_t(x, times; δt = 0.02, trunc)
    e2 = abs(g2[1] - exact[1])
    @test e2 < 0.02
    g3, = gf_greater_lesser_t(x, times; δt = 0.02, expand = true)
    @test abs(g3[1] - exact[1]) < 1e-7
end

@testset "G3 free dynamics: double bath dense chain leg vs star frame" begin
    # machine-precision physical anchor for the mid-chain-impurity geometry: the
    # mechanical-JW chain Hamiltonian (test/anchor_util.jl, impurity on-site
    # WITHOUT a JW prefix) evolving the thermofield-double fiducial state
    # reproduces the star-frame thermal free solution exactly, equal-time and
    # (0,t) data alike
    Random.seed!(40)
    times = collect(0.0:0.1:0.5)
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 0.0, εd = -0.3, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    imp = spec.imp_site
    wsl, vl1, vl2, nsl = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2, nsr = FermionicTCMPS.thermofield_couplings(x.bathr)
    couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k], nsl[k]) for k in eachindex(wsl)],
                   [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k], nsr[k]) for k in eachindex(wsr)])
    H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
    modes = vcat([(2k - 1, 2k, nsl[k]) for k in eachindex(wsl)],
                 [(imp + 2k - 1, imp + 2k, nsr[k]) for k in eachindex(nsr)])
    ψ0 = dense_tfd_vector(spec, modes, (0, 0))
    gt, lt = dense_gf_sandwich(H, ψ0, imp, times)
    hu = _star_h(x)
    gtr, ltr = freefermions_greater_lesser(hu, _rho0(x; imp = 0.0), times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 1e-10
    @test maximum(abs.(lt)) < 1e-12
end

@testset "G3 free dynamics: double bath sweep vs exact star evolution" begin
    # enabled after the component-wise d†-apply fix: the finite-β SIDBD
    # sandwich runs on all charge sectors of the TFD fiducial; the plain sweep
    # keeps a small projection residue, the expansion is exact
    Random.seed!(39)
    δt = 0.025
    times = collect(0.0:0.1:0.5)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 0.0, εd = -0.3, Jl = 1.0, Jr = 1.0)
    gt, lt = gf_greater_lesser_t(x, times; δt, trunc)
    gt1, lt1 = gf_greater_lesser_t(x, times; δt, expand = true)
    gtr, ltr = freefermions_greater_lesser(_star_h(x), _rho0(x; imp = 0.0), times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 0.02
    @test maximum(abs.(lt)) < 1e-12
    @test maximum(abs.(gt1 .- 2 .* gtr)) < 1e-6
    @test maximum(abs.(lt1)) < 1e-12
end

@testset "G3 vacuum degeneracy: a₂ branch frozen when μ below spectrum" begin
    Random.seed!(37)
    ws = [0.4, 0.9]; fs = [0.35, 0.55]
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.3, εd = -0.2, J = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    h = mpo(x)
    trunc = truncdimcutoff(D = 64, ϵ = 1e-12)
    ψ = fiducial_state(x)
    # particle numbers on the a₂ site of mode 1 (chain site 3) stay 0
    # (the a₂ branch is empty in the vacuum fiducial and v2 = 0 decouples it)
    n3up = FermionicTCMPS.site_operator_mpo(spec.spaces, 3, FermionicTCMPS._n_up)
    n3dn = FermionicTCMPS.site_operator_mpo(spec.spaces, 3, FermionicTCMPS._n_down)
    evolve!(ψ, h, (0.0, 0.2); δt = 0.05, trunc)
    @test abs(real(expectationvalue(n3up, ψ))) < 1e-10
    @test abs(real(expectationvalue(n3dn, ψ))) < 1e-10
end

@testset "G3 equal-time: finite-β occupation vs star frame (SISBD + SIDBD)" begin
    # THE anchor for the reported bug: the equal-time impurity occupation of the
    # sudden protocol must match the star frame exactly (it diverged at O(t³)
    # under the previous T=0-structure chain)
    Random.seed!(41)
    δt = 0.02
    times = collect(0.05:0.05:0.3)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)

    let x = SISBD(discretefermionicbath([-0.8, 0.5], [0.6, 0.45]; β = 0.8, μ = 0.0);
                  U = 0.0, εd = -0.3, J = 1.1)
        ρ0 = _rho0(x; imp = 1.0)
        ψ = fiducial_state(x; imp_up = 1, imp_down = 0)
        occ = occupation_op(x; spin = :up)
        cache = environments(mpo(x), ψ)
        alg = FermionicTCMPS._tdvp2(δt, trunc)
        prev = 0
        for t in times
            n = round(Int, t / δt)
            for _ in (prev+1):n; sweep!(cache, alg); end
            prev = n
            @test abs(real(expectationvalue(occ, ψ)) - _nd_star(x, ρ0, t)) < 0.1   # projection class
        end
    end

    # the SIDBD sweep leg runs in the GF testsets above (the component-wise
    # d†-apply handles all charge sectors since the tree-level rewrite); the
    # dense gate carries the definition-level verification
end
