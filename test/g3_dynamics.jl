# ==============================================================================
# G3 dynamics anchors: TDVP real-time evolution of the thermo-field chain.
#
# Reference side: exact free-fermion evolution of the STAR model with the density
# matrix ρ₀ = fixed impurity occupations ⊗ thermofield-fiducial bath occupations
# (diag(n_F(ω_k))), which is the star-frame image of the separable chain state.
# The a₂ branch PHT is again exercised through the filled modes below μ.
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

@testset "G3 free dynamics: greater/lesser GF vs exact star evolution" begin
    Random.seed!(31)
    ws = [-0.9, 0.5]; fs = [0.35, 0.55]
    β = 0.8; μ = 0.0; δt = 0.025
    times = collect(0.0:0.1:0.5)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)

    # impurity empty: greater nontrivial, lesser zero
    x = SISBD(discretefermionicbath(ws, fs; β, μ); U = 0.0, εd = -0.3, J = 1.1)
    gt, lt = gf_greater_lesser_t(x, times; δt, trunc)
    hu = _star_h(x)
    ρ0 = _rho0(x; imp = 0.0)
    gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 1e-6        # spin-summed
    @test maximum(abs.(lt)) < 1e-12                    # d annihilates the empty impurity

    # impurity up-occupied: the greater runs on the down channel only (G^>(0) = −i),
    # and the lesser of the occupied spin is analytically i·e^{−iεd·t}
    # (the annihilated state is an eigenstate of the chain)
    gt2, lt2 = gf_greater_lesser_t(x, times; imp_up = 1, δt, trunc)
    @test abs(gt2[1] + im) < 1e-8
    @test maximum(abs.(lt2 .- im .* exp.(-im .* x.εd .* times))) < 1e-8
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
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    @test all(iszero, v2) && all(!iszero, v1)
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

    ψ = separable_state(x; imp_up = 1, imp_down = 0)
    cache = environments(h, ψ)
    alg = FermionicTCMPS._tdvp2(δt, trunc)
    for n in 1:10
        sweep!(cache, alg)
        t = n * δt
        (n % 2 == 0) && (@test abs(real(expectationvalue(occ, ψ)) - _nd_star(x, ρ0, t)) < 1e-5)
    end

    # empty impurity: the separable state is an exact eigenstate — frozen
    ψ0 = separable_state(x)
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

    # Q↑ = n_imp↑ + Σ_{a₁} n↑ + Σ_{a₂} (holes↑) is conserved (the −2·Σn_h offset
#     from measuring −n_p at a₂ sites would not be)
    function qup(ψ)
        s = real(expectationvalue(FermionicTCMPS.site_operator_mpo(spec.spaces, spec.imp_site, FermionicTCMPS._n_up), ψ))
        for k in 2:length(spec.kinds)
            # conserved charge at a hole site = hole number, in hole-sector labels
            tbl = spec.kinds[k] === :hole ? FermionicTCMPS._nh_up : FermionicTCMPS._n_up
            s += real(expectationvalue(FermionicTCMPS.site_operator_mpo(spec.spaces, k, tbl), ψ))
        end
        return s
    end

    for (up, down) in ((1, 0), (0, 0))
        ψ = separable_state(x; imp_up = up, imp_down = down)
        q0 = qup(ψ)
        evolve!(ψ, h, (0.0, 0.3); δt, trunc)
        @test abs(qup(ψ) - q0) < 1e-8
    end

    # single vacuum mode (direct coupling only): the doubly-occupied case too
    x1 = SISBD(discretefermionicvacuum([0.4], [0.55]; μ = 0.0); U = 1.7, εd = -0.3, J = 1.1)
    spec1 = FermionicTCMPS.chain_spec(x1)
    h1 = mpo(x1)
    ψ = separable_state(x1; imp_up = 1, imp_down = 1)
    function qup1(ψ)
        s = real(expectationvalue(FermionicTCMPS.site_operator_mpo(spec1.spaces, spec1.imp_site, FermionicTCMPS._n_up), ψ))
        for k in 2:length(spec1.kinds)
            tbl = spec1.kinds[k] === :hole ? FermionicTCMPS._nh_up : FermionicTCMPS._n_up
            s += real(expectationvalue(FermionicTCMPS.site_operator_mpo(spec1.spaces, k, tbl), ψ))
        end
        return s
    end
    q0 = qup1(ψ)
    evolve!(ψ, h1, (0.0, 0.3); δt, trunc)
    @test abs(qup1(ψ) - q0) < 1e-8
end

@testset "G3 free dynamics: δt² convergence of the greater GF" begin
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
    g1, = gf_greater_lesser_t(x, times; δt = 0.04, trunc)
    g2, = gf_greater_lesser_t(x, times; δt = 0.02, trunc)
    e1, e2 = abs(g1[1] - exact[1]), abs(g2[1] - exact[1])
    # two-site TDVP is second order: halving δt quarters the error
    @test e2 < e1 / 2.5
    @test e2 < 5e-6
end

# exact (dense mechanical-JW) impurity greater GF of the chain-frame fiducial state:
# G^>(t) = −i·e^{iE₀t}·⟨χ|e^{−iHt}|χ⟩ summed over the spin insertions χ = d†_σ|ψ₀⟩
function _dense_greater(x, couples, times)
    spec = FermionicTCMPS.chain_spec(x)
    L = length(spec.kinds); imp = spec.imp_site
    H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
    flat(occs) = 1 + sum((occs[s][1] + 2 * occs[s][2]) * 4^(L - s) for s in 1:L)
    occs0 = ntuple(i -> i == imp ? (0, 0) : (spec.kinds[i] === :hole ? (1, 1) : (0, 0)), L)
    ψ0 = zeros(ComplexF64, 4^L); ψ0[flat(occs0)] = 1.0
    E0 = real(H[flat(occs0), flat(occs0)])
    F = eigen(Hermitian(H))
    gt = zeros(ComplexF64, length(times))
    for (tbl, step, up) in ((FermionicTCMPS._adag_up, 4^(L - imp), true),
                            (FermionicTCMPS._adag_down, 2 * 4^(L - imp), false))
        χ = zeros(ComplexF64, 4^L); found = false
        for i in eachindex(ψ0)
            iszero(ψ0[i]) && continue
            d = (i - 1) ÷ 4^(L - imp) % 4
            q_in = (d == 1 || d == 3 ? 1 : 0, d >= 2 ? 1 : 0)
            q_out = up ? (1, q_in[2]) : (q_in[1], 1)
            (q_out[1] <= 1 && q_out[2] <= 1) || continue
            amp = FermionicTCMPS._amplitude(tbl, q_out, q_in)
            amp === nothing && continue
            χ[i + step] += amp * ψ0[i]; found = true
        end
        found || continue
        g = F.vectors' * χ
        for (j, t) in enumerate(times)
            gt[j] += -im * exp(im * E0 * t) *
                     sum(abs2(g[k]) * exp(-im * F.values[k] * t) for k in eachindex(g))
        end
    end
    return gt
end

@testset "G3 free dynamics: double bath dense chain leg vs star frame" begin
    # machine-precision physical anchor for the mid-chain-impurity geometry: the
    # mechanical-JW chain Hamiltonian (test/anchor_util.jl, impurity on-site
    # WITHOUT a JW prefix) evolving the fiducial product state reproduces the
    # star-frame thermal free solution exactly
    Random.seed!(40)
    times = collect(0.0:0.1:0.5)
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 0.0, εd = -0.3, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    imp = spec.imp_site
    wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
    couples = vcat([(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in eachindex(wsl)],
                   [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in eachindex(wsr)])
    gt = _dense_greater(x, couples, times)
    gtr, = freefermions_greater_lesser(_star_h(x), _rho0(x; imp = 0.0), times, 1, 1)
    @test maximum(abs.(gt .- 2 .* gtr)) < 1e-10
end

@testset "G3 free dynamics: double bath sweep vs exact star evolution" begin
    # the mid-chain impurity exercises both channel families (coefficients at the
    # opening end on the left bath, at the closing end on the right one); the
    # star-frame free solution is the physical anchor for that geometry
    Random.seed!(39)
    δt = 0.025
    times = collect(0.0:0.1:0.5)
    trunc = truncdimcutoff(D = 256, ϵ = 1e-12)
    xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
    xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
    x = SIDBD(xl, xr; U = 0.0, εd = -0.3, Jl = 1.0, Jr = 1.0)
    gt, lt = gf_greater_lesser_t(x, times; δt, trunc)
    hu = _star_h(x)
    ρ0 = _rho0(x; imp = 0.0)
    gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
    # the chain-frame conventions are exact for the double bath (see the dense
    # chain-leg testset above); what remains here is the sweep pipeline's
    # intrinsic O(δt) splitting error, which does NOT cancel for the mid-chain
    # geometry (0.0024 at δt = 0.025, halving with δt) — unlike the single-bath
    # cases above, where the TDVP manifold contains the dynamics
    @test maximum(abs.(gt .- 2 .* gtr)) < 5e-3        # spin-summed
    @test maximum(abs.(lt)) < 1e-12                    # d annihilates the empty impurity
end

@testset "G3 vacuum degeneracy: a₂ branch frozen when μ below spectrum" begin
    Random.seed!(37)
    ws = [0.4, 0.9]; fs = [0.35, 0.55]
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.3, εd = -0.2, J = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    h = mpo(x)
    trunc = truncdimcutoff(D = 64, ϵ = 1e-12)
    ψ = separable_state(x)
    # hole numbers on the a₂ site of mode 1 (chain site 3) stay 0
    # (the a₂ branch carries no holes in the fiducial state and v2 = 0 decouples it)
    nh_up3 = FermionicTCMPS.site_operator_mpo(spec.spaces, 3, FermionicTCMPS._nh_up)
    nh_dn3 = FermionicTCMPS.site_operator_mpo(spec.spaces, 3, FermionicTCMPS._nh_down)
    evolve!(ψ, h, (0.0, 0.2); δt = 0.05, trunc)
    @test abs(real(expectationvalue(nh_up3, ψ))) < 1e-10
    @test abs(real(expectationvalue(nh_dn3, ψ))) < 1e-10
end
