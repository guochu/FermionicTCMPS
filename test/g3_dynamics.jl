# ==============================================================================
# G3 dynamics anchors: TDVP real-time evolution of the thermo-field chain.
#
# Reference side: exact free-fermion evolution of the STAR model with the density
# matrix ρ₀ = fixed impurity occupations ⊗ thermofield-fiducial bath occupations
# (diag(n_F(ω_k))), which is the star-frame image of the separable chain state.
# The a₂ branch PHT is again exercised through the filled modes below μ.
# ==============================================================================

using LinearAlgebra: eigen, diagm, tr

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
