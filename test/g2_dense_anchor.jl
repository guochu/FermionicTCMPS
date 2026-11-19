include("anchor_util.jl")

@testset "G2 dense anchor: FSM MPO vs mechanical JW" begin
    Random.seed!(20261119)

    # ---------- 1. SISBD, vacuum (μ = 0, spectrum straddling zero) ----------
    ws = [-0.9, 0.5]
    fs = [0.35, 0.55]
    x = SISBD(discretefermionicvacuum(ws, fs; μ = 0.0); U = 1.7, εd = -0.3, J = 1.1)
    spec = FermionicTCMPS.chain_spec(x)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    @test v2[1] ≈ sqrt(fs[1]) && iszero(v2[2])   # ω < μ filled: n = 1; ω > μ: n = 0
    @test iszero(v1[1]) && v1[2] ≈ sqrt(fs[2])
    couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)]
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    @test _anchor_model(x, H) < 1e-12

    # ---------- 2. SISBD, finite β ----------
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 2.3, εd = 0.4, J = 0.9)
    spec = FermionicTCMPS.chain_spec(x)
    _, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    @test all(0 .< v1 .< 1) && all(v2 .> 0)
    couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k]) for k in eachindex(ws)]
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    @test _anchor_model(x, H) < 1e-12

    # ---------- 3. SIDBD, vacuum ----------
    xl = discretefermionicvacuum([-0.6], [0.5]; μ = 0.0)
    xr = discretefermionicvacuum([0.8], [0.7]; μ = 0.0)
    x = SIDBD(xl, xr; U = 1.2, εd = 0.1, Jl = 1.3, Jr = 0.8)
    spec = FermionicTCMPS.chain_spec(x)
    imp = spec.imp_site
    @test spec.kinds[imp] === :particle
    wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
    couples = vcat(
        [(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in 1:length(vl1)],
        [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in 1:length(vr1)],
    )
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    @test _anchor_model(x, H) < 1e-12

    # ---------- 4. SIDBD, finite β (asymmetric spectra) ----------
    xl = discretefermionicbath([-0.8], [0.6]; β = 1.2, μ = 0.0)
    xr = discretefermionicbath([0.7], [0.35]; β = 0.6, μ = 0.1)
    x = SIDBD(xl, xr; U = 0.8, εd = -0.2, Jl = 1.0, Jr = 1.0)
    spec = FermionicTCMPS.chain_spec(x)
    imp = spec.imp_site
    wsl, vl1, vl2 = FermionicTCMPS.thermofield_couplings(x.bathl)
    wsr, vr1, vr2 = FermionicTCMPS.thermofield_couplings(x.bathr)
    couples = vcat(
        [(2k - 1, 2k, wsl[k], x.Jl * vl1[k], x.Jl * vl2[k]) for k in 1:length(vl1)],
        [(imp + 2k - 1, imp + 2k, wsr[k], x.Jr * vr1[k], x.Jr * vr2[k]) for k in 1:length(vr1)],
    )
    H = dense_hamiltonian(copy(spec.kinds), spec.imp_site, couples; εd = x.εd, U = x.U)
    @test _anchor_model(x, H) < 1e-12
end

@testset "G2 dense anchor: current operator MPO" begin
    Random.seed!(77)
    ws = [-0.9, 0.5]
    fs = [0.35, 0.55]
    x = SISBD(discretefermionicbath(ws, fs; β = 0.8, μ = 0.0); U = 1.7, εd = -0.3, J = 1.1)
    spec = FermionicTCMPS.chain_spec(x)
    wsf, v1, v2 = FermionicTCMPS.thermofield_couplings(x.bath)
    hoppings = vcat([(2k, x.J * v1[k]) for k in eachindex(v1)],
                    [(2k + 1, x.J * v2[k]) for k in eachindex(v2)])
    L = length(spec.kinds)

    occs = _all_occs(L)
    # dense current operator: Σ_bath V' (c_site† c_imp + h.c.) per spin, both branches
    prefixes = Vector{Any}(undef, L + 1)
    prefixes[1] = ones(1, 1)
    for s in 1:L
        prefixes[s+1] = kron(prefixes[s], _LOC.par)
    end
    imp = spec.imp_site
    imp_ops = (_kron3(prefixes[imp], _LOC.cup, _suffix_id(L - imp)),
               _kron3(prefixes[imp], _LOC.cdn, _suffix_id(L - imp)))
    Jd = zeros(4^L, 4^L)
    for (s, vc) in hoppings
        site_ops = spec.kinds[s] === :hole ?
                   (_kron3(prefixes[s], _LOC.bup, _suffix_id(L - s)),
                    _kron3(prefixes[s], _LOC.bdn, _suffix_id(L - s))) :
                   (_kron3(prefixes[s], _LOC.cup, _suffix_id(L - s)),
                    _kron3(prefixes[s], _LOC.cdn, _suffix_id(L - s)))
        Jd .+= (imp_ops[1]' * site_ops[1] .+ site_ops[1]' * imp_ops[1]) .* vc
        Jd .+= (imp_ops[2]' * site_ops[2] .+ site_ops[2]' * imp_ops[2]) .* vc
    end
    jop = current_op(x)
    err = 0.0
    for _ in 1:300
        v, w = rand(occs), rand(occs)
        _qtotal(spec, v) == _qtotal(spec, w) || continue
        ψv = _basis_state(spec, v)
        ψw = _basis_state(spec, w)
        err = max(err, abs(expectation(ψv, jop, ψw) - Jd[_flat(spec, v), _flat(spec, w)]))
    end
    @test err < 1e-12
end
