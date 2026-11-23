# Decisive: is the sparse warm-up accurate at t = δt? occupation vs dense.
using Random, TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath
using LinearAlgebra
const F = FermionicTCMPS
include(joinpath(pwd(), "test", "anchor_util.jl"))
Random.seed!(31)
ws = [-0.9, 0.5]; fs = [0.35, 0.55]
β = 0.8; δt = 0.025
x = SISBD(discretefermionicbath(ws, fs; β, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
spec = F.chain_spec(x)
L = length(spec.kinds)
imp = spec.imp_site
_, v1, v2, ns = F.thermofield_couplings(x.bath)
couples = [(2k, 2k + 1, ws[k], x.J * v1[k], x.J * v2[k], ns[k]) for k in eachindex(ws)]
H = dense_hamiltonian(copy(spec.kinds), imp, couples; εd = x.εd, U = x.U)
Fh = eigen(Hermitian(H))
modes = [(c[1], c[2], c[6]) for c in couples]
ψ0 = dense_tfd_vector(spec, modes, (0, 0))
g0 = Fh.vectors' * ψ0
ψδ_dense = Fh.vectors * (exp.(-im .* Fh.values .* δt) .* g0)

h = mpo(x)
δw = δt / 10
Wm = DMRG.MPO(DMRG.timeevompo(h, -im * δw, DMRG.WII(tol = 1.0e-12, maxiter = 10_000)))
wt = truncdimcutoff(D = 256, ϵ = 1.0e-12, add_back = 0)
occ = occupation_op(x; spin = :up)
nddiag = diagm([((i - 1) ÷ 4^(L - imp)) % 4 in (1, 3) ? 1.0 : 0.0 for i in 1:4^L])
nd_dense = real(sum(conj(ψδ_dense) .* (nddiag * ψδ_dense)))

function warmocc(h, x, wt, occ, nsteps_warm)
    ψ = fiducial_state(x)
    c = ψ.comps[1]
    for _ in 1:nsteps_warm
        c = Wm * c
        canonicalize!(c, alg = Orthogonalize(alg = DMRG.SVD(), trunc = wt))
    end
    env = DMRG.environments(h, DMRG.ExactMPS(c))
    return real(real(DMRG.expectation(env.state, occ, env.state))), [dim(space(c.data[b], 3)) for b in 1:5]
end
Wm2 = DMRG.MPO(DMRG.timeevompo(h, -im * δw, DMRG.WII(tol = 1.0e-12, maxiter = 10_000)))
n_warm, bonds_w = warmocc(h, x, wt, occ, 10)

ψp = fiducial_state(x)
cp = ψp.comps[1]
alg = F._tdvp2(δt, truncdimcutoff(D = 256, ϵ = 1e-10, add_back = 10))
cachep = DMRG.environments(h, cp)
DMRG.sweep!(cachep, alg)
n_plain = real(DMRG.expectationvalue(occ, cachep.state))

Wbig = DMRG.MPO(DMRG.timeevompo(h, -im * δt, DMRG.WII(tol = 1.0e-12, maxiter = 10_000)))
ψb = fiducial_state(x)
cb = Wbig * ψb.comps[1]
canonicalize!(cb, alg = Orthogonalize(alg = DMRG.SVD(), trunc = wt))
envb = DMRG.environments(h, DMRG.ExactMPS(cb))
n_big = real(DMRG.expectation(envb.state, occ, envb.state))

println("⟨n_up(δt)⟩: dense = ", round(nd_dense, digits = 8),
        "  10×WII δt/10 = ", round(n_warm, digits = 8), "  bonds = ", bonds_w,
        "  1×WII δt = ", round(n_big, digits = 8),
        "  plain-TDVP = ", round(n_plain, digits = 8))
