# Design anchor for the TFD equal-time fix (dense, single-particle level).
#
# Exact Bogoliubov image of the thermal star bath (one spin sector; U = 0):
#   star  h_s = [[εd, V], [V, ω]],          Γ_s = diag(imp, n_F(ω))
#   chain h_A = [[εd, v1, v2],
#                [v1, ω(1−n), ω·c],
#                [v2, ω·c, ω·n    ]],       Γ_A = diag(imp, n, n)
# with v1 = V√(1−n), v2 = V√n, c = √(n(1−n)).
# span{e_d, w} with w = √(1−n)e₁ + √n e₂ is h_A-invariant and carries h_s; PΓ_AP = Γ_s
# and (1−P)Γ_AP = 0, so equal-time AND (0,t) impurity data are reproduced exactly.
# The OLD chain (KS Eq. 9 T=0 structure: unweighted on-sites, no intra-mode pairing)
# fails the equal-time comparison at O(t³) — the reported bug.
using LinearAlgebra

const εd, J, β, μ = -0.3, 1.1, 0.8, 0.0
nF(ω) = 1 / (exp(β * (ω - μ)) + 1)

"max |star − exact-image-chain| over times for equal-time n(t) and (0,t) greater"
function check_exactimage(ws, fs, imp; times = 0.0:0.05:0.8)
    K = length(ws)
    ns = nF.(ws)
    Vs = J .* fs
    # star: imp + K diagonal bath modes per spin
    h_s = zeros(1 + K, 1 + K)
    h_s[1, 1] = εd
    for k in 1:K
        h_s[1, 1+k] = Vs[k]; h_s[1+k, 1] = Vs[k]; h_s[1+k, 1+k] = ws[k]
    end
    Γ_s = diagm(vcat([Float64(imp)], ns))

    # exact-image chain: imp + interleaved (a₁, a₂)
    L = 1 + 2K
    h_A = zeros(L, L)
    h_A[1, 1] = εd
    for k in 1:K
        s1, s2 = 2k, 2k + 1
        n, c = ns[k], sqrt(ns[k] * (1 - ns[k]))
        h_A[1, s1] = Vs[k] * sqrt(1 - n); h_A[s1, 1] = h_A[1, s1]
        h_A[1, s2] = Vs[k] * sqrt(n);     h_A[s2, 1] = h_A[1, s2]
        h_A[s1, s1] = ws[k] * (1 - n)     # weighted particle branch
        h_A[s2, s2] = ws[k] * n           # weighted branch, particle-labelled
        h_A[s1, s2] = ws[k] * c; h_A[s2, s1] = h_A[s1, s2]  # intra-mode pairing
    end
    Γ_A = diagm(vcat([Float64(imp)], repeat(ns; inner = 2)))

    F_s = eigen(h_s); F_A = eigen(h_A)
    U_s(t) = F_s.vectors * diagm(exp.(-im .* F_s.values .* t)) * F_s.vectors'
    U_A(t) = F_A.vectors * diagm(exp.(-im .* F_A.values .* t)) * F_A.vectors'

    en, eq = 0.0, 0.0
    for t in times
        en = max(en, abs(real((U_s(t) * Γ_s * U_s(t)')[1, 1]) -
                         real((U_A(t) * Γ_A * U_A(t)')[1, 1])))
        eq = max(eq, abs(real((U_s(t) * (I - Γ_s))[1, 1]) -
                         real((U_A(t) * (I - Γ_A))[1, 1])))
    end
    return en, eq
end

"OLD chain (KS Eq. 9 carried to finite T): unweighted on-sites, no pairing"
function check_old(ws, fs, imp; times = 0.0:0.05:0.8)
    K = length(ws)
    ns = nF.(ws)
    Vs = J .* fs
    h_s = zeros(1 + K, 1 + K); h_s[1, 1] = εd
    for k in 1:K
        h_s[1, 1+k] = Vs[k]; h_s[1+k, 1] = Vs[k]; h_s[1+k, 1+k] = ws[k]
    end
    Γ_s = diagm(vcat([Float64(imp)], ns))

    L = 1 + 2K
    h_old = zeros(L, L); h_old[1, 1] = εd
    for k in 1:K
        s1, s2 = 2k, 2k + 1
        n = ns[k]
        h_old[1, s1] = Vs[k] * sqrt(1 - n); h_old[s1, 1] = h_old[1, s1]
        h_old[1, s2] = Vs[k] * sqrt(n);     h_old[s2, 1] = h_old[1, s2]
        h_old[s1, s1] = ws[k]
        h_old[s2, s2] = -ws[k]
    end
    Γ_old = diagm(vcat([Float64(imp)], zeros(K), ones(K)))   # a₂ particles full

    Fs, Fo = eigen(h_s), eigen(h_old)
    Us(t) = Fs.vectors * diagm(exp.(-im .* Fs.values .* t)) * Fs.vectors'
    Uo(t) = Fo.vectors * diagm(exp.(-im .* Fo.values .* t)) * Fo.vectors'
    return maximum(abs(real((Us(t) * Γ_s * Us(t)')[1, 1]) -
                       real((Uo(t) * Γ_old * Uo(t)')[1, 1])) for t in times)
end

let
    println("config                 maxerr n(t)   maxerr greater(t)")
    for (name, ws, fs, imp) in (
        ("K=1  ω=0.4    imp=1", [0.4], [0.55], 1),
        ("K=1  ω=0.4    imp=0", [0.4], [0.55], 0),
        ("K=1  ω=-0.5   imp=1", [-0.5], [0.6], 1),
        ("K=2  -0.8/0.5 imp=1", [-0.8, 0.5], [0.6, 0.45], 1),
    )
        en, eq = check_exactimage(ws, fs, imp)
        println("$name   $(round(en, digits=12))    $(round(eq, digits=12))")
    end
    println()
    println("OLD chain equal-time deviation (the reported bug):")
    e = check_old([0.4], [0.55], 1)
    println("K=1  ω=0.4    imp=1   ", round(e, digits=6))
end
