# WI/WII warm-up vs plain vs expand: finite-β K=2 accuracy comparison
using Random, TensorKit, DMRG, FermionicTCMPS
import ImpurityModelBase
using ImpurityModelBase: discretefermionicbath, freefermions_greater_lesser
using LinearAlgebra
const F = FermionicTCMPS
Random.seed!(31)
ws = [-0.9, 0.5]; fs = [0.35, 0.55]
β = 0.8; δt = 0.025
times = collect(0.0:0.1:0.5)
x = SISBD(discretefermionicbath(ws, fs; β, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
hu = zeros(3,3); hu[1,1] = -0.3
hu[1,2] = hu[2,1] = 1.1*sqrt(0.35); hu[2,2] = -0.9
hu[1,3] = hu[3,1] = 1.1*sqrt(0.55); hu[3,3] = 0.5
ρ0 = diagm([0.0, 1/(exp(-0.8*(-0.9))+1), 1/(exp(-0.8*(0.5))+1)])
gtr, ltr = freefermions_greater_lesser(hu, ρ0, times, 1, 1)
exact = 2 .* gtr

gt0, lt0 = gf_greater_lesser_t(x, times; δt)
println("plain sweep      : max |gt−2gtr| = ", round(maximum(abs.(gt0 .- exact)), digits = 6),
        "  per-t: ", [round(abs(gt0[i]-exact[i]), digits=4) for i in eachindex(times)])
gt1, lt1 = gf_greater_lesser_t(x, times; δt, warmup = :WI)
println("WI warm-up       : max |gt−2gtr| = ", round(maximum(abs.(gt1 .- exact)), digits = 6),
        "  per-t: ", [round(abs(gt1[i]-exact[i]), digits=4) for i in eachindex(times)])
gt2, lt2 = gf_greater_lesser_t(x, times; δt, warmup = :WII)
println("WII warm-up      : max |gt−2gtr| = ", round(maximum(abs.(gt2 .- exact)), digits = 6),
        "  per-t: ", [round(abs(gt2[i]-exact[i]), digits=4) for i in eachindex(times)])
gt3, lt3 = gf_greater_lesser_t(x, times; δt, expand = true)
println("expand           : max |gt−2gtr| = ", round(maximum(abs.(gt3 .- exact)), digits = 8))
