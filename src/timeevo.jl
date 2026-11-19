# ==============================================================================
# Real-time evolution driven by two-site TDVP.
# ==============================================================================

function _tdvp2(δt::Real, trunc::DMRG.TruncationDimCutoff)
    return TDVP2(stepsize = -im * δt, ishermitian = true, trunc = trunc)
end

"""
    evolve!(ψ, h; δt, nsteps, trunc=truncdimcutoff(D=256, ϵ=1e-10))
    evolve!(ψ, h, tspan; δt, trunc)

Evolve `ψ` under the MPO Hamiltonian `h` in real time (|ψ⟩ → exp(−iHt)|ψ⟩) with
two-site TDVP, one full left+right sweep per step, truncating with `trunc`.
`tspan = (t₀, t₁)` evolves for `(t₁ − t₀)/δt` steps; the MPS is modified in place
and returned.
"""
function evolve!(ψ, h; δt::Real, nsteps::Integer,
                 trunc::DMRG.TruncationDimCutoff = truncdimcutoff(D = 256, ϵ = 1.0e-10))
    alg = _tdvp2(δt, trunc)
    cache = environments(h, ψ)
    for _ in 1:nsteps
        sweep!(cache, alg)
    end
    return ψ
end

function evolve!(ψ, h, tspan::Tuple{<:Real,<:Real}; δt::Real, kwargs...)
    nsteps = round(Int, (tspan[2] - tspan[1]) / δt)
    return evolve!(ψ, h; δt = δt, nsteps = nsteps, kwargs...)
end
