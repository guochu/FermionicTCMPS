# ==============================================================================
# Real-time evolution driven by two-site TDVP.
#
# Projection note (seeding): with coupling terms spanning beyond the current
# sweep pair, the two-site projection zeroes the corresponding amplitudes when
# the bond spaces carry no directions for them yet — from a sharp product state
# the dynamics then gets stuck (Kohn–Santoro 2107.02807, Appendix C, "TDVP
# beyond nearest neighbor"). The prescribed cure is a minimum bond dimension:
# keep discarded singular directions via the truncation's add_back, so the
# tangent space retains room for the not-yet-seeded amplitudes. The pipeline
# defaults therefore use truncdimcutoff(D, ϵ, add_back = MIND) with MIND = 10
# (the value validated in the KS appendix).
# ==============================================================================

const MIND = 10   # minimum bond dimension kept for the projection directions

# Truncation for equivalently expanded initial states: the expansion directions
# acquire physical weight at O(δt·v) from the first step, so neither the KS MinD
# crutch nor a large bond cap is needed — discarded (zero-amplitude) directions
# are simply truncated away.
const TRUNC_EXPANDED = truncdimcutoff(D = 10, ϵ = 1.0e-10, add_back = 0)

# Bond-dimension budget of the equivalent initial-state expansion: two slots per
# admitted far-pair direction on top of the base channel's intrinsic charge
# branches. The directions are admitted best-first over their coupling strength,
# so the cap bounds the expansion polynomially (the states themselves stay exact).
# Slot floor of the seeding breadth: at least EXPANDED_MIN_DIR direction slots
# (± pairs over the direction descriptors) are admitted even when the bond
# budget D is exhausted by the base channel's intrinsic charge branches (large
# K), so the tangent space keeps its far-pair generators at O(K) extra slots.
const EXPANDED_MIN_DIR = 64

const EXPANDED_DBUDGET = 160

function _tdvp2(δt::Real, trunc::DMRG.TruncationDimCutoff)
    return TDVP2(stepsize = -im * δt, ishermitian = true, trunc = trunc)
end

"""
    evolve!(ψ, h; δt, nsteps, trunc=truncdimcutoff(D=256, ϵ=1e-10, add_back=MIND))
    evolve!(ψ, h, tspan; δt, trunc)

Evolve `ψ` under the MPO Hamiltonian `h` in real time (|ψ⟩ → exp(−iHt)|ψ⟩) with
two-site TDVP, one full left+right sweep per step, truncating with `trunc`.
`tspan = (t₀, t₁)` evolves for `(t₁ − t₀)/δt` steps; the MPS is modified in place
and returned.
"""
function evolve!(ψ, h; δt::Real, nsteps::Integer,
                 trunc::DMRG.TruncationDimCutoff = truncdimcutoff(D = 256, ϵ = 1.0e-10, add_back = MIND))
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
