# ==============================================================================
# Initial states.
# ==============================================================================

"""
    separable_state(x::Union{SISBD,SIDBD}; imp_up=0, imp_down=0, T=ComplexF64)

Separable product state of the thermo-field chain: the impurity orbital carries the
given up/down occupations and every bath branch is in its fiducial state
(a₁ sites empty, a₂ sites fully occupied = zero holes), i.e. |∅₁⟩⊗|F₂⟩.
"""
function separable_state(x::Union{SISBD,SIDBD}; imp_up::Int = 0, imp_down::Int = 0,
                         T::Type{<:Number} = ComplexF64)
    (imp_up in (0, 1) && imp_down in (0, 1)) || throw(ArgumentError("impurity occupations must be 0 or 1"))
    spec = chain_spec(x)
    sectors = [_sector(k == spec.imp_site ? (imp_up, imp_down) : (0, 0)) for k in 1:length(spec.kinds)]
    return _widen_svectors(prodmps(T, spec.spaces, sectors))
end

# DMRG's TDVP stores plain (1,1) tensors into `mps.s`, while `prodmps` leaves the
# container typed to DiagonalTensorMap; widen it so the assignment type-checks.
function _widen_svectors(ψ::DMRG.MPS)
    B = DMRG.MPSBondTensor{spacetype(ψ[1])}
    raw = getfield(ψ, :s)
    sv = Vector{Union{Missing, B}}(missing, length(raw))
    for i in eachindex(raw)
        raw[i] isa B && (sv[i] = raw[i])
    end
    return DMRG.MPS(ψ.data, sv)
end
