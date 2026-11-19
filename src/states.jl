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
    return prodmps(T, spec.spaces, sectors)
end
