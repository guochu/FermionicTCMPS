# generate separable state


const ThermofieldOneBathDIM{C} = OneImpurityOneBathDIM{C} where {C<:ThermofieldConfiguration}
const VacuumOneBathDIM{C} = OneImpurityOneBathDIM{C} where {C<:StarChainVacuumConfiguration}
# const ThermalOneBathDIM{C} = OneImpurityOneBathDIM{C} where {C<:StarChainBathConfiguration}




"""
	thermal_state(x::SISBD{C}, sys_states::Vector{Int}, symmetry; kwargs...)
	
Computing a thermal state for the system plus bath.

Different algorithms are used depending on the type of C:
* If C <: Union{ThermofieldConfiguration, StarVacuumConfiguration}, a quench is used
from the separable of a system state plus the thermal state of the bath
* If C <: StarChainBathConfiguration, then a purified thermal state of the system plus
bath is returned, in this case symmetry can only be NoSymmetry()
"""
function thermal_state(x::ThermofieldOneBathDIM; symmetry::FermionicSymmetry=SpinCharge(), sys_states::Vector{Int}=default_sys_states(x),
	quench::QuenchScheme=KohnSantoroQuench(ramp=4., flat=4.), 
	stepper::AbstractStepper=TDVPStepper(alg=TDVP2(stepsize=0.01, ishermitian=true, trunc=truncdimcutoff(D=200, ϵ=1.0e-6)))) 
	return qunched_thermal_state(x, sys_states, symmetry, quench, stepper)
end


