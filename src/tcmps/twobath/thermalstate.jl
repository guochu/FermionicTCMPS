const ThermofieldTwoBathDIM{C1, C2} = TwoBathDIM{C1, C2} where {C1 <: ThermofieldConfiguration, C2 <: ThermofieldConfiguration}
const VacuumTwoBathDIM{C1, C2} = TwoBathDIM{C1, C2} where {C1 <: StarChainVacuumConfiguration, C2 <: StarChainVacuumConfiguration}
# const ThermalTwoBathDIM{C1, C2} = TwoBathDIM{C1, C2} where {C1<:StarChainBathConfiguration, C2<:StarChainBathConfiguration}

# what is the equilibrium state for two bath, if they do not have the same temperatue
# or the same chemical potential?
function thermal_state(x::TwoBathDIM; symmetry::FermionicSymmetry=SpinCharge(), 
	sys_states::Vector{Int}=default_sys_states(x),
	quench::QuenchScheme=KohnSantoroQuench(ramp=4., flat=4.), 
	stepper::AbstractStepper=TDVPStepper(alg=TDVP2(stepsize=0.01, ishermitian=true, trunc=truncdimcutoff(D=200, ϵ=1.0e-6))))
	check_same_beta_mu(x)
	return qunched_thermal_state(x, sys_states, symmetry, quench, stepper)
end

function check_same_beta_mu(x::TwoBathDIM)
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected"))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected"))	
end