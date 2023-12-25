# all the bath are in the thermofield configuration
const ThermofieldImpurityModel = Union{ThermofieldOneBathDIM, ThermofieldTwoBathDIM}

function _gf_util(model, i, j, symmetry)
	mapping = default_sys_index_mapping(model)
	ii = mapping[i]
	jj = mapping[j]
	adag = creation(ii, symmetry)
	a = creation(jj, symmetry)'
	ham = consolidate(FermionicHamiltonian(model), symmetry)
	return ham, prodmpo(physical_spaces(ham), adag), prodmpo(physical_spaces(ham), a)
end

function gf_greater_t(model::AbstractDIM, 
	times::Vector{<:Real}, i::Int=1, j::Int=i;
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),	
	th_state::AbstractMPS=thermal_state(model, symmetry=symmetry, sys_states = default_sys_states(model), quench=KohnSantoroQuench(ramp=4), stepper=gf_stepper))	

	ham, adag, a = _gf_util(model, j, i, symmetry)

	obs = correlation_2op_1t(ham, a, adag, th_state, times, stepper=gf_stepper, reverse=false)
	sca = -im 
	return obs .* sca
end

function gf_greater_τ(model::AbstractDIM, 
	times::Vector{<:Real}, i::Int=1, j::Int=i;
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	th_state::AbstractMPS=thermal_state(model, symmetry=symmetry, sys_states = default_sys_states(model), quench=KohnSantoroQuench(ramp=4), stepper=gf_stepper))
	
	ham, adag, a = _gf_util(model, j, i, symmetry)

	obs = correlation_2op_1τ(ham, a, adag, th_state, times, stepper=gf_stepper, reverse=false)
	sca = -1 
	return obs .* sca
end

function gf_lesser_t(model::AbstractDIM, 
	times::Vector{<:Real}, i::Int=1, j::Int=i;
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	th_state::AbstractMPS=thermal_state(model, symmetry=symmetry, sys_states = default_sys_states(model), quench=KohnSantoroQuench(ramp=4), stepper=gf_stepper))

	ham, adag, a = _gf_util(model, j, i, symmetry)

	obs = correlation_2op_1t(ham, adag, a, th_state, times, stepper=gf_stepper, reverse=true)
	sca = im 
	return obs .* sca
end

function gf_greater_lesser_t(model::AbstractDIM, 
	times::Vector{<:Real}, i::Int=1, j::Int=i;
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	th_state::AbstractMPS=thermal_state(model, symmetry=symmetry, sys_states = default_sys_states(model), quench=KohnSantoroQuench(ramp=4), stepper=gf_stepper))

	obs_greater = @spawn gf_greater_t(model, times, i, j, symmetry=symmetry, th_state=th_state, gf_stepper=gf_stepper)
	obs_lesser = gf_lesser_t(model, times, i, j, symmetry=symmetry, th_state=th_state, gf_stepper=gf_stepper)

	return fetch(obs_greater), obs_lesser
end


