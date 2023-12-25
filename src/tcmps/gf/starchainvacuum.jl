# the star/chain configuration with vacuum bath.

const VacuumImpurityModel = Union{VacuumOneBathDIM, VacuumTwoBathDIM}

# gf for StarChainVacuumConfiguration, not implemented for two baths
function gs_gf_greater_t(model::VacuumImpurityModel, times::Vector{<:Real}, i::Int=1, j::Int=i; 
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	gs::AbstractMPS=complex(ground_state(model, symmetry=symmetry, alg=DMRG2())[2])) 

	ham, adag, a = _gf_util(model, j, i, symmetry)

	obs = gs_correlation_2op_1t(ham, a, adag, gs, times, stepper=gf_stepper, reverse=false)
	sca = -im 
	return obs .* sca
end

function gs_gf_greater_τ(model::VacuumImpurityModel, times::Vector{<:Real}, i::Int=1, j::Int=i; 
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	gs::AbstractMPS=ground_state(model, symmetry=symmetry, alg=DMRG2())[2])

	ham, adag, a = _gf_util(model, j, i, symmetry)

	obs = gs_correlation_2op_1τ(ham, a, adag, gs, times, stepper=gf_stepper, reverse=false)
	sca = -1 
	return obs .* sca
end

function gs_gf_lesser_t(model::VacuumImpurityModel, times::Vector{<:Real}, i::Int=1, j::Int=i; 
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	gs::AbstractMPS=complex(ground_state(model, symmetry=symmetry, alg=DMRG2())[2]))

	ham, adag, a = _gf_util(model, j, i, symmetry)

	obs = gs_correlation_2op_1t(ham, adag, a, gs, times, stepper=gf_stepper, reverse=true)
	sca = im 
	return obs .* sca
end

function gs_gf_greater_lesser_t(model::VacuumImpurityModel, times::Vector{<:Real}, i::Int=1, j::Int=i; 
	symmetry::FermionicSymmetry=SpinCharge(), 
	gf_stepper::AbstractStepper=TEBDStepper(stepsize=0.05, order=4),
	gs::AbstractMPS=complex(ground_state(model, symmetry=symmetry, alg=DMRG2())[2])) 

	obs_greater = @spawn gs_gf_greater_t(model, times, i, j, symmetry=symmetry, gs=gs, gf_stepper=gf_stepper)
	obs_lesser = gs_gf_lesser_t(model, times, i, j, symmetry=symmetry, gs=gs, gf_stepper=gf_stepper)

	return fetch(obs_greater), obs_lesser
end


DMRG.ground_state(x::VacuumImpurityModel;symmetry::FermionicSymmetry=SpinCharge(), alg::MPSAlgorithm=DMRG2()) = _ground_state(x, symmetry, alg)

function _ground_state(x, symmetry::FermionicSymmetry, alg::MPSAlgorithm)
	h = consolidate(FermionicHamiltonian(x), symmetry=symmetry)
	mpo = MPO(h)
	# init_state = separable_initial_state(x, sys_states)
	# init_state = [Int(item) for item in init_state]

	dim = CartesianIndices(Tuple(2 for i in 1:length(default_sys_states(x))))

	local energy
	local state
	energies = Float64[]
	for index in dim
		sys_states = [index[i]-1 for i in 1:length(index)]
		_energy, _state =  _ground_state_single(x, symmetry, sys_states=sys_states, alg=alg)
		if !(@isdefined energy)
			energy = _energy
			state = _state
		else
			if _energy < energy
				energy = _energy
				state = _state
			end
		end
		push!(energies, _energy)
	end
	# println("all energies $energies")
	return energy, state	
end

function _ground_state_single(x, symmetry::FermionicSymmetry=SpinCharge(); sys_states::Vector{Int}=default_sys_states(x), alg::MPSAlgorithm=DMRG2()) 
	h = consolidate(FermionicHamiltonian(x), symmetry)
	mpo = MPO(h)
	init_state = separable_initial_state(x, sys_states)

	init_state_tmp = separable_state_util(scalartype(mpo), init_state, symmetry)

	return ground_state(mpo, alg, right=space_r(init_state_tmp)')
end	

