

"""
	separable_state(x::SISBD{C}; symmetry::FermionicSymmetry=SpinCharge(), sys_states::Vector{Int}=default_sys_states(x)) 
	where {C <: Union{ThermofieldConfiguration, StarVacuumConfiguration}}

Return a separable state for quench
"""
separable_state(x::AbstractDIM; symmetry::FermionicSymmetry=SpinCharge(), sys_states::Vector{Int}=default_sys_states(x)) = separable_state_util(
	scalartype(x), separable_initial_state(x, sys_states), symmetry)  


function qunched_thermal_state(x::AbstractDIM, sys_states::Vector{Int}=default_sys_states(x), symmetry::FermionicSymmetry=SpinCharge(), 
	quench::QuenchScheme=KohnSantoroQuench(ramp=4., flat=4.), 
	stepper::AbstractStepper=TDVPStepper(alg=TDVP2(stepsize=0.01, ishermitian=true, trunc=truncdimcutoff(D=200, ϵ=1.0e-6))))

	state = complex(separable_state(x, symmetry=symmetry, sys_states=sys_states))
	h = consolidate(quench_hamiltonian(x, quench), symmetry=symmetry)

	stepper = similar(stepper, tspan=(0, -im * quench_time(quench)))
	state, cache = timeevo!(state, h, stepper)

	return state	
end

function free_separable_state(x::AbstractDIM; sys_states::Vector{Int}=default_sys_states(x)) 
	m = separable_initial_state(x, sys_states)
	ρ = diagm(m)
	return convert(Matrix{ComplexF64}, ρ)
end

# generate initial states

# is_vaild_state(x::String) = (x == "0") || (x == "↑↓")

# utility functions
is_vaild_state(x::Int) = (x == 0) || (x == 1)
function check_states(states::Vector{Int})
	for (i, s) in enumerate(states)
		is_vaild_state(s) || throw(ArgumentError("state $s (on site $i) is not valid (should be 0 or 1)"))
	end
end

function set_product_env_state_util!(init_state::Vector, x::StarConfiguration, envsites)
	m_freqs, _couplings = couplings(x)
	for (wj, site) in zip(m_freqs, envsites)
		nj = thermaloccupation(x.bath, wj)
		init_state[site] = nj
	end
	return init_state	
end
# set_product_env_state_util!(init_state::Vector, x::ChainConfiguration, envsites) = error("chain configuration does not have a simple product thermal state")
function set_product_env_state_util!(init_state::Vector, x::ThermofieldConfiguration, env1sites, env2sites)
	for site in env1sites
		init_state[site] = 0
	end
	for site in env2sites
		init_state[site] = 1
	end	
	return init_state
end
function set_product_env_state_util!(init_state::Vector, x::ChainVacuumConfiguration, envsites)
	# @warn "chain configuration does not have a simple product thermal state, this function should only be used as a DMRG initial guess"
	r = chain_vacuum_env_state_guess(x)
	for (site, state) in zip(envsites, r)
		init_state[site] = state
	end
	return init_state	
end

function chain_vacuum_env_state_guess(x::ChainVacuumConfiguration) 
	m_freqs, _couplings = couplings(x.star)
	alphas, betas, U = star_to_chain_full(m_freqs, _couplings, x.Nchain)
	f_occupations = [thermaloccupation(x.bath, item) for item in m_freqs]
	mat = U' * Diagonal(f_occupations) * U
	ncounts = round(Int, tr(mat))
	env_sites = collect(default_env_sites(x))
	@assert ncounts <= length(default_env_sites(x))
	r = zeros(Int, length(default_env_sites(x)))
	r[1:ncounts] .= 1
	return r
end
function set_product_sys_state_util!(init_state::Vector, x::AbstractDIM, sys_states)
	for (site, state) in zip(default_sys_sites(x), sys_states)
		init_state[site] = state
	end
	return init_state
end


function separable_state_util(::Type{T}, states::Vector{Int}, symmetry::ChargeCharge) where T
	check_states(states)
	ph = physical_space(symmetry)
	init_state = [(item==0) ? (0, 0) : (1, 1) for item in states]
	n1 = sum([item[1] for item in init_state])
	n2 = sum([item[2] for item in init_state])
	sector = ChargeChargeSector(up=n1, down=n2)
	return prodmps(T, [ph for i in 1:length(init_state)], init_state, right=space(sector) )
end
function separable_state_util(::Type{T}, states::Vector{Int}, symmetry::SpinCharge) where T
	check_states(states)
	ph = physical_space(symmetry)
	init_state = [(item==0) ? (0, 0) : (2, 0) for item in states]
	n = sum([item[1] for item in init_state])
	sector = SpinChargeSector(spin=0, charge=n)
	return prodmps(T, [ph for i in 1:length(init_state)], init_state, right=space(sector) )
end
function separable_state_util(::Type{T}, states::Vector{Int}, symmetry::Charge) where T
	check_states(states)
	ph = physical_space(symmetry)
	n = sum(states)
	sector = ChargeSector(charge=n)
	return prodmps(T, [ph for i in 1:length(states)], states, right=space(sector) )
end
separable_state_util(::Type{T}, states::Vector{<:Real}, symmetry::FermionicSymmetry) where T = separable_state_util(T, convert(Vector{Int}, states), symmetry)
