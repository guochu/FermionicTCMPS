

"""
	SISBD{C<:BathConfiguration} <: OneImpurityOneBathDIM{C}

Discrete Single Impurity Single Bath Anderson Model
"""
struct SISBD{C<:BathConfiguration} <: OneImpurityOneBathDIM{C}
	env::C
	U::Float64
	μ::Float64
	J::Float64
end

"""
	SingleImpurityModel(env::BathConfiguration; U::Real=0.)
"""
SISBD(env::BathConfiguration; U::Real=0., μ::Real=0., J::Real=1.) = SISBD(env, convert(Float64, U), convert(Float64, μ), convert(Float64, J) )

sysbath_tunneling(x::OneImpurityOneBathDIM{<:StarChainConfiguration}) = current_op(x.env, default_env_sites(x), only(default_sys_sites(x)))
sysbath_tunneling(x::OneImpurityOneBathDIM{<:ThermofieldConfiguration}) = current_op(x.env, default_env_a1_sites(x), default_env_a2_sites(x), only(default_sys_sites(x)))


function sys_ham(x::SISBD)
	ham = FermionicHamiltonian()
	pos = only(default_sys_sites(x))
	push!(ham, TwoBodyTerm(pos, pos, coeff=x.μ))
	if x.U != 0.
		push!(ham, FourBodyTerm(pos, pos, pos, pos, coeff=x.U))
	end
	return ham
end


function FermionicHamiltonian(x::OneImpurityOneBathDIM{<:StarConfiguration}) 
	m_freqs, _couplings = couplings(x.env)
	ham = sys_ham(x)
	sys_site = only(default_sys_sites(x))
	for (j, wj, fj) in zip(default_env_sites(x), m_freqs, _couplings)
		# on-site 
		push!(ham, TwoBodyTerm(j, j, coeff=wj))
		# sys-bath coupling
		t = TwoBodyTerm(sys_site, j, coeff=x.J*fj)
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end


function FermionicHamiltonian(x::OneImpurityOneBathDIM{<:ChainConfiguration})
	alphas, betas = couplings(x.env)
	ham = sys_ham(x)
	sites = collect(default_env_sites(x))
	for (j, alpha) in zip(sites, alphas)
		push!(ham, TwoBodyTerm(j, j, coeff=alpha))
	end
	for i in 1:length(sites)-1
		t = TwoBodyTerm(sites[i], sites[i+1], coeff=betas[i+1])
		push!(ham, t)
		push!(ham, t')
	end
	# sys-bath coupling
	j = sites[1]
	sys_site = only(default_sys_sites(x))
	t = TwoBodyTerm(sys_site, j, coeff=x.J*betas[1])
	push!(ham, t)
	push!(ham, t')
	return ham
end

"""
	The two baths c₁ and c₂ are interlacing each other
"""
function FermionicHamiltonian(x::OneImpurityOneBathDIM{<:StarThermofieldConfiguration})
	m_freqs_1, couplings_1, m_freqs_2, couplings_2 = couplings(x.env)
	ham = sys_ham(x)
	sys_site = only(default_sys_sites(x))
	for (j1, wj1, fj1, j2, wj2, fj2) in zip(default_env_a1_sites(x), m_freqs_1, couplings_1, default_env_a2_sites(x), m_freqs_2, couplings_2)
		push!(ham, TwoBodyTerm(j1, j1, coeff=wj1))
		push!(ham, TwoBodyTerm(j2, j2, coeff=wj2))
		t = TwoBodyTerm(sys_site, j1, coeff=x.J*fj1)
		push!(ham, t)
		push!(ham, t')
		t = TwoBodyTerm(sys_site, j2, coeff=x.J*fj2)
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end

function FermionicHamiltonian(x::OneImpurityOneBathDIM{<:ChainThermofieldConfiguration})
	alpha1s, beta1s, alpha2s, beta2s = couplings(x.env)
	ham = sys_ham(x)
	for (j1, alpha1, j2, alpha2) in zip(default_env_a1_sites(x), alpha1s, default_env_a2_sites(x), alpha2s)
		push!(ham, TwoBodyTerm(j1, j1, coeff=alpha1))
		push!(ham, TwoBodyTerm(j2, j2, coeff=alpha2))
	end
	a1_sites = collect(default_env_a1_sites(x))
	a2_sites = collect(default_env_a2_sites(x))
	for i in 1:length(a1_sites)-1
		t = TwoBodyTerm(a1_sites[i], a1_sites[i+1], coeff=beta1s[i+1])
		push!(ham, t)
		push!(ham, t')
		t = TwoBodyTerm(a2_sites[i], a2_sites[i+1], coeff=beta2s[i+1])
		push!(ham, t)
		push!(ham, t')
	end
	# sys-bath coupling
	j1 = a1_sites[1]
	j2 = a2_sites[1]
	sys_site = only(default_sys_sites(x))

	t = TwoBodyTerm(sys_site, j1, coeff=x.J*beta1s[1])
	push!(ham, t)
	push!(ham, t')
	t = TwoBodyTerm(sys_site, j2, coeff=x.J*beta2s[1])
	push!(ham, t)
	push!(ham, t')
	return ham
end



"""
	separable_initial_state(x::SISBD{<:ThermofieldConfiguration}, sys_states::Vector{Int}=default_sys_states(x))

Return a separable initial state in case of thermofield configuration.

The initial state of the bath is a proper thermal state, while the system 
state is a given product state. To obtain a global thermal state one still 
needs to do a quench.

For the quench see Reference "Efficient mapping for Anderson impurity problems with matrix product states".
"""
function separable_initial_state(x::OneImpurityOneBathDIM{<:ThermofieldConfiguration}, sys_states::Vector{Int}=default_sys_states(x))
	init_state = zeros(Int, length(x))
	set_product_sys_state_util!(init_state, x, sys_states)
	set_product_env_state_util!(init_state, x.env, default_env_a1_sites(x), default_env_a2_sites(x))
	return init_state
end 
function separable_initial_state(x::OneImpurityOneBathDIM{<:StarChainConfiguration}, sys_states::Vector{Int}=default_sys_states(x))
	init_state = zeros(Float64, length(x))
	set_product_sys_state_util!(init_state, x, sys_states)
	set_product_env_state_util!(init_state, x.env, default_env_sites(x))
	return init_state
end
