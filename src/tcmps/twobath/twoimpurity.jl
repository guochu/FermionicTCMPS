
struct SIDBD{L<:BathConfiguration, R<:BathConfiguration} <: TwoBathDIM{L, R}
	leftenv::L
	rightenv::R
	U::Float64
	μ::Float64	
end
SIDBD(leftenv::BathConfiguration, rightenv::BathConfiguration; U::Real, μ::Real) = SIDBD(leftenv, rightenv, convert(Float64, U), convert(Float64, μ))
function sys_ham(x::SIDBD)
	h = FermionicHamiltonian()
	push!(h, TwoBodyTerm(1, 1, coeff=x.μ))
	push!(h, FourBodyTerm(1,1,1,1, coeff=x.U))
	return changepositions(h, default_sys_index_mapping(x))
end
sys_size(x::SIDBD) = 1

"""
	GIDBD{C1<:BathConfiguration, C2<:BathConfiguration}

Discrete Multiple Impurity Two Bath Model
The left bath couples to the left most impurity 
the right bath couples to the right most impurity
The final D means discrete
"""
struct GIDBD{C1<:BathConfiguration, C2<:BathConfiguration} <: TwoBathDIM{C1, C2}
	leftenv::C1
	rightenv::C2
	sys::FermionicHamiltonian
	nsys::Int
end

GIDBD(sys::FermionicHamiltonian, leftenv::BathConfiguration, rightenv::BathConfiguration) = GIDBD(leftenv, rightenv, sys, length(sys))
GIDBD(nsys::Int, leftenv::BathConfiguration, rightenv::BathConfiguration) = GIDBD(leftenv, rightenv, FermionicHamiltonian(), nsys)
sys_size(x::GIDBD) = x.nsys
# utility functions
sys_ham(x::GIDBD) = changepositions(x.sys, default_sys_index_mapping(x))
"""
	add_sys!(m::TwoFBathFSysImpurityModel, t::FermionicTerm) 
	add system operators
"""
function add_sys!(m::GIDBD, t::AbstractFermionicTerm)
	check_sys_positions(t, sys_size(m))
	push!(m.sys, t)
end 

function check_sys_positions_util(pos::Tuple, nsys::Int)
	for item in pos
		(1 <= item <= nsys) || error("system operator out of range (should be in [1,$(nsys)])")
	end
end
check_sys_positions(x::AbstractFermionicTerm, nsys::Int) = check_sys_positions_util(positions(x), nsys)
function check_sys_positions(x::AbstractFermionicHamiltonian, nsys::Int) 
	for item in x.data
		check_sys_positions_util(positions(item), nsys)
	end
end

left_sysbath_tunneling(x::TwoBathDIM{<:StarChainConfiguration}) = current_op(x.leftenv, default_leftenv_sites(x), first(default_sys_sites(x)))
left_sysbath_tunneling(x::TwoBathDIM{<:ThermofieldConfiguration}) = current_op(x.leftenv, default_leftenv_a1_sites(x), default_leftenv_a2_sites(x), first(default_sys_sites(x)))
right_sysbath_tunneling(x::TwoBathDIM{C, <:StarChainConfiguration}) where {C} = current_op(x.rightenv, default_rightenv_sites(x), last(default_sys_sites(x)))
right_sysbath_tunneling(x::TwoBathDIM{C, <:ThermofieldConfiguration}) where {C} = current_op(x.rightenv, default_rightenv_a1_sites(x), default_rightenv_a2_sites(x), last(default_sys_sites(x)))


"""
	struct IRLMD{L<:BathConfiguration, R<:BathConfiguration}

Two impurity two bath spinless fermionic model. 
The second last "S" means spinless/
"""
struct IRLMD{L<:BathConfiguration, R<:BathConfiguration} <: TwoBathDIM{L, R}
	leftenv::L
	rightenv::R
	μ::Float64	
	J::Float64
	U::Float64
end
IRLMD(leftenv::BathConfiguration, rightenv::BathConfiguration; μ::Real, J::Real, U::Real) = IRLMD(
	leftenv, rightenv, convert(Float64, μ), convert(Float64, J), convert(Float64, U))

function sys_ham(x::IRLMD)
	h = FermionicHamiltonian()
	push!(h, TwoBodyTerm(2, 2, coeff=x.μ-x.U))
	t = TwoBodyTerm(2, 1, coeff=x.J)
	push!(h, t)
	push!(h, t')
	t = TwoBodyTerm(2, 3, coeff=x.J)
	push!(h, t)
	push!(h, t')
	push!(h, FourBodyTerm(1,2,2,1, coeff=x.U))
	push!(h, FourBodyTerm(3,2,2,3, coeff=x.U))
	return changepositions(h, default_sys_index_mapping(x))
end
sys_size(x::IRLMD) = 3

function _add_leftenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{<:StarConfiguration})
	m_freqs_left, _couplings_left = couplings(x.leftenv)
	pos1, pos2 = default_boundary_sys_sites(x)
	for (j, wj, fj) in zip(default_leftenv_sites(x), m_freqs_left, _couplings_left)
		push!(ham, TwoBodyTerm(j, j, coeff=wj))
		# sys-bath interaction
		t = TwoBodyTerm(pos1, j, coeff=fj)
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
function _add_rightenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{C1, <:StarConfiguration}) where {C1}
	m_freqs_right, _couplings_right = couplings(x.rightenv)
	pos1, pos2 = default_boundary_sys_sites(x)
	for (j, wj, fj) in zip(default_rightenv_sites(x), m_freqs_right, _couplings_right)
		push!(ham, TwoBodyTerm(j, j, coeff=wj))
		# sys-bath interaction
		t = TwoBodyTerm(pos2, j, coeff= fj)
		push!(ham, t)
		push!(ham, t')
	end	
	return ham
end
function _add_leftenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{<:ChainConfiguration})
	alphas, betas = couplings(x.leftenv)

	sites = collect(default_leftenv_sites(x))
	for (j, alpha) in zip(sites, alphas)
		push!(ham, TwoBodyTerm(j,j,coeff=alpha))
	end
	for i in 1:length(sites)-1
		t = TwoBodyTerm(sites[i], sites[i+1], coeff=betas[i+1])
		push!(ham, t)
		push!(ham, t')
	end
	# sys-bath coupling
	j = sites[1]
	pos1, pos2 = default_boundary_sys_sites(x)
	t = TwoBodyTerm(pos1, j, coeff=betas[1])
	push!(ham, t)
	push!(ham, t')

	return ham
end
function _add_rightenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{C1, <:ChainConfiguration}) where {C1}
	alphas, betas = couplings(x.rightenv)
	sites = collect(default_rightenv_sites(x))
	for (j, alpha) in zip(sites, alphas)
		push!(ham, TwoBodyTerm(j,j,coeff=alpha))
	end
	for i in 1:length(sites)-1
		t = TwoBodyTerm(sites[i], sites[i+1], coeff=betas[i+1])
		push!(ham, t)
		push!(ham, t')
	end
	# sys-bath coupling
	j = sites[1]
	pos1, pos2 = default_boundary_sys_sites(x)
	t = TwoBodyTerm(pos2, j, coeff=betas[1])
	push!(ham, t)
	push!(ham, t')
	return ham	
end

function _add_leftenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{<:StarThermofieldConfiguration})
	m_freqs_1, couplings_1, m_freqs_2, couplings_2 = couplings(x.leftenv)
	pos1, pos2 = default_boundary_sys_sites(x)
	for (j1, wj1, fj1, j2, wj2, fj2) in zip(default_leftenv_a1_sites(x), m_freqs_1, couplings_1, default_leftenv_a2_sites(x), m_freqs_2, couplings_2)
		push!(ham, TwoBodyTerm(j1,j1,coeff=wj1))
		push!(ham, TwoBodyTerm(j2,j2,coeff=wj2))

		t = TwoBodyTerm(pos1, j1, coeff=fj1)
		push!(ham, t)
		push!(ham, t')
		t = TwoBodyTerm(pos1, j2, coeff=fj2)
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
function _add_rightenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{C1, <:StarThermofieldConfiguration}) where {C1}
	m_freqs_1, couplings_1, m_freqs_2, couplings_2 = couplings(x.rightenv)
	pos1, pos2 = default_boundary_sys_sites(x)
	for (j1, wj1, fj1, j2, wj2, fj2) in zip(default_rightenv_a1_sites(x), m_freqs_1, couplings_1, default_rightenv_a2_sites(x), m_freqs_2, couplings_2)
		push!(ham, TwoBodyTerm(j1,j1,coeff=wj1))
		push!(ham, TwoBodyTerm(j2,j2,coeff=wj2))

		t = TwoBodyTerm(pos2, j1, coeff=fj1)
		push!(ham, t)
		push!(ham, t')
		t = TwoBodyTerm(pos2, j2, coeff=fj2)		
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
function _add_leftenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{<:ChainThermofieldConfiguration})
	alpha1s, beta1s, alpha2s, beta2s = couplings(x.leftenv)
	pos1, pos2 = default_boundary_sys_sites(x)
	for (j1, alpha1, j2, alpha2) in zip(default_leftenv_a1_sites(x), alpha1s, default_leftenv_a2_sites(x), alpha2s)
		push!(ham, TwoBodyTerm(j1,j1,coeff=alpha1))
		push!(ham, TwoBodyTerm(j2,j2, coeff=alpha2))
	end
	a1_sites = collect(default_leftenv_a1_sites(x))
	a2_sites = collect(default_leftenv_a2_sites(x))
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

	t = TwoBodyTerm(pos1, j1, coeff=beta1s[1])
	push!(ham, t)
	push!(ham, t')
	t = TwoBodyTerm(pos1, j2, coeff=beta2s[1])
	push!(ham, t)
	push!(ham, t')
	return ham
end
function _add_rightenv!(ham::AbstractFermionicHamiltonian, x::TwoBathDIM{C1, <:ChainThermofieldConfiguration}) where {C1}
	alpha1s, beta1s, alpha2s, beta2s = couplings(x.rightenv)
	pos1, pos2 = default_boundary_sys_sites(x)
	for (j1, alpha1, j2, alpha2) in zip(default_rightenv_a1_sites(x), alpha1s, default_rightenv_a2_sites(x), alpha2s)
		push!(ham, TwoBodyTerm(j1,j1,coeff=alpha1))
		push!(ham, TwoBodyTerm(j2,j2, coeff=alpha2))
	end
	a1_sites = collect(default_rightenv_a1_sites(x))
	a2_sites = collect(default_rightenv_a2_sites(x))
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

	t = TwoBodyTerm(pos2, j1, coeff=beta1s[1])
	push!(ham, t)
	push!(ham, t')
	t = TwoBodyTerm(pos2, j2, coeff=beta2s[1])
	push!(ham, t)
	push!(ham, t')

	return ham
end

function FermionicHamiltonian(x::TwoBathDIM) 
	ham = sys_ham(x)
	_add_leftenv!(ham, x)
	_add_rightenv!(ham, x)
	return ham
end


function separable_initial_state(x::TwoBathDIM, sys_states::Vector{Int}=default_sys_states(x))
	check_states(sys_states)
	rho = zeros(length(x))
	set_product_sys_state_util!(rho, x, sys_states)
	_set_leftrho!(rho, x)
	_set_rightrho!(rho, x)
	return rho
end

function _set_leftrho!(rho::Vector, x::TwoBathDIM{<:StarChainConfiguration})
	return set_product_env_state_util!(rho, x.leftenv, default_leftenv_sites(x)) 
end
function _set_rightrho!(rho::Vector, x::TwoBathDIM{C1, <:StarChainConfiguration}) where {C1}
	return set_product_env_state_util!(rho, x.rightenv, default_rightenv_sites(x))
end

function _set_leftrho!(rho::Vector, x::TwoBathDIM{<:ThermofieldConfiguration}) 
	return set_product_env_state_util!(rho, x.leftenv, default_leftenv_a1_sites(x), default_leftenv_a2_sites(x) )
end
function _set_rightrho!(rho::Vector, x::TwoBathDIM{C1, <:ThermofieldConfiguration}) where {C1}
	return set_product_env_state_util!(rho, x.rightenv, default_rightenv_a1_sites(x), default_rightenv_a2_sites(x))
end


