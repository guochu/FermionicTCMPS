abstract type AbstractDIM end
abstract type OneBathDIM{B<:BathConfiguration} <: AbstractDIM end
abstract type TwoBathDIM{L<:BathConfiguration, R<:BathConfiguration} <: AbstractDIM end


DMRG.scalartype(x::AbstractDIM) = Float64
# interfaces
FermionicHamiltonian(x::AbstractDIM) = error("FermionicHamiltonian not implemented for model type $(typeof(x))")
default_env_sites(x::AbstractDIM) = error("default_env_sites not implemented for model type $(typeof(x))")
sys_size(x::AbstractDIM) = error("sys_size not implemented for model type $(typeof(x))")
sys_ham(x::AbstractDIM) = error("sys_ham not implemented for model type $(typeof(x))")

default_sys_states(x::AbstractDIM) = [0 for site in default_sys_sites(x)]
default_sys_index_mapping(m::AbstractDIM) = Dict(i=>pos for (i, pos) in enumerate(default_sys_sites(m)))
separable_initial_state(x::AbstractDIM, args...) = error("separable_initial_state is not defined for model type $(typeof(x))")

# general functions
function quench_hamiltonian(x::AbstractDIM, qs::QuenchScheme)
	hsys, henv = _split_env_hamiltonian(x)
	coeff = Coefficient(t->qs(t))
	return coeff * hsys + henv
end

function free_quench(x::AbstractDIM, qs::QuenchScheme)
	hsys, henv = _split_env_hamiltonian(x)
	hsys = coefficient_matrix(FreeFermionicHamiltonian(hsys), length(x))
	henv = coefficient_matrix(FreeFermionicHamiltonian(henv), length(x))
	return FermionicCommutator(henv, hsys, t->qs(t))
end

function is_env_term(x::AbstractDIM, m::AbstractFermionicTerm)
	env_sites = default_env_sites(x)
	for item in positions(m)
		(item in env_sites) || return false
	end
	return true	
end

function _split_env_hamiltonian(x::AbstractDIM)
	h = FermionicHamiltonian(x)
	hsys = FermionicHamiltonian()
	henv = FermionicHamiltonian()
	for item in h.data
		if is_env_term(x, item)
			push!(henv, item)
		else
			push!(hsys, item)
		end
	end
	return hsys, henv
end

function hybridization_hamiltonian(x::AbstractDIM)
	h = FermionicHamiltonian(x)
	hybrid_h = FermionicHamiltonian()
	for item in h.data
		if is_hybrid_term(x, item)
			push!(hybrid_h, item)
		end
	end
	return hybrid_h
end

function is_hybrid_term(x::AbstractDIM, m::AbstractFermionicTerm)
	has_a = false
	has_b = false
	sys_sites = default_sys_sites(x)
	env_sites = default_env_sites(x)
	for item in positions(m)
		if item in sys_sites
			has_a = true
		end
		if item in env_sites
			has_b = true
		end
	end
	return has_a && has_b	
end

function environment_hamiltonian(x::AbstractDIM)
	h = FermionicHamiltonian(x)
	henv = FermionicHamiltonian()
	for item in h.data
		if is_env_term(x, item)
			push!(henv, item)
		end
	end	
	return henv
end
