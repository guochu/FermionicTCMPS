abstract type AbstractFermionicHamiltonian end

Hamiltonians.qterms(x::AbstractFermionicHamiltonian) = x.data
Base.isempty(x::AbstractFermionicHamiltonian) = isempty(x.data)
Base.length(x::AbstractFermionicHamiltonian) = maximum_site(x)

DMRG.scalartype(x::AbstractFermionicHamiltonian) = DMRG.compute_scalartype(x.data)
function Hamiltonians.isconstant(x::AbstractFermionicHamiltonian)
	for item in x.data
		isconstant(item) || return false
	end
	return true
end

function maximum_site(x::AbstractFermionicHamiltonian)
	@assert !isempty(x)
	m = 0
	for item in x.data
		m = max(m, maximum(positions(item)))
	end
	return m
end

function consolidate(x::AbstractFermionicHamiltonian, symmetry::FermionicSymmetry)
	terms = []
	for item in x.data
		_add!(terms, consolidate(item, symmetry))
	end
	isempty(terms) && error("vanishing hamiltonian")
	ph = physical_space(symmetry)
	physpaces = [ph for i in 1:maximum_site(x)]
	return QuantumOperator(physpaces, terms)
end
consolidate(x::AbstractFermionicHamiltonian; symmetry::FermionicSymmetry=SpinCharge()) = consolidate(x, symmetry)

_add!(terms::Vector, m::QTerm) = push!(terms, m)
_add!(terms::Vector, m::Vector) = append!(terms, m)

