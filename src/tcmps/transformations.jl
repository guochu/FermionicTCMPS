abstract type BathConfiguration{B<:AbstractFermionicBath} end

default_env_sites(x::BathConfiguration) = 1:length(x)

# 4 configurations for bosons and fermions

# star and star thermofield

"""
	StarConfiguration{B<:AbstractFermionicBath, D<:DiscretizationScheme} <: BathConfiguration{B}

Fields:
* Ws: the discretized on-site frequencies
* Vs: the discretized hybridization
"""
struct StarConfiguration{B<:AbstractFermionicBath, D<:DiscretizationScheme} <: BathConfiguration{B}
	bath::B
	discretization::D
	Ws::Vector{Float64}
	Vs::Vector{Float64}
end

function StarConfiguration(bath::AbstractFermionicBath, D::DiscretizationScheme; kwargs...)
	omegas, couplings = spectrum_couplings(frequencies(D), bath.spectrum; kwargs...)
	return StarConfiguration(bath, D, omegas, couplings)
end
function StarConfiguration(bath::AbstractFermionicBath; dw::Real, kwargs...)
	f = bath.spectrum
	lb, ub = lowerbound(f), upperbound(f)
	(isinf(lb) || isinf(ub)) && throw(ArgumentError("integrand can not contain Inf for LinearDiscretization"))
	return StarConfiguration(bath, LinearDiscretization(dw=dw, wmin=lb, wmax=ub); kwargs...)
end 
star(bath::AbstractFermionicBath, d::DiscretizationScheme; kwargs...) = StarConfiguration(bath, d; kwargs...)
star(bath::AbstractFermionicBath; kwargs...) = StarConfiguration(bath; kwargs...)
nfrequencies(x::StarConfiguration) = length(x.Ws)
Base.length(x::StarConfiguration) = nfrequencies(x)

struct StarThermofieldConfiguration{B<:AbstractFermionicBath, D<:DiscretizationScheme} <: BathConfiguration{B}
	star::StarConfiguration{B, D}
end

thermofield(m::StarConfiguration) = StarThermofieldConfiguration(m)
Base.length(x::StarThermofieldConfiguration) = 2*nfrequencies(x.star)


# chain and chain thermofield
struct ChainConfiguration{B<:AbstractFermionicBath, D<:DiscretizationScheme} <: BathConfiguration{B}
	star::StarConfiguration{B, D}
	Nchain::Int
end

chainmapping(m::StarConfiguration; Nchain::Int=nfrequencies(m)-1) = ChainConfiguration(m, Nchain)
Base.length(x::ChainConfiguration) = x.Nchain

struct ChainThermofieldConfiguration{B<:AbstractFermionicBath, D<:DiscretizationScheme} <: BathConfiguration{B}
	star::StarConfiguration{B, D}
	Nchain::Int
end

chain_thermofield(star::StarConfiguration; Nchain::Int=nfrequencies(star)-1) = ChainThermofieldConfiguration(star, Nchain)
thermofield(m::ChainConfiguration) = chain_thermofield(m.star, Nchain=m.Nchain)
chainmapping(m::StarThermofieldConfiguration; Nchain::Int=nfrequencies(m)-1) = chain_thermofield(m.star, Nchain=Nchain)
Base.length(x::ChainThermofieldConfiguration) = 2*x.Nchain


const ThermofieldConfiguration{B, D} = Union{StarThermofieldConfiguration{B, D}, ChainThermofieldConfiguration{B, D}} where {B<:AbstractFermionicBath, D<:DiscretizationScheme}
const StarChainConfiguration{B, D} = Union{StarConfiguration{B, D}, ChainConfiguration{B, D}} where {B<:AbstractFermionicBath, D<:DiscretizationScheme}
const StarVacuumConfiguration{B, D} = StarConfiguration{B, D} where {B<:FermionicVacuum, D}
const StarBathConfiguration{B, D} = StarConfiguration{B, D} where {B<:FermionicBath, D}
const ChainVacuumConfiguration{B, D} = ChainConfiguration{B, D} where {B<:FermionicVacuum, D}
const StarChainVacuumConfiguration{B, D} = StarChainConfiguration{B, D} where {B<:FermionicVacuum, D}
const StarChainBathConfiguration{B, D} = StarChainConfiguration{B, D} where {B<:FermionicBath, D}
const TransformedConfiguration{B, D} = Union{ChainConfiguration{B, D}, StarThermofieldConfiguration{B, D}, ChainThermofieldConfiguration{B, D}} where {B<:AbstractFermionicBath, D<:DiscretizationScheme}

default_env_a1_sites(x::ThermofieldConfiguration) = 1:2:length(x)
default_env_a2_sites(x::ThermofieldConfiguration) = 2:2:length(x)


nfrequencies(x::TransformedConfiguration) = nfrequencies(x.star)
function Base.getproperty(x::TransformedConfiguration, s::Symbol)
	if s == :bath
		return x.star.bath
	elseif s == :discretization
		return x.star.discretization
	else
		getfield(x, s)
	end
end

"""	
	couplings(m::StarConfiguration)
	couplings(m::ChainConfiguration)
	couplings(m::StarThermofieldConfiguration)
	couplings(m::ChainThermofieldConfiguration)

Return all the coupling coefficients for the impurity-bath couplings

For fermionic thermofield configuration, we have used particle-hole transformation transformation 
on top of the thermofield transformantion, which perserves the symmetry, in this case, even for vacuum, 
one can do a nontrivial thermofield transformantion plus particle-hole transformation. 

Reference arXiv:2012.01424v1.
"""
couplings(m::StarConfiguration) = m.Ws, m.Vs
function couplings(m::ChainConfiguration)
	m_freqs, _couplings = couplings(m.star)
	return star_to_chain(m_freqs, _couplings, m.Nchain)
end
function couplings(m::StarThermofieldConfiguration{<:AbstractFermionicBath}) 
	m_freqs, couplings_1, couplings_2 = _thermofield_couplings(m)
	return m_freqs, couplings_1, m_freqs, couplings_2
end 
function couplings(m::ChainThermofieldConfiguration{<:AbstractFermionicBath})
	m_freqs, couplings_1, couplings_2  = _thermofield_couplings(m)
	alpha1s, beta1s = star_to_chain_safe(m_freqs, couplings_1, m.Nchain)
	alpha2s, beta2s = star_to_chain_safe(m_freqs, couplings_2, m.Nchain)
	return alpha1s, beta1s, alpha2s, beta2s
end


function _thermofield_couplings(m::ThermofieldConfiguration)
	m_freqs, _couplings = couplings(m.star)

	couplings_1 = Float64[]
	couplings_2 = Float64[]
	for (wj, fj) in zip(m_freqs, _couplings)
		n_k = thermaloccupation(m.bath, wj)
		push!(couplings_1, fj * sqrt((1 - n_k)))
		push!(couplings_2, fj * sqrt(n_k))
	end

	return m_freqs, couplings_1, couplings_2	
end