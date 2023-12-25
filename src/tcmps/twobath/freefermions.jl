struct FreeSIDBD{L<:BathConfiguration, R<:BathConfiguration} <: TwoBathDIM{L, R}
	leftenv::L
	rightenv::R
	μ::Float64	
end
FreeSIDBD(leftenv::BathConfiguration, rightenv::BathConfiguration; μ::Real) = FreeSIDBD(leftenv, rightenv, convert(Float64, μ))
function sys_ham(x::FreeSIDBD)
	h = FermionicHamiltonian()
	push!(h, TwoBodyTerm(1, 1, coeff=x.μ))
	return changepositions(h, default_sys_index_mapping(x))
end
sys_size(x::FreeSIDBD) = 1

function FreeSIDBD(x::SIDBD) 
	(x.U == 0.) || throw(ArgumentError("cannot convert an interacting model into a free model"))
	FreeSIDBD(x.leftenv, x.rightenv, x.μ)
end 
free(x::SIDBD) = FreeSIDBD(x)
SIDBD(x::FreeSIDBD; U::Real=0) = SIDBD(x.leftenv, x.rightenv, U=U, μ=x.μ)
FreeFermionicHamiltonian(x::FreeSIDBD) = FreeFermionicHamiltonian(FermionicHamiltonian(x))

quench(x::FreeSIDBD, qs::QuenchScheme) = free_quench(x, qs)
separable_state(x::FreeSIDBD; kwargs...) = free_separable_state(x; kwargs...)

# only in case of star or chain configuration we can initialize a sys+bath wise initial thermal state
function thermal_state(x::FreeSIDBD{<:StarChainConfiguration, <:StarChainConfiguration})
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected."))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected.")) 
	h = free_hamiltonian_matrix(x)
	return free_thermal_state_util(h, x.leftbath)
end

function gf_greater_lesser_t(x::FreeSIDBD, ρ₀::AbstractMatrix) 
	ham = free_hamiltonian_matrix(x)
	sys_site = only(default_sys_sites(x))
	return free_gf_real(ham, ρ₀, sys_site, sys_site)
end
function gf_greater_lesser_t(x::FreeSIDBD, ρ₀::AbstractMatrix, ts::Vector{<:Real})
	f1, f2 = gf_greater_lesser_t(x, ρ₀)
	return f1.(ts), f2.(ts)
end
function gf_greater_lesser_t(x::FreeSIDBD{<:StarChainConfiguration, <:StarChainConfiguration}) 
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected."))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected.")) 
	ham = free_hamiltonian_matrix(x)
	sys_site = only(default_sys_sites(x))
	return free_gf_util(ham, x.leftbath, sys_site, sys_site)
end
function gf_greater_lesser_t(x::FreeSIDBD{<:StarChainConfiguration, <:StarChainConfiguration}, ts::Vector{<:Real}) 
	f = gf_greater_lesser_t(x)
	tmp = f.(ts)
	return [item[1] for item in tmp], [item[2] for item in tmp]
end

struct FreeGIDBD{C1<:BathConfiguration, C2<:BathConfiguration} <: TwoBathDIM{C1, C2}
	leftenv::C1
	rightenv::C2
	sys::FreeFermionicHamiltonian
	nsys::Int
end
sys_size(x::FreeGIDBD) = x.nsys
sys_ham(x::FreeGIDBD) = changepositions(x.sys, default_sys_index_mapping(x))
function add_sys!(m::FreeGIDBD, t::TwoBodyTerm)
	check_sys_positions(t, sys_size(m))
	push!(m.sys, t)
end 


FreeGIDBD(sys::FreeFermionicHamiltonian, leftenv::BathConfiguration, rightenv::BathConfiguration) = FreeGIDBD(leftenv, rightenv, sys, length(sys))
FreeGIDBD(nsys::Int, leftenv::BathConfiguration, rightenv::BathConfiguration) = FreeGIDBD(leftenv, rightenv, FreeFermionicHamiltonian(), nsys)
GIDBD(x::FreeGIDBD) = GIDBD(FermionicHamiltonian(x.sys), x.leftenv, x.rightenv)
FreeGIDBD(x::GIDBD) = FreeGIDBD(x.leftenv, x.rightenv, FreeFermionicHamiltonian(x.sys), x.nsys)
free(x::GIDBD) = FreeGIDBD(x)
FreeFermionicHamiltonian(x::FreeGIDBD) = FreeFermionicHamiltonian(FermionicHamiltonian(x))


quench(x::FreeGIDBD, qs::QuenchScheme) = free_quench(x, qs)
separable_state(x::FreeGIDBD; kwargs...) = free_separable_state(x; kwargs...)


struct FreeIRLMD{L<:BathConfiguration, R<:BathConfiguration} <: TwoBathDIM{L, R}
	leftenv::L
	rightenv::R
	μ::Float64
	J::Float64
end
FreeIRLMD(leftenv::BathConfiguration, rightenv::BathConfiguration; μ::Real, J::Real) = FreeIRLMD(
	leftenv, rightenv, convert(Float64, μ), convert(Float64, J))
function FreeIRLMD(x::IRLMD)
	(x.U == 0.) || throw(ArgumentError("cannot convert an interacting model into a free model"))
	return FreeIRLMD(x.leftenv, x.rightenv, μ=x.μ, J=x.J)
end
sys_size(x::FreeIRLMD) = 3
function sys_ham(x::FreeIRLMD)
	h = FermionicHamiltonian()
	push!(h, TwoBodyTerm(2, 2, coeff=x.μ))
	t = TwoBodyTerm(2, 1, coeff=x.J)
	push!(h, t)
	push!(h, t')
	t = TwoBodyTerm(2, 3, coeff=x.J)
	push!(h, t)
	push!(h, t')
	return changepositions(h, default_sys_index_mapping(x))
end
IRLMD(x::FreeIRLMD; U::Real=0) = IRLMD(x.leftenv, x.rightenv, μ=x.μ, J=x.J, U=U)
free(x::IRLMD) = FreeIRLMD(x)
FreeFermionicHamiltonian(x::FreeIRLMD) = FreeFermionicHamiltonian(FermionicHamiltonian(x))

quench(x::FreeIRLMD, qs::QuenchScheme) = free_quench(x, qs)
separable_state(x::FreeIRLMD; kwargs...) = free_separable_state(x; kwargs...)

# only in case of star or chain configuration we can initialize a sys+bath wise initial thermal state
function thermal_state(x::FreeIRLMD{<:StarChainConfiguration, <:StarChainConfiguration})
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected."))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected.")) 
	h = free_hamiltonian_matrix(x)
	return free_thermal_state_util(h, x.leftbath)
end

function gf_greater_lesser_t(x::FreeIRLMD, ρ₀::AbstractMatrix, i::Int, j::Int) 
	@assert (1 <= i <= sys_size(x)) && (1 <= j <= sys_size(x))
	ham = free_hamiltonian_matrix(x)
	m = default_sys_index_mapping(x)
	return free_gf_real(ham, ρ₀, m[i], m[j])
end
function gf_greater_lesser_t(x::FreeIRLMD, ρ₀::AbstractMatrix,  i::Int, j::Int, ts::Vector{<:Real})
	f1, f2 = gf_greater_lesser_t(x, ρ₀, i, j)
	return f1.(ts), f2.(ts)
end
function gf_greater_lesser_t(x::FreeIRLMD{<:StarChainConfiguration, <:StarChainConfiguration}, i::Int, j::Int) 
	@assert (1 <= i <= sys_size(x)) && (1 <= j <= sys_size(x))
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected."))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected.")) 
	ham = free_hamiltonian_matrix(x)
	m = default_sys_index_mapping(x)
	return free_gf_util(ham, x.leftbath, m[i], m[j])
end
function gf_greater_lesser_t(x::FreeIRLMD{<:StarChainConfiguration, <:StarChainConfiguration}, i::Int, j::Int, ts::Vector{<:Real}) 
	f = gf_greater_lesser_t(x, i, j)
	tmp = f.(ts)
	return [item[1] for item in tmp], [item[2] for item in tmp]
end

# only in case of star or chain configuration we can initialize a sys+bath wise initial thermal state
function thermal_state(x::FreeGIDBD{<:StarChainConfiguration, <:StarChainConfiguration})
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected."))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected.")) 
	h = free_hamiltonian_matrix(x)
	return free_thermal_state_util(h, x.leftbath)
end

function gf_greater_lesser_t(x::FreeGIDBD{<:StarChainConfiguration, <:StarChainConfiguration}, ρ₀::AbstractMatrix, i::Int, j::Int=i) 
	ham = free_hamiltonian_matrix(x)
	m = default_sys_index_mapping(x)
	return free_gf_real(ham, ρ₀, m[i], m[j])
end
function gf_greater_lesser_t(x::FreeGIDBD{<:StarChainConfiguration, <:StarChainConfiguration}, ρ₀::AbstractMatrix, i::Int, j::Int, ts::Vector{<:Real})
	f1, f2 = gf_greater_lesser_t(x, ρ₀, i, j)
	return f1.(ts), f2.(ts)
end
function gf_greater_lesser_t(x::FreeGIDBD{<:StarChainConfiguration, <:StarChainConfiguration}, i::Int, j::Int=i) 
	(x.leftbath.β == x.rightbath.β) || throw(ArgumentError("same temperature for two baths is expected."))
	(x.leftbath.μ == x.rightbath.μ) || throw(ArgumentError("same chemical potential for two baths is expected.")) 
	check_sys_positions_util((i, j), sys_size(x))
	ham = free_hamiltonian_matrix(x)
	m = default_sys_index_mapping(x)
	return free_gf_util(ham, x.leftbath, m[i], m[j])
end
function gf_greater_lesser_t(x::FreeGIDBD{<:StarChainConfiguration, <:StarChainConfiguration}, i::Int, j::Int, ts::Vector{<:Real}) 
	f = gf_greater_lesser_t(x)
	tmp = f.(ts)
	return [item[1] for item in tmp], [item[2] for item in tmp]
end
