

"""
	SingleImpurityModel{C<:BathConfiguration} <: OneImpurityOneBathDIM{C}

Single Impurity Anderson Model
"""
struct FreeSISBD{C<:BathConfiguration} <: OneImpurityOneBathDIM{C}
	env::C
	μ::Float64
	J::Float64
end

FreeSISBD(env::BathConfiguration; μ::Real=0., J::Real=1.) = FreeSISBD(env, convert(Float64, μ), convert(Float64, J) )
SISBD(x::FreeSISBD; U::Real=0.) = SISBD(x.env, U=U, μ=x.μ, J=x.J)
function FreeSISBD(x::SISBD)
	(x.U == 0.) || throw(ArgumentError("U=0 is required by free fermionic model"))
	return FreeSISBD(x.env, μ=x.μ, J=x.J)
end
free(x::SISBD) = FreeSISBD(x)
FreeFermionicHamiltonian(x::FreeSISBD) = FreeFermionicHamiltonian(FermionicHamiltonian(x))
function thermal_state(x::FreeSISBD{<:StarChainConfiguration})
	h = free_hamiltonian_matrix(x)
	return free_thermal_state_util(h, x.bath)
end
quench(x::FreeSISBD, qs::QuenchScheme) = free_quench(x, qs)
function sys_ham(x::FreeSISBD)
	ham = FermionicHamiltonian()
	pos = only(default_sys_sites(x))
	push!(ham, TwoBodyTerm(pos, pos, coeff=x.μ))
	return ham
end
separable_state(x::FreeSISBD; kwargs...) = free_separable_state(x; kwargs...)


function free_thermal_state_util(h::AbstractMatrix, bath::AbstractFermionicBath)
	@assert ishermitian(h)
	h = convert(Matrix{eltype(h)}, h)
	evals, U = eigen(Hermitian(h))
	evals = [thermaloccupation(bath, item) for item in evals]
	return U * Diagonal(evals) * U'
end

free_hamiltonian_matrix(x::AbstractDIM) = coefficient_matrix(FreeFermionicHamiltonian(FermionicHamiltonian(x)))

function gf_greater_τ(x::FreeSISBD{<:StarChainConfiguration})
	h = free_hamiltonian_matrix(x)
	β, μ = x.bath.β, x.bath.μ
	c = μ .* one(h)
	sys_site = only(default_sys_sites(x))
	c[sys_site, sys_site] = 0
	ham = h - c
	λs, U = eigen(Hermitian(ham))
	ns = [fermidirac(β, 0., λs[k]) for k in 1:length(λs)]
	return τ -> _free_gf_imag(U, λs, ns, 1, 1, τ)
end
function gf_greater_τ(x::FreeSISBD{<:StarChainConfiguration}, τs::Vector{<:Real})
	f = gf_greater_τ(x)
	return f.(τs)
end 
function gf_greater_lesser_t(x::FreeSISBD, ρ₀::AbstractMatrix) 
	ham = free_hamiltonian_matrix(x)
	sys_site = only(default_sys_sites(x))
	return free_gf_real(ham, ρ₀, sys_site, sys_site)
end

function gf_greater_lesser_t(x::FreeSISBD, ρ₀::AbstractMatrix, ts::Vector{<:Real})
	f1, f2 = gf_greater_lesser_t(x, ρ₀)
	return f1.(ts), f2.(ts)
end


"""
	f_greater_lesser(x::FreeSISBD{<:StarConfiguration})
return free greater green function -i⟨cᵢ(t)cⱼ†⟩ and lesser green function -i⟨cⱼ†cᵢ(t)⟩, over the equilibrium state.
"""
function gf_greater_lesser_t(x::FreeSISBD{<:StarChainConfiguration}) 
	ham = free_hamiltonian_matrix(x)
	sys_site = only(default_sys_sites(x))
	return free_gf_util(ham, x.bath, sys_site, sys_site)
end

function gf_greater_lesser_t(x::FreeSISBD{<:StarChainConfiguration}, ts::Vector{<:Real}) 
	f = gf_greater_lesser_t(x)
	tmp = f.(ts)
	return [item[1] for item in tmp], [item[2] for item in tmp]
end

function free_gf_util(h::AbstractMatrix, bath::AbstractFermionicBath, i::Int, j::Int)
	@assert ishermitian(h)
	ham = convert(Matrix{eltype(h)}, h)
	evals, evecs = eigen(Hermitian(ham))
	L = size(ham, 1)
	function f(t::Number)
		r_g = zero(ComplexF64)
		r_l = zero(ComplexF64)
		for k in 1:L
			ss = evecs[i, k] * conj(evecs[j, k])
			exp_t = exp(-im * evals[k] * t)
			n_k = thermaloccupation(bath, evals[k])
			r_g += ss * (1 - n_k) * exp_t
			r_l += ss * n_k * exp_t
		end
		return -im * r_g, im * r_l
	end
	return f	
end

