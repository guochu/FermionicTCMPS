

# utilities for exact time evolution

# do h^t ρ - ρ h^t
# time evolution for the coefficient matrix rho of free fermions is dρ/dt = -i [h^t, ρ]
# h is the coefficient matrix
# h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]
struct FermionicCommutatorBase
	h::Matrix{ComplexF64}
	U::Matrix{ComplexF64}
	λs::Vector{Float64}
end

function FermionicCommutatorBase(h::AbstractMatrix)
	h2 = convert(Matrix{ComplexF64}, h)
	λs, U = eigen(Hermitian(h2))
	return FermionicCommutatorBase(h2, U, λs)
end

function (op::FermionicCommutatorBase)(rho::AbstractMatrix)
	r = op.h * rho
	mul!(r, rho, op.h, -1., 1.)
	return r
end
(x::FermionicCommutatorBase)(t::Number, rho::AbstractMatrix) = x(rho)

struct TimeDependentFermionicCommutatorBase{F}
	h0::Matrix{ComplexF64}
	h1::Matrix{ComplexF64}
	f::F

function TimeDependentFermionicCommutatorBase{F}(h0::AbstractMatrix, h1::AbstractMatrix, f) where {F}
	((size(h0) == size(h1)) && (size(h0, 1)==size(h0, 2))) || throw(ArgumentError("square matrices expected."))
	return new{F}(convert(Matrix{ComplexF64}, h0), convert(Matrix{ComplexF64}, h1), f)
end

end


function (op::TimeDependentFermionicCommutatorBase)(t::Number, rho::AbstractMatrix)
	r = op.h0 * rho
	mul!(r, rho, op.h0, -1., 1.)
	v = op.f(t)
	mul!(r, op.h1, rho, v, 1.)
	mul!(r, rho, op.h1, -v, 1.)
	return r
end

# initializers
# do transpose here
FermionicCommutator(h::AbstractMatrix) = FermionicCommutatorBase(convert(Matrix{ComplexF64}, transpose(h)))
function FermionicCommutator(h0::AbstractMatrix, h1::AbstractMatrix, f)
	return TimeDependentFermionicCommutatorBase{typeof(f)}(transpose(h0), transpose(h1), f)
end

function free_gf_real(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int)
	ham = convert(Matrix{eltype(h)}, h)
	λs, U = eigen(Hermitian(ham))
	ns = [real(ρ₀[k, k]) for k in 1:size(U, 1)]
	ρ₁ = U' * Diagonal(1 .- ns) * U
	ρ₂ = U' * Diagonal(ns) * U
	return t -> -im*_free_gf_real(U, λs, ρ₁, i, j, t), t -> im*_free_gf_real(U, λs, ρ₂, i, j, t)
end

# HU = Uλs
# U brings H into diagonal
# ns are diagonal terms of ρ in the diagonal representation of H
function _free_gf_real(U, λs, ρ, i::Int, j::Int, t::Real)
	r_g = zero(eltype(U))
	L = size(U, 1)
	for k in 1:L
		for k′ in 1:L
			r_g += U[i, k] * conj(U[j, k′]) * ρ[k, k′] * exp(-im * λs[k] * t)
		end
	end
	return r_g
end

# function _free_gf_real(U, λs, ρ, i::Int, j::Int, t::Real)
# 	phase = [exp(-im * λs[k] * t) for k in 1:length(λs)]
# 	return sum(U[i, :] .* (Diagonal(phase) * ρ * conj(U[j, :])) )
# end



function _free_gf_imag(U, λs, ns, i::Int, j::Int, t::Real)
	r_g = zero(eltype(U))
	for k in 1:size(U, 1)
		ss = U[i, k] * conj(U[j, k])
		exp_t = exp(-λs[k] * t)
		n_k = ns[k]
		r_g += ss * (1 - n_k) * exp_t
	end
	return r_g
end




