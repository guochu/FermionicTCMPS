abstract type OneImpurityOneBathDIM{B} <: OneBathDIM{B} end
sys_size(m::OneImpurityOneBathDIM) = 1

function Base.getproperty(m::OneImpurityOneBathDIM, s::Symbol)
	if s == :bath
		return m.env.bath
	elseif s == :discretization
		return m.env.discretization
	else
		return getfield(m, s)
	end
end


default_sys_sites(m::OneBathDIM) = 1:sys_size(m)
default_env_sites(m::OneBathDIM) = default_env_sites(m.env) .+ sys_size(m)

"""
	env_a1_sites(m::AbstractSingleBathImpurityModel{<:StarThermofieldConfiguration})
	env_a1_sites(m::AbstractSingleBathImpurityModel{<:ChainThermofieldConfiguration})
	positions for the first auxiliary bath a1 after thermofield transformation, we chose an interlacing convention.
	From bosoinc bath, both a1 and a2 should be initialized as 0
	From fermionic bath, a1 is initialized as 0 and a2 is initialized as 1 due to the PHT for a2
"""
default_env_a1_sites(m::OneBathDIM{<:ThermofieldConfiguration}) = default_env_a1_sites(m.env) .+ sys_size(m)
default_env_a2_sites(m::OneBathDIM{<:ThermofieldConfiguration}) = default_env_a2_sites(m.env) .+ sys_size(m)
Base.length(m::OneBathDIM) = sys_size(m) + length(m.env)




