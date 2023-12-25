function Base.getproperty(m::TwoBathDIM, s::Symbol)
	if s == :leftbath
		return m.leftenv.bath
	elseif s == :rightbath
		return m.rightenv.bath
	elseif s == :leftdiscretization
		return m.leftenv.discretization
	elseif s == :rightdiscretization
		return m.rightenv.discretization
	else
		return getfield(m, s)
	end
end


default_leftenv_sites(m::TwoBathDIM) = Iterators.reverse(default_env_sites(m.leftenv))
default_leftenv_a1_sites(m::TwoBathDIM{<:ThermofieldConfiguration}) = Iterators.reverse(default_env_a1_sites(m.leftenv))
default_leftenv_a2_sites(m::TwoBathDIM{<:ThermofieldConfiguration}) = Iterators.reverse(default_env_a2_sites(m.leftenv))


default_sys_sites(m::TwoBathDIM) = length(m.leftenv)+1:length(m.leftenv)+sys_size(m)
function default_boundary_sys_sites(m::TwoBathDIM)
	sites = default_sys_sites(m)
	return first(sites), last(sites)
end
default_env_sites(m::TwoBathDIM) = Iterators.flatten((default_leftenv_sites(m), default_rightenv_sites(m)))


_l_size(m::TwoBathDIM) = length(m.leftenv) + sys_size(m)

default_rightenv_sites(m::TwoBathDIM) = default_env_sites(m.rightenv) .+ _l_size(m)
default_rightenv_a1_sites(m::TwoBathDIM) = default_env_a1_sites(m.rightenv) .+ _l_size(m)
default_rightenv_a2_sites(m::TwoBathDIM) = default_env_a2_sites(m.rightenv) .+ _l_size(m)

default_env_a1_sites(m::TwoBathDIM) = Iterators.flatten((default_leftenv_a1_sites(m), default_rightenv_a1_sites(m)))
default_env_a2_sites(m::TwoBathDIM) = Iterators.flatten((default_leftenv_a2_sites(m), default_rightenv_a2_sites(m)))

Base.length(m::TwoBathDIM) = length(m.leftenv) + length(m.rightenv) + sys_size(m)
