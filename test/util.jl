
f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = SpectrumFunction(ϵ->f(D, ϵ), lb=-D, ub=D)
spectrum_func() = spectrum_func(10)

function star_config(β, μ, dw)
	bath = fermionicbath(spectrum_func(), β=β, μ=μ)
	return star(bath, dw=dw)
end


function _free_model(config)
	ϵ_d = -1.25*pi
	model = FreeSISBD(config, μ=ϵ_d, J=1.1)
	return model
end
function free_model_thermal_occupy(config)
	model = _free_model(config)
	rho = thermal_state( model )
	sys_site = first(default_sys_sites(model))
	return real(rho[sys_site, sys_site])
end