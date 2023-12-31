
abstract type DiscretizationScheme end

# interface, DiscretizationScheme should support .freqs, which is sorted from small to large
frequencies(x::DiscretizationScheme) = x.freqs

# num_sites(x::DiscretizationScheme) = length(frequencies(x)) - 1

"""
	spectrum_couplings(freqs::Vector{<:Real}, f; kwargs...)

Arguments
* f: pectrum function
* freqs: L frequences {w1, w2, ... , wL}, sorted from small to large

Return discretized frequencies Xn and hybridizations Vn of size at most L-1,
where |Vn|^2 = ∫_{w_{n}}^{w_{n+1}} dx f^2(x), and 
Xn = ∫_{w_{n}}^{w_{n+1}} dx x f^2(x) / |Vn|^2.
If |Vn| = 0 for some n, the corresponding Vn and Xn are removed.


Reference "How to discretize a quantum bath for real-time evolution".
Eqs.(10a, 10b) 
"""
function spectrum_couplings(freqs::Vector{<:Real}, f; atol::Real=1.0e-6, kwargs...)
	@assert atol > 0
	L = length(freqs)
	omegas = Float64[]
	couplings = Float64[]
	for i in 1:L-1
		v1, err1 = quadgk(f, freqs[i], freqs[i+1]; kwargs...)
		v2, err2 = quadgk(x -> x*f(x), freqs[i], freqs[i+1]; kwargs...)
		if v1 <= zero(v1)
			@warn "couplings is nonpositive for frequency interval $((round(freqs[i], digits=6), round(freqs[i+1], digits=6)))"
			v1 = zero(v1)
		end
		if v1 > atol
			omega = v2 / v1
			@assert (freqs[i] <= omega <= freqs[i+1]) 
			push!(omegas, omega)
			push!(couplings, sqrt(v1))
		end
	end
	return omegas, couplings
end
spectrum_couplings(x::DiscretizationScheme, f; kwargs...) = spectrum_couplings(frequencies(x), f; kwargs...)


struct LinearDiscretization <: DiscretizationScheme
	freqs::Vector{Float64}
	dw::Float64
end

function LinearDiscretization(;dw::Real, wmax::Real, wmin::Real)
	(wmax > wmin) || error("wmax must be larger than wmin")
	(dw > 0.) || error("dw must be positive.")
	return LinearDiscretization(collect(wmin:dw:wmax), dw)
end 

