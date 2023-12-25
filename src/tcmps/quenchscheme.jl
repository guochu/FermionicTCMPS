
abstract type QuenchScheme end

# interface to be implemented
quench_time(x::QuenchScheme) = error("quench_time not implemented for $(typeof(x))")
(op::QuenchScheme)(t::Real) = error("typeof(op)(t) not implemented.")
# work around for hamiltonian evolution with -im * t
(op::QuenchScheme)(t::Complex) = op(abs(t))

struct KohnSantoroQuench <: QuenchScheme
	ramp::Float64
	flat::Float64	
end

KohnSantoroQuench(;ramp::Real, flat::Real=ramp) = KohnSantoroQuench(convert(Float64, ramp), convert(Float64, flat))
KohnSantoroQuench(t::Real) = KohnSantoroQuench(ramp=t/2)

ramp_time(x::KohnSantoroQuench) = x.ramp
flat_time(x::KohnSantoroQuench) = x.flat
quench_time(x::KohnSantoroQuench) = ramp_time(x) + flat_time(x)

function (op::KohnSantoroQuench)(t::Real)
	qt = quench_time(op)
	if (t >= qt) && (t - qt) <= 1.0e-14
		t = qt
	end
	# ((t >= 0.) && (t <= qt) ) || error("time out of range.")
	if t >= ramp_time(op)
		return 1.
	else
		return t / ramp_time(op) 
	end
end



