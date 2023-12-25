

function _first_order!(x, f, t, dt)
	# x .+= f(t + dt/2, x) .* dt
	axpy!(dt, f(t + dt/2, x), x)
end

function _rungekutta2order!(x, f, t, dt)
	k1 = f(t, x)
	k2 = f(t+dt, x + dt*k1)
	h = 0.5*dt
	axpy!(h, k1, x)
	axpy!(h, k2, x)
end


function _rungekutta4order!(x, f, t, dt)
	k1 = rmul!(f(t, x), dt)
	k2 = rmul!(f(t+dt/2, x + k1 / 2), dt)
	k3 = rmul!(f(t+dt/2, x + k2 / 2), dt)
	k4 = rmul!(f(t+dt, x + k3), dt)
	axpy!(1/6, k1, x)
	axpy!(1/3, k2, x)
	axpy!(1/3, k3, x)
	axpy!(1/6, k4, x)
end

struct RK2Stepper{T} <: AbstractStepper
	tspan::Tuple{T, T}
	stepsize::T
	n::Int
end

function RK2Stepper(;stepsize::Number, tspan::Tuple{<:Number, <:Number}=(0., stepsize))
	ti, tf = tspan
	delta = tf - ti
	n, stepsize = TEBD.compute_step_size(delta, stepsize)
	T = promote_type(typeof(ti), typeof(tf), typeof(stepsize))
	return RK2Stepper((convert(T, ti), convert(T, tf)), stepsize, n)
end

function Base.similar(x::RK2Stepper; tspan::Tuple{<:Number, <:Number}, stepsize::Number=x.stepsize)
	return RK2Stepper(tspan=tspan, stepsize=stepsize)
end

struct RK4Stepper{T} <: AbstractStepper
	tspan::Tuple{T, T}
	stepsize::T
	n::Int
end

function RK4Stepper(;stepsize::Number, tspan::Tuple{<:Number, <:Number}=(0., stepsize))
	ti, tf = tspan
	delta = tf - ti
	n, stepsize = TEBD.compute_step_size(delta, stepsize)
	T = promote_type(typeof(ti), typeof(tf), typeof(stepsize))
	return RK4Stepper((convert(T, ti), convert(T, tf)), stepsize, n)
end

function Base.similar(x::RK4Stepper; tspan::Tuple{<:Number, <:Number}, stepsize::Number=x.stepsize)
	return RK4Stepper(tspan=tspan, stepsize=stepsize)
end

function TEBD.timeevo!(state::AbstractMatrix, f::Union{FermionicCommutatorBase, TimeDependentFermionicCommutatorBase}, stepper::RK2Stepper)
	t_start = stepper.tspan[1]
	for i in 1:stepper.n
		_rungekutta2order!(state, f, t_start, stepper.stepsize)
		t_start += stepper.stepsize
	end
	return state
end

function TEBD.timeevo!(state::AbstractMatrix, f::Union{FermionicCommutatorBase, TimeDependentFermionicCommutatorBase}, stepper::RK4Stepper)
	t_start = stepper.tspan[1]
	for i in 1:stepper.n
		_rungekutta4order!(state, f, t_start, stepper.stepsize)
		t_start += stepper.stepsize
	end
	return state
end
# function timeevo(state::AbstractMatrix, f::FermionicCommutatorBase, stepper::ExactStepper)
# 	driver = stepper.ishermitian ? Lanczos() : Arnoldi()
# 	rho, info = exponentiate(f, stepper.δ, state, driver)
# 	(info.converged >= 1) || error("fail to converge")
# 	return rho
# end
function timeevo(state::AbstractMatrix, f::FermionicCommutatorBase, stepper::ExactStepper)
	dt = stepper.δ
	@assert real(dt) ≈ 0.
	λs = [exp(λ*dt) for λ in f.λs]
	exp_h = f.U * Diagonal(λs) * adjoint(f.U)
	return exp_h * state * exp_h'
end
function TEBD.timeevo!(state::AbstractMatrix, f::FermionicCommutatorBase, stepper::ExactStepper)
	rho = timeevo(state, f, stepper)
	return copyto!(state, rho)
end

TEBD.timeevo!(state::AbstractMatrix, f::Union{FermionicCommutatorBase, TimeDependentFermionicCommutatorBase}; stepper::AbstractStepper) = timeevo!(
	state, f, stepper)
