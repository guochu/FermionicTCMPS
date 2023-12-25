
const LANCZOS_ORTH_TOL = 1.0e-12


# ****************************************************************
# first step of lanczos iteration. v is the initial vector.
# the output a, b, are the first element of alpha, beta. qiminus1
# is the normalized initial vector, and q1 is the second lanczos
# vector
# the norm of initial input vector is outputed as vm
# ****************************************************************
function lanczos_make_first_step(A, v)
	vm = norm(v)
	qiminus1 = v / vm
	qi = A * qiminus1
	a = dot(qiminus1, qi)
	if abs(imag(a)) > LANCZOS_ORTH_TOL
		@warn "imaginary part the a is $(abs(imag(a))), larger than tolarence $LANCZOS_ORTH_TOL."
	end
	a = real(a)
	qi -= qiminus1 * a
	b = norm(qi)
	if b != 0.
	    qi /= b
	end
	return a, b, qiminus1, qi, vm
end

# ****************************************************************
# Assume we are in the i-th iteration. Then the input
# qiminus1 = V[i-1], qi=V[i], bold = beta[i-1],
# the output are alpha[i+1], beta[i], V[i+1], V[i]
# ***************************************************************
function lanczos_make_step(A, qiminus1, qi, bold)
	qiplus1 = A * qi
	# println(qi)
	# println(qiplus1)
	# println("----------------------------")
	a = dot(qi, qiplus1)

	if abs(imag(a)) > LANCZOS_ORTH_TOL
		@warn "imaginary part the a is $(abs(imag(a))), larger than tolarence $LANCZOS_ORTH_TOL."
	end
	a = real(a)
	qiplus1 -= (qi*a + qiminus1 * bold)
	b = norm(qiplus1)
	if b != 0.
	    qiplus1 /= b
	end
	qiminus1 = qi
	qi = qiplus1
	return a, b, qiminus1, qi
end

function fixed_step_lanczos_iteration_full(A, v0, kmax::Int=10)
	a, b, qiminus1, qi, vm = lanczos_make_first_step(A, v0)
	V = []
	push!(V, qiminus1)
	alphas = Vector{Float64}()
	betas = Vector{Float64}()
	push!(alphas, a)
	k = 1
	while k < kmax
		push!(betas, b)
	    a, b, qiminus1, qi = lanczos_make_step(A, qiminus1, qi, b)
	    push!(V, qiminus1)
	    push!(alphas, a)
	    k += 1
	end
	return alphas, betas, V
end

function star_to_chain_full(ws::Vector{Float64}, star_coupling::Vector{Float64}, Nchain::Int)
	(length(ws)==length(star_coupling)) || error("ws and start coupling size mismatch.")
	nrm = norm(star_coupling)
	(nrm == 0.) && error("the system-bath coupling is 0.")
	alphas, betas, vv = fixed_step_lanczos_iteration_full(Diagonal(ws), star_coupling/nrm, Nchain)
	pushfirst!(betas, nrm)
	return alphas, betas, hcat(vv...)
end
star_to_chain(ws::Vector{Float64}, star_coupling::Vector{Float64}, Nchain::Int) = star_to_chain_full(ws, star_coupling, Nchain)[1:2]

function star_to_chain_safe(ws::Vector{Float64}, star_coupling::Vector{Float64}, Nchain::Int)
	if norm(star_coupling) == 0.
		@warn "no system-bath coupling."
		return zeros(Nchain), zeros(Nchain)
	else
		return star_to_chain(ws, star_coupling, Nchain)
	end
end
