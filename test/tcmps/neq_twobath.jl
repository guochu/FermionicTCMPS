println("----------------------------------------------------")
println("|         Two Bath Non-Equilibrium Dynamics        |")
println("----------------------------------------------------")

f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = SpectrumFunction(ϵ->f(D, ϵ), lb=-D, ub=D)


@testset "NEQ dynamics: benchmark ED with thermofield configurations" begin

	function _time_evo_quench(model1, model2, sys_states, symmetry)
		hybrid_ham = coefficient_matrix(FreeFermionicHamiltonian(hybridization_hamiltonian(model1)), length(model1))
		energies_free = Float64[]

		qs = KohnSantoroQuench(ramp=4)
		ham = quench(model1, qs)
		state = separable_state(model1, sys_states = sys_states)

		delta = 0.1
		push!(energies_free, real(sum(hybrid_ham .* state)) )
		for i in 1:5
			stepper=RK4Stepper(tspan=(-im*(i-1)*delta , -im*i*delta), stepsize=0.05)
			state = timeevo!(state, ham, stepper=stepper)
			push!(energies_free, real(sum(hybrid_ham .* state)) )
		end

		hybrid_ham = consolidate(hybridization_hamiltonian(model2), symmetry=symmetry)
		energies = Float64[]

		ham = consolidate(quench_hamiltonian(model2, qs), symmetry=symmetry)
		state = complex(separable_state(model2, symmetry=symmetry, sys_states=sys_states))
		push!(energies, real(expectation(hybrid_ham, state) ) )

		local cache
		for i in 1:5
			# println("we are at step****************$i")
			stepper=TEBDStepper(tspan=(-im*(i-1)*delta , -im*i*delta), stepsize=0.05, order=4, trunc=truncdimcutoff(D=50, ϵ=1.0e-6))
			if !(@isdefined cache)
				state, cache = timeevo!(state, ham, stepper)
			else
				state, cache = timeevo!(state, ham, stepper, cache)
			end
			push!(energies, real(expectation(hybrid_ham, state) ) )
		end
		
		return energies_free, energies / 2
	end

	function _time_evo(model1, model2, sys_states, symmetry)
		hybrid_ham = coefficient_matrix(FreeFermionicHamiltonian(hybridization_hamiltonian(model1)), length(model1))
		energies_free = Float64[]

		ham = FermionicCommutator( coefficient_matrix(FreeFermionicHamiltonian(model1)) )
		state = separable_state(model1, sys_states = sys_states)

		delta = 0.1
		push!(energies_free, real(sum(hybrid_ham .* state)) )
		for i in 1:5
			stepper=ExactStepper(tspan=(-im*(i-1)*delta , -im*i*delta), ishermitian=true)
			state = timeevo!(state, ham, stepper=stepper)
			push!(energies_free, real(sum(hybrid_ham .* state)) )
		end

		hybrid_ham = consolidate(hybridization_hamiltonian(model2), symmetry=symmetry)
		energies = Float64[]

		ham = consolidate(FermionicHamiltonian(model2), symmetry=symmetry)
		state = complex(separable_state(model2, symmetry=symmetry, sys_states=sys_states))
		push!(energies, real(expectation(hybrid_ham, state) ) )

		local cache
		for i in 1:5
			stepper = ExactStepper(tspan=(-im*(i-1)*delta , -im*i*delta), ishermitian=true)
			# println("we are at step****************$i")
			if !(@isdefined cache)
				state, cache = timeevo!(state, ham, stepper)
			else
				state, cache = timeevo!(state, ham, stepper, cache)
			end
			push!(energies, real(expectation(hybrid_ham, state) ) )
		end
		
		return energies_free, energies / 2
	end

	epsilon_d = -1.25  * pi
	leftmu = -2.
	rightmu = 5.
	leftbeta = 0.2
	rightbeta = 0.1
	Dl = 1.
	Dr = 10.
	Nl = 4
	leftdw = 0.25
	rightdw = 2.5
	Nr = 5
	tol = 1.0e-3

	leftbath = fermionicbath(spectrum_func(Dl), β=leftbeta, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=rightbeta, μ=rightmu)
	leftconfig = thermofield(chainmapping(star(leftbath, dw=leftdw), Nchain=Nl))
	rightconfig = thermofield(chainmapping(star(rightbath, dw=rightdw), Nchain=Nr))

	# two impurity two bath
	model1 = FreeGIDBD(2, leftconfig, rightconfig)
	add_sys!(model1, TwoBodyTerm(1,1,coeff=epsilon_d))
	add_sys!(model1, TwoBodyTerm(2,2,coeff=epsilon_d * 0.8))

	tunneling = TwoBodyTerm(1,2,coeff=0.5+0.6*im)
	add_sys!(model1, tunneling)
	add_sys!(model1, tunneling')

	model2 = GIDBD(model1)
	for sys_states in ([0,0], [0,1], [1,0], [1,1])
		# println("1 sys states is ", sys_states)
		for symmetry in (SpinCharge(), ChargeCharge())
			energies_free, energies = _time_evo_quench(model1, model2, sys_states, symmetry)
			@test maximum(abs.(energies_free - energies )) < tol
		end
	end

	# single impurity left bath + right vacuum
	leftbath = fermionicbath(spectrum_func(Dl), β=leftbeta, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=Inf, μ=rightmu)
	leftconfig = thermofield(chainmapping(star(leftbath, dw=leftdw), Nchain=2))
	rightdw = 5
	rightconfig = star(rightbath, dw=rightdw)
	model1 = FreeSIDBD(leftconfig, rightconfig, μ=epsilon_d)
	model2 = SIDBD(model1)

	
	for sys_states in ([0], [1])
		# println("2 sys states is ", sys_states)
		for symmetry in (SpinCharge(), ChargeCharge())
			energies_free, energies = _time_evo(model1, model2, sys_states, symmetry)
			@test maximum(abs.(energies_free - energies )) < 1.0e-8
		end
	end	

	# single impurity left vacum + right vacuum
	leftdw = 0.5
	leftbath = fermionicbath(spectrum_func(Dl), β=Inf, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=Inf, μ=rightmu)
	leftconfig = star(leftbath, dw=leftdw)
	rightconfig = star(rightbath, dw=rightdw)

	model1 = FreeSIDBD(leftconfig, rightconfig, μ=epsilon_d)
	model2 = SIDBD(model1)

	for sys_states in ([0], [1])
		# println("3 sys states is ", sys_states)
		for symmetry in (SpinCharge(), ChargeCharge())
			energies_free, energies = _time_evo(model1, model2, sys_states, symmetry)
			@test maximum(abs.(energies_free - energies )) < 1.0e-8
		end
	end	

end

@testset "NEQ GF: benchmark ED with thermofield configurations" begin

	function _time_evo(model1, model2, sys_states, symmetry)
		ts = collect(0:0.1:0.5)
		ρ₀ = separable_state(model1, sys_states=sys_states)
		a1, a2 = gf_greater_lesser_t(model1, ρ₀, ts)
		a3 = a1 - a2

		state = complex(separable_state(model2, symmetry=symmetry, sys_states=sys_states))
		# gf_stepper=TDVPStepper(tspan=(-im*(i-1)*0.1 , -im*i*0.1), alg = TDVP2(stepsize=0.05, ishermitian=true, trunc=truncdimcutoff(D=50, ϵ=1.0e-6)))
		gf_stepper=ExactStepper(tspan=(0 , -im*0.1), ishermitian=true)
		b1, b2 = gf_greater_lesser_t(model2, ts, symmetry=symmetry, gf_stepper=gf_stepper, th_state=state)
		b3 = b1 - b2
		return a1, a2, a3, b1, b2, b3
	end

	epsilon_d = -1.25  * pi
	leftmu = -2.
	rightmu = 5.
	leftbeta = 0.2
	rightbeta = 0.1
	Dl = 1.
	Dr = 10.
	leftdw = 1
	rightdw = 5
	tol = 1.0e-4

	# single impurity left bath + right vacuum
	leftbath = fermionicbath(spectrum_func(Dl), β=leftbeta, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=Inf, μ=rightmu)
	leftconfig = thermofield(chainmapping(star(leftbath, dw=leftdw)))
	rightconfig = star(rightbath, dw=rightdw)
	model1 = FreeSIDBD(leftconfig, rightconfig, μ=epsilon_d)
	model2 = SIDBD(model1)


	for sys_state in (0, 1)
		for symmetry in (SpinCharge(), ChargeCharge())
			a1, a2, a3, b1, b2, b3 = _time_evo(model1, model2, [sys_state], symmetry)
			if sys_state == 0
				@test norm(a2) < 1.0e-8
				@test norm(b2) < 1.0e-8
			else
				@test norm(a1) < 1.0e-8
				@test norm(b1) < 1.0e-8				
			end
			@test norm(a3 - b3) / norm(a3) < 1.0e-6
		end
	end

	# single impurity left vacum + right vacuum
	leftbath = fermionicbath(spectrum_func(Dl), β=Inf, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=Inf, μ=rightmu)
	leftconfig = star(leftbath, dw=leftdw)
	rightconfig = star(rightbath, dw=rightdw)

	model1 = FreeSIDBD(leftconfig, rightconfig, μ=epsilon_d)
	model2 = SIDBD(model1)

	for sys_state in (0, 1)
		for symmetry in (SpinCharge(), ChargeCharge())
			a1, a2, a3, b1, b2, b3 = _time_evo(model1, model2, [sys_state], symmetry)
			if sys_state == 0
				@test norm(a2) < 1.0e-8
				@test norm(b2) < 1.0e-8
			else
				@test norm(a1) < 1.0e-8
				@test norm(b1) < 1.0e-8				
			end
			@test norm(a3 - b3) / norm(a3) < 1.0e-6
		end
	end

	# spinless fermion, two bath
	function _time_evo(model1, model2, sys_states, i, j)
		ts = collect(0:0.1:0.5)
		ρ₀ = separable_state(model1, sys_states=sys_states)
		a1, a2 = gf_greater_lesser_t(model1, ρ₀, i, j, ts)
		a3 = a1 - a2

		state = complex(separable_state(model2, symmetry=Charge(), sys_states=sys_states))
		# gf_stepper=TDVPStepper(tspan=(-im*(i-1)*0.1 , -im*i*0.1), alg = TDVP2(stepsize=0.05, ishermitian=true, trunc=MPSTruncation(D=50, ϵ=1.0e-6)))
		gf_stepper=ExactStepper(tspan=(0 , -im*0.1), ishermitian=true)
		b1, b2 = gf_greater_lesser_t(model2, ts, i, j, symmetry=Charge(), gf_stepper=gf_stepper, th_state=state)
		b3 = b1 - b2
		return a3, b3
	end
	μ = -1
	model1 = FreeIRLMD(leftconfig, rightconfig, μ=μ, J=1)
	model2 = IRLMD(model1)

	for sys_states in ([0,0,0], [0,1,0], [1,0,1], [1,1,1])
		for i in 1:2, j in 1:2
			a, b = _time_evo(model1, model2, sys_states, i, j)
			@test norm(a - b) / norm(a) < 1.0e-6
		end
	end
end

