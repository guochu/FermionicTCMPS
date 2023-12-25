println("----------------------------------------------------")
println("|       Single Bath Non-Equilibrium Dynamics       |")
println("----------------------------------------------------")

@testset "NEQ dynamics: benchmark ED with thermofield configurations" begin

	β = 0.15
	dw = 0.2
	mu = 0.
	nchain = 4
	config = thermofield(chainmapping(star_config(β, mu, dw), Nchain=nchain))

	# evo_t = 0.5
	model = _free_model(config)
	hybrid_ham = coefficient_matrix(FreeFermionicHamiltonian(hybridization_hamiltonian(model)), length(model))
	energies_free = Float64[]

	qs = KohnSantoroQuench(ramp=4)
	ham = quench(model, qs)
	state = separable_state(model, sys_states = [0])


	delta = 0.1
	push!(energies_free, real(sum(hybrid_ham .* state)) )
	for i in 1:10
		stepper=RK4Stepper(tspan=(-im*(i-1)*delta , -im*i*delta), stepsize=0.05)
		state = timeevo!(state, ham, stepper=stepper)
		push!(energies_free, real(sum(hybrid_ham .* state)) )
	end

	tol = 1.0e-4

	model = SISBD(model)
	for symmetry in (ChargeCharge(), SpinCharge())
		hybrid_ham = consolidate(hybridization_hamiltonian(model), symmetry=symmetry)
		energies = Float64[]

		ham = consolidate(quench_hamiltonian(model, qs), symmetry=symmetry)
		state = complex(separable_state(model, symmetry=symmetry, sys_states=[0]))
		push!(energies, real(expectation(hybrid_ham, state) ) )

		local cache
		for i in 1:10
			# println("we are at step****************$i")
			stepper=TEBDStepper(tspan=(-im*(i-1)*delta , -im*i*delta), stepsize=0.05, order=4, trunc=truncdimcutoff(D=100, ϵ=1.0e-6))
			if !(@isdefined cache)
				state, cache = timeevo!(state, ham, stepper)
			else
				state, cache = timeevo!(state, ham, stepper, cache)
			end
			push!(energies, real(expectation(hybrid_ham, state) ) )
		end		

		@test maximum(abs.(energies_free - energies / 2)) < tol
	end
end

@testset "NEQ GF: benchmark ED with thermofield configurations" begin
	β = 1.
	δt = 0.02
	N = 20
	t = N * δt
	ts = collect(0:δt:t)

	μ = 0.
	dw = 2.
	ϵ_d = 1.25*pi

	trunc = truncdimcutoff(D=100, ϵ=1.0e-6)
	gf_stepper = TEBDStepper(tspan=(0 , -im*δt), stepsize=δt, order=2, trunc=trunc)
	
	bath = fermionicbath(spectrum_func(), β=β, μ=μ)
	config = star(bath, dw=dw)
	model = FreeSISBD(config, μ=ϵ_d)


	for sys_state in (0, 1)
		ρ₀ = separable_state(model, sys_states=[sys_state])
		a1, a2 = gf_greater_lesser_t(model, ρ₀, ts)
		a3 = a1 - a2
		if sys_state == 0
			@test norm(a2) < 1.0e-8
		else
			@test norm(a1) < 1.0e-8
		end
		
		for symmetry in (SpinCharge(), ChargeCharge())
			config1 = thermofield( chainmapping(config) )
			model1 = SISBD(config1, μ=ϵ_d, U=0)
			state = complex(separable_state(model1, symmetry=symmetry, sys_states=[sys_state]))
			b1, b2 = gf_greater_lesser_t(model1, ts, 1, 1, symmetry=symmetry, gf_stepper=gf_stepper, th_state=state)
			b3 = b1 - b2
			if sys_state == 0
				@test norm(b2) < 1.0e-8
			else
				@test norm(b1) < 1.0e-8
			end
			@test norm(a3 - b3) / norm(a3) < 1.0e-3
		end		
	end

end





