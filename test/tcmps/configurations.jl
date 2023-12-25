println("------------------------------------")
println("|      Bath Configurations         |")
println("------------------------------------")

@testset "Star configuration" begin
	# scaling with number of sites
	β = 0.25
	bath = fermionicbath(spectrum_func(), β=β)

	config1 = star(bath, dw=0.1)
	@test isa(config1, StarConfiguration)
	n1 = free_model_thermal_occupy(config1)
	config2 = star(bath, dw=0.02)	
	n2 = free_model_thermal_occupy(config2)
	@test abs((n1 - n2) / n1) < 1.0e-4
	config3 = star(bath, dw=0.01)
	n3 = free_model_thermal_occupy(config3)
	@test abs((n1 - n3) / n1) < 1.0e-4

	# versus chain configuration
	β = 0.2
	config1 = star_config(0.2, -0.5, 0.1)
	n1 = free_model_thermal_occupy(config1)
	config2 = chainmapping(config1)
	n2 = free_model_thermal_occupy(config2)
	@test abs((n1 - n2) / n1) < 1.0e-4
end


@testset "Particle Current in different configurations" begin

	function compute_current_1(model)
		δt = 0.02
		N = 10
		t = N * δt
		ham = FermionicCommutator(coefficient_matrix(FreeFermionicHamiltonian(model)))
		observer = coefficient_matrix(FreeFermionicHamiltonian(sysbath_tunneling(model)), length(model))

		ρ₀ = separable_state(model, sys_states=[0])
		currents = [sum(observer .* ρ₀)]
		for i in 1:N
			stepper = ExactStepper(tspan=(-im*(i-1)*δt, -im*i*δt), ishermitian=true)
			timeevo!(ρ₀, ham, stepper)
			push!(currents, sum(observer .* ρ₀))
		end	
		return currents	
	end

	β = 1.
	D = 10.
	μ = 5
	dw = 0.2
	ϵ_d = 1.25*pi
	tol = 1.0e-10

	# single bath
	bath = fermionicbath(spectrum_func(D), β=β, μ=μ)
	config = star(bath, dw=dw)
	model = FreeSISBD(config, μ=ϵ_d)
	currents_1 = compute_current_1(model)
	model = FreeSISBD(thermofield(config), μ=ϵ_d)
	currents_2 = compute_current_1(model)
	@test norm(currents_1 - currents_2) / norm(currents_1) < tol
	model = FreeSISBD(chainmapping(thermofield(config)), μ=ϵ_d)
	currents_3 = compute_current_1(model)
	@test norm(currents_1 - currents_3) / norm(currents_1) < tol

	function compute_current_2(model)
		δt = 0.02
		N = 10
		t = N * δt
		ham = FermionicCommutator(coefficient_matrix(FreeFermionicHamiltonian(model)))
		observer1 = coefficient_matrix(FreeFermionicHamiltonian(left_sysbath_tunneling(model)), length(model))
		observer2 = coefficient_matrix(FreeFermionicHamiltonian(right_sysbath_tunneling(model)), length(model))

		ρ₀ = separable_state(model, sys_states=[1])
		currents1 = [sum(observer1 .* ρ₀)]
		currents2 = [sum(observer2 .* ρ₀)]
		for i in 1:N
			stepper = ExactStepper(tspan=(-im*(i-1)*δt, -im*i*δt), ishermitian=true)
			timeevo!(ρ₀, ham, stepper)
			push!(currents1, sum(observer1 .* ρ₀))
			push!(currents2, sum(observer2 .* ρ₀))
		end	
		return currents1, currents2	
	end

	# two baths
	epsilon_d = -1.25  * pi
	U = 2.5 * pi
	leftmu = -2.
	rightmu = 5.
	leftbeta = 0.2
	rightbeta = 0.1
	Dl = 1.
	Dr = 10.
	leftdw = 0.02
	rightdw = 0.2

	leftbath = fermionicbath(spectrum_func(Dl), β=leftbeta, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=rightbeta, μ=rightmu)
	leftconfig = star(leftbath, dw=leftdw)
	rightconfig = star(rightbath, dw=rightdw)
	model = FreeSIDBD(leftconfig, rightconfig, μ=epsilon_d)
	a1, a2 = compute_current_2(model)
	leftconfigs = [leftconfig, thermofield(leftconfig), chainmapping(thermofield(leftconfig))]
	rightconfigs = [rightconfig, thermofield(rightconfig), chainmapping(thermofield(rightconfig))]

	for l in leftconfigs, r in rightconfigs
		m = FreeSIDBD(l, r, μ=epsilon_d)
		t1, t2 = compute_current_2(model)
		@test norm(t1-a1) / norm(a1) < tol
		@test norm(t2-a2) / norm(a2) < tol
	end

end
