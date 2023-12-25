println("----------------------------------------------")
println("|         Single Bath Equilibrium GF         |")
println("----------------------------------------------")


@testset "Ground state Occupation: benchmark ED with star and chain configurations" begin
	β = Inf
	μ = 0.
	config1 = star_config(Inf, 0., 4)
	config2 = chainmapping(config1)
	fms = [_free_model(config1), _free_model(config2)]

	D = 40
	tol = 1.0e-6

	for _m in fms
		rho = thermal_state(_m)
		n1s = [rho[i, i] for i in 1:size(rho, 1)]
		# TCMPS
		model = SISBD(_m)
		for symmetry in (ChargeCharge(), SpinCharge())
			observers = [consolidate(TwoBodyTerm(i, i), symmetry=symmetry) for i in 1:length(model)]
			for alg in (DMRG2(trunc=truncdimcutoff(D=D, ϵ=1.0e-8)), ED())
				energy, state = ground_state(model, symmetry=symmetry, alg = alg)
				n2s = [real(expectation(item, state)) / 2 for item in observers]
				@test maximum(abs.(n1s - n2s)) < tol	
			end
		end
	end
end

@testset "Ground state GF: benchmark ED with star and chain configurations" begin
	dw = 5.
	ϵ_d = -1.25*pi
	μ = 0.
	β = Inf

	# check ED, noninteracting
	bath = fermionicbath(spectrum_func(), β=β, μ=μ)
	config1 = star(bath, dw=dw)
	config2 = chainmapping(config1)
	times = collect(0:0.1:1)
	tol = 1.0e-6
	for config in (config1, config2)
		model = SISBD(config, μ=ϵ_d, J=1, U=0.)
		a2, b2 = gf_greater_lesser_t(free(model), times)
		for symmetry in (ChargeCharge(), SpinCharge())
			gs = complex(ground_state(model, symmetry=symmetry, alg=ED())[2])
			a1, b1 = gs_gf_greater_lesser_t(model, times, 1, 1, symmetry=symmetry, gs=gs, gf_stepper=ExactStepper(tspan=(0, 0.01), ishermitian=true))
			@test max(maximum(abs.(a1 - a2)), maximum(abs.(b1 - b2))) < tol
		end
	end

	# interacting
	ϵ_d = -1.
	bath = fermionicbath(spectrum_func(), β=β, μ=μ)
	config1 = star(bath, dw=dw)
	config2 = chainmapping(config1)
	times = collect(0:0.1:1.3)
	stepper=TDVPStepper(alg = TDVP2(stepsize=0.05, ishermitian=true, trunc=truncdimcutoff(D=50, ϵ=1.0e-5)))
	gs_alg = DMRG2(trunc=truncdimcutoff(D=100, ϵ=1.0e-8))
	for config in (config1, config2)
		model = SISBD(config, μ=ϵ_d, J=1, U=1.5)
		for symmetry in (SpinCharge(), ChargeCharge())
			gs = complex(ground_state(model, symmetry=symmetry, alg=gs_alg)[2])
			a1, b1 = gs_gf_greater_lesser_t(model, times, 1, 1, symmetry=symmetry, gs=gs, gf_stepper=stepper)
			gs = complex(ground_state(model, symmetry=symmetry, alg=ED())[2])
			a2, b2 = gs_gf_greater_lesser_t(model, times, 1, 1, symmetry=symmetry, gs=gs, gf_stepper=ExactStepper(tspan=(0, 0.01), ishermitian=true))
			@test max(maximum(abs.(a1 - a2)), maximum(abs.(b1 - b2))) < tol
		end
	end
end
