println("----------------------------------------------")
println("|          Two Bath Equilibrium GF           |")
println("----------------------------------------------")


@testset "Interacting resonant level model: benchmark ED with star and chain configurations" begin
	# this test fails for leftmu = rightmu !=0, ther may be some issue to check for eq gf.
	μ = -1.25  * pi
	leftmu = -0.
	rightmu = leftmu
	Dl = 1.
	Dr = 10.
	leftdw = 1
	rightdw = 5
	tol = 1.0e-4
	symmetry = Charge()

	leftbath = fermionicbath(spectrum_func(Dl), β=Inf, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=Inf, μ=rightmu)
	leftconfig1 = star(leftbath, dw=leftdw)
	rightconfig1 = star(rightbath, dw=rightdw)

	leftconfig2 = chainmapping(leftconfig1)
	rightconfig2 = chainmapping(rightconfig1)
	times = collect(0:0.1:0.5)

	model = IRLMD(leftconfig1, rightconfig1, μ=μ, J=1, U=0.)
	gs = complex(ground_state(model, symmetry=symmetry, alg=ED())[2])
	for i in 1:3, j in 1:3
		a2, b2 = gf_greater_lesser_t(free(model), i, j, times)
		a1, b1 = gs_gf_greater_lesser_t(model, times, i, j, symmetry=symmetry, gs=gs, gf_stepper=ExactStepper(tspan=(0, 0.1), ishermitian=true))
		@test max(maximum(abs.(a1 - a2)), maximum(abs.(b1 - b2))) < tol
	end


	stepper=TDVPStepper(alg = TDVP2(stepsize=0.05, ishermitian=true, trunc=truncdimcutoff(D=50, ϵ=1.0e-5)))
	gs_alg = DMRG2(trunc=truncdimcutoff(D=100, ϵ=1.0e-8))
	
	for c1 in (leftconfig1, leftconfig2), c2 in (rightconfig1, rightconfig2)
		model = IRLMD(c1, c2, μ=μ, J=1, U=1.3)
		gs_ed = complex(ground_state(model, symmetry=symmetry, alg=ED())[2])
		gs_dmrg = complex(ground_state(model, symmetry=symmetry, alg=gs_alg)[2])
		for i in 1:2, j in 1:2
			a2, b2 = gs_gf_greater_lesser_t(model, times, i, j, symmetry=symmetry, gs=gs_dmrg, gf_stepper=stepper)
			a1, b1 = gs_gf_greater_lesser_t(model, times, i, j, symmetry=symmetry, gs=gs_ed, gf_stepper=ExactStepper(tspan=(0, 0.01), ishermitian=true))
			@test max(maximum(abs.(a1 - a2)), maximum(abs.(b1 - b2))) < tol
		end
	end	
end