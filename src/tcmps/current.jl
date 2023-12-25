
function current_op(config::StarConfiguration, envsites, sys_site)
	m_freqs, _couplings = couplings(config)
	h = FermionicHamiltonian()
	for (site, c) in zip(envsites, _couplings)
		push!(h, TwoBodyTerm(site, sys_site, coeff=c))
	end
	return h
end

function current_op(config::ChainConfiguration, envsites, sys_site)
	alphas, betas = couplings(config)
	h = FermionicHamiltonian()
	env_site = first(envsites)
	push!(h, TwoBodyTerm(env_site, sys_site, coeff=betas[1]))
	return h
end

function current_op(config::StarThermofieldConfiguration, env1sites, env2sites, sys_site)
	m_freqs_1, couplings_1, m_freqs_2, couplings_2 = couplings(config)
	h = FermionicHamiltonian()
	for (j1, wj1, fj1, j2, wj2, fj2) in zip(env1sites, m_freqs_1, couplings_1, env2sites, m_freqs_2, couplings_2)
		t = TwoBodyTerm(j1, sys_site, coeff=fj1)
		push!(h, t)
		t = TwoBodyTerm(j2, sys_site, coeff=fj2)
		push!(h, t)
	end
	return h
end

function current_op(config::ChainThermofieldConfiguration, env1sites, env2sites, sys_site)
	alpha1s, beta1s, alpha2s, beta2s = couplings(config)
	j1 = first(env1sites)
	j2 = first(env2sites)
	h = FermionicHamiltonian()
	t = TwoBodyTerm(j1, sys_site, coeff=beta1s[1])
	push!(h, t)
	t = TwoBodyTerm(j2, sys_site, coeff=beta2s[1])
	push!(h, t)
	return h
end	