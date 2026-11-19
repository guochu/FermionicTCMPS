# FermionicTCMPS

Thermo-field chain-mapping solver for fermionic Anderson impurity models on a
TensorKit U₁×U₁ symmetric tensor stack. The single-impurity Anderson model with a
discrete fermionic bath is doubled via the Kohn–Santoro thermo-field mapping
([arXiv:2012.01424](https://arxiv.org/abs/2012.01424)); the impurity dynamics is
then obtained by real-time evolution of a separable product state with two-site
TDVP, with the impurity state stored as an MPS.

## Scope

```julia
using FermionicTCMPS

ws = [-0.9, -0.3, 0.5]           # bath mode frequencies
fs = [0.35, 0.45, 0.5]           # spectral values (couplings = √fs)
bath = discretefermionicbath(ws, fs; β = 1.0, μ = 0.0)   # or discretefermionicvacuum
model = SISBD(bath; U = 1.7, εd = -0.3, J = 1.1)         # single impurity, single bath
model = SIDBD(bathl, bathr; U, εd, Jl, Jr)               # single impurity, two baths

h  = mpo(model)                  # FSM MPO (bond dimension 6, constant along the chain)
ψ0 = separable_state(model)      # |imp⟩ ⊗ |∅₁⟩ ⊗ |F₂⟩ thermofield product state
evolve!(ψ0, h, (0.0, 1.0); δt = 0.02, trunc = truncdimcutoff(D = 256, ϵ = 1e-10))

occ = occupation_op(model; spin = :up)
⟨n_d↑⟩ = real(expectationvalue(occ, ψ0))
jop = current_op(model)          # J(t) = 2 Im⟨J⟩
gt, lt = gf_greater_lesser_t(model, 0:0.1:1.0; δt = 0.02)
```

Exports: `SISBD`, `SIDBD`, `mpo`, `separable_state`, `evolve!`,
`occupation_op`, `current_op`, `gf_greater_lesser_t`.

## Constructions

* **FSM MPO** (`src/mpoconstruction.jl`): each hopping direction is one finite-state
  channel; coefficients sit at the non-shared end (the impurity is shared by all
  terms). Channel bond slots carry single U₁×U₁ charges, so the bond dimension is
  constant (vacuum + 4 channels + final = 6). Every term carries the explicit
  Jordan–Wigner string `c̃_a† (∏_{min ≤ t < max} P_t) c̃_b` with intra-site phases
  folded into the local tables (`src/siteops.jl`).
* **Thermo-field mapping** (`src/models.jl`): each bath mode splits into an a₁
  branch (coupling `V√(1−n)`) and an a₂ branch in the particle–hole representation
  (coupling `V√n`, ordinary hopping in hole variables, `+ω n_h` on-site), sites
  interleaved `(a₁, a₂)` per mode. Exact per Kohn–Santoro; the supplemented chain
  Hamiltonian is diagonal in the fiducial bath state.
* **Green's functions** (`src/observables.jl`): the fiducial state is an
  eigenstate of the chain Hamiltonian, so the real-time impurity correlators
  reduce to same-charge overlaps of freely evolved insertion chains,
  `G^>(t) = −i e^{iEt} ⟨χ|χ(t)⟩` with `χ = d†|ψ₀⟩` and
  `G^<(t) = i e^{−iEt} conj⟨η|η(t)⟩` with `η = d|ψ₀⟩`.

## Correctness gates

* **G2** (`test/g2_dense_anchor.jl`): the MPO is compared element-wise against a
  dense Hamiltonian assembled by mechanical Jordan–Wigner transformation
  (independent plain-matrix code) for one/two baths × vacuum/finite-β, at machine
  precision. This pins the intra-site JW phases, the string rule, the coefficient
  placement and the a₂ particle–hole convention.
* **G3** (`test/g3_dynamics.jl`): finite-β greater/lesser functions vs exact
  star-frame free-fermion evolution, occupied-spin lesser analytic value,
  impurity occupation vs exact star evolution, conserved-charge and δt²-convergence
  checks.
* **G4** (`test/gtempo_crosscheck.jl` + `scripts/`): crosscheck against GTEMPO.
  The two stacks pin different TensorKit majors (0.13 vs 0.17) and cannot share an
  environment, so the GTEMPO-side values are produced by
  `scripts/gtempo_benchmark.jl` (run inside the GTEMPO project) into a CSV that
  the test consumes; `scripts/compare_gtempo.jl` prints the full grid comparison.
  The single-mode zero-temperature config agrees on the full complex greater
  function to ≤ 5e-3 (thermal impurity weight `(1−⟨n↑⟩)` per band, our spin-summed
  gt halved). The multi-mode finite-β configs show deviations growing with time
  from the sweep-TDVP cascade limitation below — compare with that envelope in
  mind (|Δ| ≈ 1e-3 up to t ≈ 0.3, O(0.1) by t = 1).

## Known limitations

The sweep-based two-site TDVP grows entanglement pair-by-pair: from a sharp
product state it cannot seed hopping terms whose support spans more than the
current pair when the intervening sites carry no direct coupling (e.g. the
escape into bath modes other than the first one at zero temperature, where
`v₂ = 0` removes the a₂-branch seeds). This is precisely the situation the
Kohn–Santoro equilibration ramp addresses; the ramp is intentionally out of
scope here (rev-2 decision). Consequences: zero-temperature multi-mode baths
and the thermal lesser function / `n_d(t)` / `J(t)` of the thermal impurity are
not captured by the bare separable evolution — the greater function at finite
temperature is fully validated (the a₂ branch seeds the cascade).

## Running tests

```sh
julia --project=test test/runtests.jl
```

Dependencies: `TensorKit` (pinned 0.13), `DMRG` (dev path `../DMRG`, the Missile-level engine snapshot),
`ImpurityModelBase` (dev path), `LinearAlgebra`.
