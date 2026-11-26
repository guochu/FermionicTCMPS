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
ψ0 = fiducial_state(model)       # thermofield-double fiducial (sector-sum state)
evolve!(ψ0, h, (0.0, 1.0); δt = 0.02, trunc = truncdimcutoff(D = 256, ϵ = 1e-10))

occ = occupation_op(model; spin = :up)
⟨n_d↑⟩ = real(expectationvalue(occ, ψ0))
jop = current_op(model)          # J(t) = 2 Im⟨J⟩
gt, lt = gf_greater_lesser_t(model, 0:0.1:1.0; δt = 0.02)
```

Exports: `SISBD`, `SIDBD`, `mpo`, `fiducial_state`, `evolve!`,
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
* **Green's functions** (`src/observables.jl`): the thermofield-double fiducial is
  NOT a Hamiltonian eigenstate (the branch hoppings act on the doubles), so the
  correlators are sampled as sandwiches of forward-evolved states,
  `G^>(t) = −i ⟨ψ(t)| d_σ |χ(t)⟩` with `χ = d†|ψ₀⟩` and
  `G^<(t) = +i ⟨η(t)| d_σ |ψ(t)⟩` with `η = d|ψ₀⟩`, all chains evolved in
  lockstep. The insertions are applied mechanically per charge sector (exact
  per-configuration Jordan–Wigner strings, tree-level `_apply_dag`).

## Correctness gates

* **G2** (`test/g2_dense_anchor.jl`): the MPO is compared element-wise against a
  dense Hamiltonian assembled by mechanical Jordan–Wigner transformation
  (independent plain-matrix code) for one/two baths × vacuum/finite-β, at machine
  precision. This pins the intra-site JW phases, the string rule, the coefficient
  placement and the a₂ particle–hole convention.
* **G3** (`test/g3_dynamics.jl`): finite-β greater/lesser functions vs exact
  star-frame free-fermion evolution, occupied-spin lesser analytic value,
  impurity occupation vs exact star evolution, conserved-charge checks, and the
  sweep-drift regression gate of the plain pipeline.
* **Exact dynamics validation** (`scripts/exact_validation.jl`): three independent
  legs per config — dense mechanical-JW evolution, `DMRG.exact_timeevolution` on
  our MPO (ExactMPS, lossless), and the sweep-TDVP pipeline. The two exact legs
  agree to machine precision on every config; the sweep leg is quantified
  against them — see the accuracy table below.
* **G4** (`test/gtempo_crosscheck.jl` + `scripts/`): crosscheck against GTEMPO.
  The GTEMPO-side values are produced by
  `scripts/gtempo_benchmark.jl` (run inside the GTEMPO project) into a CSV that
  the test consumes; `scripts/compare_gtempo.jl` prints the full grid comparison.
  (Kept as two processes — the GTEMPO project pins its own dependency set, and
  mixing dev'd stacks in one environment is fragile.)
  The single-mode zero-temperature config agrees on the full complex greater
  function to ≤ 5e-3 (thermal impurity weight `(1−⟨n↑⟩)` per band, our spin-summed
  gt halved). The multi-mode finite-β configs show deviations growing with time
  from the sweep-TDVP cascade limitation below — compare with that envelope in
  mind (|Δ| ≈ 1e-3 up to t ≈ 0.3, O(0.1) by t = 1).

## Sweep accuracy (plain route)

The pipeline is plain-only: sharp thermofield-double fiducial, mechanical
insertions, two-site TDVP with the Kohn–Santoro minimum-bond crutch
(`add_back = MIND = 10`). The two-site projection seeds the far-pair amplitudes
only at O(δt·v), so the scheme is effectively FIRST order in δt (halving δt
halves the error — measured on the double-bath insertion overlap, and flat
against fresh vs incremental environments). The exact legs stay at machine
precision on every config class; the numbers below are the sweep-level
systematic only.

Maximum sweep error over t ≤ 1.0 at δt = 0.02, D = 256, ϵ = 1e-10
(`scripts/exact_validation.jl`, configs A–G, dense ED reference):

| Config class | max\|ΔG^>\| | max\|ΔG^<\| | max\|Δn_d\| |
|---|---|---|---|
| Vacuum above μ, cold impurity (C) | 0.0 (exact) | 0 | 0.0 |
| Vacuum above μ, hot impurity (F) | 1.2e-5 | 0 | 0.0 |
| Finite-β single bath, U = 0 (A) | 4.9e-3 | 0 | 7.1e-4 |
| Finite-β single bath, U = 1.7 (B) | 5.2e-3 | 0 | 7.0e-4 |
| Finite-β double bath, cold impurity (E) | 2.1e-3 | 0 | 5.0e-4 |
| Finite-β double bath, hot impurity (G) | 5.7e-4 | 1.1e-3 | 3.9e-4 |
| Zero-temperature Fermi sea straddling μ, U = 1.7 (D) | 1.8e-2 | 0 | 0.0 |

G^< is identically zero for a cold impurity (the annihilation kills the state)
and exact on every class. The vacuum-above-μ class is exact because its
fiducial reduces to the product vacuum; the straddling zero-temperature Fermi
sea is the worst class because the saturated \|11⟩ thermofield branches put the
largest thermal weight exactly on the far-pair directions the two-site window
cannot seed.

## Known limitations

The sweep-based two-site TDVP grows entanglement pair-by-pair: from a sharp
product state it cannot seed hopping terms whose support spans more than the
current pair when the intervening sites carry no direct coupling (e.g. the
escape into bath modes other than the first one at zero temperature, where
`v₂ = 0` removes the a₂-branch seeds). This is precisely the situation the
Kohn–Santoro equilibration ramp addresses; the ramp is intentionally out of
scope here (rev-2 decision). The mitigations that were implemented against it —
the equivalent initial-state expansion and the WI/WII sparse warm-up — were
REMOVED (2026-11-26 decision): the package is plain-route only, with the error
envelope of the accuracy table above. Consequences: zero-temperature multi-mode
baths and the thermal lesser function / `n_d(t)` / `J(t)` of the thermal
impurity are not captured by the bare separable evolution — the greater
function at finite temperature is fully validated (the a₂ branch seeds the
cascade). A practical size limit: the plain fiducial's bond carries
~4(K/2)² charge-bookkeeping labels (measured 400 at K = 10), which collides
with the default D = 256 truncation — beyond K ≈ 6–8 the cap would cut charge
sectors and rewrite the thermal weight distribution.

## Running tests

```sh
julia --project=test test/runtests.jl
```

Dependencies: `TensorKit` (0.17 stack via the ported engine), `DMRG` (dev path
`../DMRG`, the Missile-level engine snapshot — the TensorKit 0.17 port with
`MatrixAlgebraKit` and an `ExpExp` dependency, both dev'd from the Missile tree),
`ImpurityModelBase` (dev path `../ImpurityModelBase`), `LinearAlgebra`.

If the shared `~/.julia/environments/v1.10` pins an incompatible TensorKit,
isolate Pkg operations: `JULIA_LOAD_PATH="@:@stdlib" julia --project=. -e '...'`.
