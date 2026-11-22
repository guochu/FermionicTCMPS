# Verify: expanded states are right-canonical on sites 2..L (site 1 = mixed
# carrier, the form the sweep expects), norm 1, exactly equivalent to the base
# chain, and print the bond dimensions.

using TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath, discretefermionicvacuum
using FermionicTCMPS: _sector, _amplitude, _adag_up, Charge
include(joinpath(@__DIR__, "..", "test", "anchor_util.jl"))

function check(name, x, phis_fn)
    spec = FermionicTCMPS.chain_spec(x)
    L = length(spec.kinds)
    imp = spec.imp_site
    worst = 0.0
    for spin in (:up, :down)
        delta = 1
        χ = FermionicTCMPS._shifted_state(x, 0, 0, spin, delta)
        χ === nothing && continue
        q_out = spin === :up ? (1, 0) : (0, 1)
        occs0 = ntuple(i -> i == imp ? q_out : (0, 0), L)
        exp = FermionicTCMPS.expanded_shifted_state(x, 0, 0, spin, delta, phis_fn(occs0))
        canon = all(DMRG.isrightcanonical, exp.data[2:end])
        d0 = abs(dot(χ, exp) - 1.0)
        dn = abs(dot(exp, exp) - 1.0)
        dφ = maximum(abs(dot(_basis_state(spec, Tuple(φ); T = ComplexF64), exp)) for φ in phis_fn(occs0))
        nrm = norm(exp)
        println(name, " spin=$spin: sites2:L right-canonical = $canon,  norm-1 = ",
                round(abs(nrm - 1); digits = 12), "  ⟨χ|exp⟩-1 = ", d0,
                "  ⟨exp|exp⟩-1 = ", dn, "  max⟨φ|exp⟩ = ", dφ)
        worst = max(worst, d0, dn, dφ, abs(nrm - 1))
        canon || (worst = Inf)
    end
    return worst
end

# bond dimensions: recompute from far_pair_directions count
function bond_dims(name, x)
    spec = FermionicTCMPS.chain_spec(x)
    imp = spec.imp_site
    L = length(spec.kinds)
    occs0 = ntuple(i -> i == imp ? (1, 0) : (0, 0), L)
    phis = FermionicTCMPS.far_pair_directions(spec)(occs0)
    println(name, ": directions = ", length(phis), "  interior bond dim = ", 1 + 2 * length(phis),
            "  (boundary bonds = 1)")
end

xc = SISBD(discretefermionicvacuum([0.4, 0.9], [0.55, 0.35]; μ = 0.0); U = 1.7, εd = 0.2, J = 1.0)
check("C", xc, FermionicTCMPS.far_pair_directions(FermionicTCMPS.chain_spec(xc)))
bond_dims("C", xc)

xl = discretefermionicbath([-0.8], [0.6]; β = 0.8, μ = 0.0)
xr = discretefermionicbath([0.5], [0.45]; β = 0.8, μ = 0.0)
xg = SIDBD(xl, xr; U = 1.7, εd = -0.3, Jl = 1.0, Jr = 1.0)
check("G", xg, FermionicTCMPS.far_pair_directions(FermionicTCMPS.chain_spec(xg)))
bond_dims("G", xg)
