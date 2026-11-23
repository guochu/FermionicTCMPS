# Dissect the W sparse tensors: are the opening/closing/propagator slots populated?
using Random, TensorKit, DMRG, FermionicTCMPS
using ImpurityModelBase: discretefermionicbath
using LinearAlgebra
const F = FermionicTCMPS
x = SISBD(discretefermionicbath([-0.9, 0.5], [0.35, 0.55]; β = 0.8, μ = 0.0); U = 0.0, εd = -0.3, J = 1.1)
h = mpo(x)
W = DMRG.timeevompo(h, -im * 0.0025, DMRG.WI())

function slotinfo(t::Union{DMRG.SchurMPOTensor, DMRG.SparseMPOTensor}, name)
    Os = t.Os
    s1, s2 = size(Os)
    println(name, ": Os size ", size(Os))
    for i in 1:s1, j in 1:s2
        e = Os[i, j]
        if e isa Number
            abs(e) > 1e-14 && println("   [$i,$j] scalar = ", round(e, digits = 6))
        else
            nrm = sqrt(real(TensorKit.inner(e, e)))
            nrm > 1e-14 && println("   [$i,$j] tensor norm = ", round(nrm, digits = 6))
        end
    end
end
slotinfo(h.data[1], "H site 1")
slotinfo(W.data[1], "W site 1")
slotinfo(W.data[2], "W site 2")
slotinfo(W.data[5], "W site 5")
