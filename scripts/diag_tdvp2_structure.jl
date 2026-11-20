using Random, LinearAlgebra
function randh(n); X = randn(n,n) + im*randn(n,n); return (X+X')/2; end
Random.seed!(1)
I2 = Matrix(1.0*I, 2, 2); herm(X) = (X+X')/2
A1 = herm(randh(2)); B1 = randn(2,2)+im*randn(2,2); B2 = randn(2,2)+im*randn(2,2)
C2 = randn(2,2)+im*randn(2,2); D3 = herm(randh(2)); E3 = herm(randh(2))
H12 = herm(kron(kron(A1, I2), I2) + kron(kron(B1, B2), I2))
H23 = herm(kron(I2, kron(C2, D3)) + kron(I2, kron(I2, E3)))
H = H12 + H23
ψ1 = randn(ComplexF64,2); ψ2 = randn(ComplexF64,2); ψ3 = randn(ComplexF64,2)
ψ1 /= norm(ψ1); ψ2 /= norm(ψ2); ψ3 /= norm(ψ3)
# TRUE projected single-site H at site 2 (product-state neighbors):
w12a = real(dot(ψ1, A1, ψ1)); w12b = dot(ψ1, B1, ψ1)
w23c = dot(ψ3, C2, ψ3); w23e = real(dot(ψ3, E3, ψ3))
h2 = herm(w12a*I2 + w12b*B2 + w23c*C2 + w23e*I2)
h2_full = kron(I2, kron(h2, I2))
B23full = H23
ψ0 = kron(kron(ψ1, ψ2), ψ3)
upstream_step(dt) = begin
    h = dt/2
    ψ = exp(-im*h*H12) * ψ0                  # L-pass bond 1
    ψ = exp(-im*h*h2_full) * ψ               # L-pass S2 (−h)
    ψ = exp(-im*h*B23full) * ψ               # L-pass bond 2 (last, no reset)
    ψ = exp(-im*h*B23full) * ψ               # R-pass bond 2 (first)
    ψ = exp(-im*h*h2_full) * ψ               # R-pass S2 (−h)
    ψ = exp(-im*h*H12) * ψ                   # R-pass bond 1
    return ψ
end
println("--- upstream TDVP2 structure, TRUE projected h2 (NN) ---")
for dt in (0.02, 0.01, 0.005, 0.0025)
    println("dt=", dt, ":  one-step err = ", norm(upstream_step(dt) - exp(-im*dt*H)*ψ0))
end
# reference: plain Suzuki S2 = e^{-ihH12} e^{-idt H23} e^{-ihH12}
suzuki_step(dt) = exp(-im*(dt/2)*H12) * exp(-im*dt*H23) * exp(-im*(dt/2)*H12) * ψ0
println("--- plain Suzuki reference ---")
for dt in (0.02, 0.01, 0.005, 0.0025)
    println("dt=", dt, ":  one-step err = ", norm(suzuki_step(dt) - exp(-im*dt*H)*ψ0))
end
