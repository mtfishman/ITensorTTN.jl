using ITensors
using ITensorTTN

N = 3
z = 3
χ = 5

@show N, z, χ

s = siteinds("S=1/2", N; conserve_qns=true)
ψ = InfiniteTTN(s, n -> isodd(n) ? "↑" : "↓"; z=z, linkdims=χ)

@assert degree(ψ) == z

for g in 0:(N + 2)
  println()
  for b in 1:(z - 1)
    @show g, b
    @show commonind(ψ[g, 1], ψ[g + 1, b])
  end
end

# Transfer matrix at generation 2 transforming
# from branch 1 to 2
t12 = TransferMatrix(ψ, 2, 1 => 2)
@show input_inds(t12)
@show output_inds(t12)

v1 = ITensor(input_inds(t12))
v2 = ITensor(output_inds(t12))
@show (dag(v2) * (t12 * v1))[]

# Transfer matrix at generation 2 transforming
# from branch 1 to 2 (towards the center)
t12 = TransferMatrix(ψ, 2, 1 => 2)

@show input_inds(t12)
@show output_inds(t12)

# Transfer matrix at generation 1 transforming
# from branch 2 to 1 (towards the center)
t21 = TransferMatrix(ψ, 1, 2 => 1)

@show input_inds(t21)
@show output_inds(t21)

# Composition of transfer matrices
t_tot = t21 * t12

@show input_inds(t_tot)
@show output_inds(t_tot)

v1 = ITensor(input_inds(t_tot))
@show inds(t_tot * v1)

