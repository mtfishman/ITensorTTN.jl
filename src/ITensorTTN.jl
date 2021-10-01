module ITensorTTN

using ITensors
using ITensorInfiniteMPS
using Reexport

@reexport using ITensorInfiniteMPS

using ITensorInfiniteMPS: celltags, CelledVector, tag_starting_with

export InfiniteTTN, degree

struct InfiniteTTN
  data::CelledVector{ITensor}
end

function degree(ψ::InfiniteTTN)
  zs = order.(ψ.data[Cell(1)])
  z1 = first(zs)
  @assert all(==(z1), zs)
  # Remove one for site index
  return z1 - 1
end

childtagprefix() = "d="
childtags(n::Integer) = TagSet(string(childtagprefix(), n))
childtags(::Nothing) = ts""

generationtagprefix() = "g="
function generationtags(n1n2::Pair{<:Integer,<:Integer})
  n1, n2 = n1n2
  return TagSet(string(generationtagprefix(), n1, "↔", n2))
end
generationtags(::Nothing) = ts""

function ttntags(; c=nothing, d=nothing, g=nothing)
  return reduce(addtags, [celltags(c), childtags(d), generationtags(g)])
end

function init_space(::Type{Index{Int}}, linkdims::Integer)
  return linkdims
end

function init_space(::Type{ITensors.QNIndex}, linkdims::Integer)
  return (QN() => linkdims)
end

function InfiniteTTN(s::Vector{IndexT}, st; z::Integer, linkdims) where {IndexT<:Index}
  N = length(s)
  s∞ = infsiteinds(s)
  product_state = [state(s∞[n], st(n)) for n in 1:N]
  ttn = Vector{ITensor}(undef, N)
  link_space = init_space(IndexT, linkdims)
  l_parents = Vector{IndexT}(undef, N)
  l_children = Vector{Vector{IndexT}}(undef, N)
  n = 1
  if N == 1
    l_parents[n] = Index(link_space; tags=ttntags(; c=0, d=1, g=(mod1(n - 1, N) => n)))
    l_child = translatecell(l_parents[n], 1)
    l_children[n] = [replacetags(l_child, childtags(1) => childtags(d)) for d in 1:(z - 1)]
    ttn[n] = randomITensor(dag(l_parents[n]), l_children[n]) * product_state[n]
    return InfiniteTTN(ttn)
  end
  l_parents[n] = Index(link_space; tags=ttntags(; c=0, d=1, g=(mod1(n - 1, N) => n)))
  l_child = Index(link_space; tags=ttntags(; c=1, g=(n => mod1(n + 1, N))))
  l_children[n] = [addtags(l_child, childtags(d)) for d in 1:(z - 1)]
  ttn[n] = randomITensor(dag(l_parents[n]), l_children[n]) * product_state[n]
  for n in 2:(N - 1)
    l_parents[n] = first(l_children[n - 1])
    l_child = Index(link_space; tags=ttntags(; c=1, g=(n => mod1(n + 1, N))))
    l_children[n] = [addtags(l_child, childtags(d)) for d in 1:(z - 1)]
    ttn[n] = randomITensor(dag(l_parents[n]), l_children[n]) * product_state[n]
  end
  n = N
  l_parents[n] = first(l_children[n - 1])
  l_child = translatecell(l_parents[1], 1)
  l_children[n] = [replacetags(l_child, childtags(1) => childtags(d)) for d in 1:(z - 1)]
  ttn[n] = randomITensor(dag(l_parents[n]), l_children[n]) * product_state[n]
  return InfiniteTTN(ttn)
end

getchildtag(ts::TagSet) = tag_starting_with(ts, childtagprefix())
function getchild(ts::TagSet)
  childtag = getchildtag(ts)
  return parse(Int, childtag[(length(childtagprefix()) + 1):end])
end
getchildtags(ts::TagSet) = childtags(getchild(ts))
getchildtags(i::Index) = getchildtags(tags(i))

# Get the tensor at generation `g` and branch `b`
function Base.getindex(ψ::InfiniteTTN, g::Integer, b::Integer)
  # TODO: don't allow `g < 0`, make `g == 0` a special case
  @assert 1 ≤ b ≤ (degree(ψ) - 1)
  # Tensor at generation `g` on the first branch
  ψg1 = ψ.data[g]
  b == 1 && return ψg1
  l_parent_1 = commonind(ψ.data[g - 1], ψg1)
  l_parent_b = replacetags(l_parent_1, getchildtags(l_parent_1) => childtags(b))
  return replaceinds(ψg1, l_parent_1 => l_parent_b)
end

function ITensorInfiniteMPS.TransferMatrix(
  ψ::InfiniteTTN, g::Integer, d::Pair{<:Integer,<:Integer}; dir="in"
)
  @assert dir == "in"
  d1, d2 = d
  A = ψ[g, d2]
  l_in = commonind(A, ψ[g + 1, d1])
  l_out = commonind(A, ψ[g - 1, 1])
  Adag = prime(dag(A); inds=(l_in, l_out))
  return ITensorMap(
    [A, Adag]; input_inds=(dag(l_in), l_in'), output_inds=(l_out, dag(l_out)')
  )
end

end
