module NumericalQuadrature

using LinearAlgebra

export QuadratureMethod
export Uniformity
export Uniform
export NonUniform
export QuadratureMethodOrType
export Trapezoid
export Simpson
export integration_weights!
export integration_weights
export integrate!
export integrate

abstract type Uniformity end
struct Uniform <: Uniformity end
struct NonUniform <: Uniformity end

abstract type QuadratureMethod end

QuadratureMethodOrType = Union{QuadratureMethod,Type{<:QuadratureMethod}}

function integration_weights!(
    weights::AbstractVector, x::AbstractVector, ::Type{M}
) where {M<:QuadratureMethod}
    length(x) <= 2 && return integration_weights!(weights, x, M{Uniform}())
    uniformity =
        (x[begin + 1] - x[begin]) ≈ (x[begin + 2] - x[begin + 1]) ? Uniform : NonUniform
    return integration_weights!(weights, x, M{uniformity}())
end

function integration_weights(x::AbstractVector, method::QuadratureMethodOrType)
    weights = similar(x)
    return integration_weights!(weights, x, method)
end

function integrate!(
    weights::AbstractVector,
    x::AbstractVector,
    y::AbstractVector,
    method::QuadratureMethodOrType,
)
    integration_weights!(weights, x, method)
    return dot(weights, y)
end

function integrate(x::AbstractVector, y::AbstractVector, method::QuadratureMethodOrType)
    weights = integration_weights(x, method)
    return dot(weights, y)
end

include("trapezoid.jl")
include("simpson.jl")
# include("qe_simpson.jl")
# include("abinit_corrected_trapezoid.jl")
end
