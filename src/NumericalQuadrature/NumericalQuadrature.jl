module NumericalQuadrature

using LinearAlgebra
using ..AtomicPotentials.MeshClassification

export QuadratureMethod
export Trapezoid
export Simpson
export AbinitCorrectedTrapezoid
export QESimpson
export integration_weights!
export integration_weights
export integrate!
export integrate

abstract type QuadratureMethod end

Base.Broadcast.broadcastable(m::QuadratureMethod) = Ref(m)

function integration_weights(x::AbstractVector, method::QuadratureMethod)
    weights = similar(x)
    return integration_weights!(weights, x, method)
end

function integrate!(
    weights::AbstractVector, x::AbstractVector, y::AbstractVector, method::QuadratureMethod
)
    integration_weights!(weights, x, method)
    return dot(weights, y)
end

function integrate(x::AbstractVector, y::AbstractVector, method::QuadratureMethod)
    weights = integration_weights(x, method)
    return dot(weights, y)
end

include("trapezoid.jl")
include("simpson.jl")
include("qe_simpson.jl")
include("abinit_corrected_trapezoid.jl")
end
