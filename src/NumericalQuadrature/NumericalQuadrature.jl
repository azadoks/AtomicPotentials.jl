module NumericalQuadrature

using LinearAlgebra
using Adapt
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
    # GPU is not supported for weight vector generation because scalar indexing is needed
    # We bring `x` to the CPU before generating the weights then put the weights on the
    # device after generation
    x = adapt(Array, x)  # x -> CPU
    weights = similar(x)
    integration_weights!(weights, x, method)
    return adapt(T, weights)  # weights -> GPU/CPU
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
