module NumericalQuadrature

using LinearAlgebra

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

function is_uniform(x::AbstractVector)
    length(x) ≤ 2 && return true
    return (x[begin + 1] - x[begin]) ≈ (x[begin + 2] - x[begin + 1])
end

function fit_log_without_zero(x::AbstractVector)
    b = first(x)
    a = log(last(x) / first(x)) / (length(x) - 1)
    return (; a, b)
end

compute_log_without_zero(a::Real, b::Real, i::Integer) = b * exp(a * (i - 1))
compute_log_without_zero_deriv(a::Real, b::Real, i::Integer) = a * b * exp(a * (i - 1))

function is_log_without_zero(x::AbstractVector)
    first(x) ≈ 0 && return false
    a, b = fit_log_without_zero(x)
    return all(x .≈ compute_log_without_zero.(a, b, eachindex(x)))
end

function fit_log_with_zero(x::AbstractVector)
    x₂, x₃ = x[(begin + 1):(begin + 2)]
    b = x₂^2 / (x₃ - 2x₂)
    a = log((x₂ + b) / b)
    return (; a, b)
end

compute_log_with_zero(a::Real, b::Real, i::Integer) = b * (exp(a * (i - 1)) - 1)
function compute_log_with_zero_deriv(a::Real, b::Real, i::Integer)
    return a * b * (exp(a * (i - 1)) - 1) + a * b
end

function is_log_with_zero(x::AbstractVector)
    !(first(x) ≈ 0) && return false
    a, b = fit_log_with_zero(x)
    return all(x .≈ compute_log_with_zero.(a, b, eachindex(x)))
end

function mesh_derivative(x::AbstractVector)
    if is_uniform(x)
        return fill(x[begin + 1] - x[begin], length(x))
    elseif is_log_with_zero(x)
        a, b = fit_log_with_zero(x)
        return compute_log_with_zero_deriv.(a, b, eachindex(x))
    elseif is_log_without_zero(x)
        a, b = fit_log_without_zero(x)
        return compute_log_without_zero_deriv.(a, b, eachindex(x))
    else
        return error("Unknown mesh type")
    end
end

include("trapezoid.jl")
include("simpson.jl")
include("qe_simpson.jl")
include("abinit_corrected_trapezoid.jl")
end
