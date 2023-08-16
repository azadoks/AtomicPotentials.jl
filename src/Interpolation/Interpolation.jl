module Interpolation

using BSplineKit: BSplineKit
using Interpolations: Interpolations
using ..AtomicPotentials.MeshClassification

export InterpolationMethod
export Spline
export Linear
export evaluate
export resample

abstract type InterpolationMethod end

function resample(
    xp::V, x::AbstractVector, y::AbstractVector, method::InterpolationMethod; kwargs...
)::V where {V<:AbstractVector}
    T = eltype(x)
    yp = similar(xp)
    itp = interpolate(x, y, method; kwargs...)
    yp .= evaluate.(itp, xp)
    return yp
end

struct Linear <: InterpolationMethod end

function interpolate(x::AbstractVector, y::AbstractVector, ::Linear; kwargs...)
    if typeof(Mesh(x)) <: Uniform
        return Interpolations.scale(
            Interpolations.interpolate(y, Interpolations.BSpline(Interpolations.Linear())),
            x,
        )
    else
        return Interpolations.interpolate(
            (x,), y, Interpolations.Gridded(Interpolations.Linear())
        )
    end
end

function evaluate(itp::Interpolations.ScaledInterpolation, x::T)::T where {T<:Real}
    return itp(x)
end

struct Spline <: InterpolationMethod
    order::Int
    bc::Union{Nothing,BSplineKit.Natural}
end

function Spline(order::Integer; bc=BSplineKit.Natural())
    return Spline(order, bc)
end

function interpolate(x::AbstractVector, y::AbstractVector, method::Spline)
    return BSplineKit.interpolate(x, y, BSplineKit.BSplineOrder(method.order), method.bc)
end

function evaluate(itp::BSplineKit.SplineWrapper, x::T)::T where {T<:Real}
    return itp.spline(x)
end
end
