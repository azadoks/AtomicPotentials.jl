module Interpolation

using CubicSplines: CubicSplines
using Dierckx: Dierckx
using BSplineKit: BSplineKit
using Interpolations: Interpolations
using ..AtomicPotentials.MeshClassification

export InterpolationMethod
export Spline
export Linear
export CubicSpline
export DierckxSpline
export resample

abstract type InterpolationMethod end

function resample(
    xp::V, x::AbstractVector, y::AbstractVector, method::InterpolationMethod; kwargs...
)::V where {V<:AbstractVector}
    T = eltype(x)
    yp = similar(xp)
    itp = interpolate(x, y, method; kwargs...)
    yp .= itp.(xp)
    return yp
end

struct Linear <: InterpolationMethod end

function interpolate(x::AbstractVector, y::AbstractVector, ::Linear; kwargs...)
    # if typeof(Mesh(x)) <: Uniform
    if (x[2] - x[1]) ≈ (x[3] - x[2])
        return Interpolations.scale(
            Interpolations.interpolate(y, Interpolations.BSpline(Interpolations.Linear())),
            range(minimum(x), maximum(x), length(x)),
        )
    else
        return Interpolations.interpolate(
            (x,), y, Interpolations.Gridded(Interpolations.Linear())
        )
    end
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

struct CubicSpline <: InterpolationMethod end

function interpolate(x::AbstractVector, y::AbstractVector, ::CubicSpline)
    return CubicSplines.CubicSpline(x, y)
end

struct DierckxSpline <: InterpolationMethod
    k::Int
    bc::String
end

function DierckxSpline(; k=3, bc="nearest")
    return DierckxSpline(k, bc)
end

function interpolate(x::AbstractVector, y::AbstractVector, method::DierckxSpline)
    return Dierckx.Spline1D(x, y; k=method.k, bc=method.bc)
end

end
