module Interpolation

using BSplineKit: BSplineKit
using Interpolations: Interpolations
using ..AtomicPotentials.MeshClassification

export InterpolationMethod
export InterpolationsLinear
export InterpolationsCubicSpline
export BSplineKitCubicSpline
export Interpolator
export resample

Interpolator = Union{<:Interpolations.ScaledInterpolation,<:BSplineKit.SplineInterpolation}

abstract type InterpolationMethod end

function resample(
    xp::V, x::AbstractVector, y::AbstractVector, method::InterpolationMethod; kwargs...
)::V where {V<:AbstractVector}
    T = eltype(x)
    yp = similar(xp, T)
    itp = interpolate(x, y, method; kwargs...)
    yp .= itp.(xp)
    return yp
end

struct InterpolationsLinear <: InterpolationMethod end

function interpolate(x::AbstractVector, y::AbstractVector, ::InterpolationsLinear)
    if std(diff(x)) ≈ 0
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

struct InterpolationsCubicSpline <: InterpolationMethod
    bc
    function InterpolationsCubicSpline(bc=Interpolations.Line(Interpolations.OnGrid()))
        return new(bc)
    end
end

function interpolate(
    x::AbstractVector, y::AbstractVector, method::InterpolationsCubicSpline
)
    return Interpolations.scale(
        Interpolations.interpolate(
            y, Interpolations.BSpline(Interpolations.Cubic(method.bc))
        ),
        range(minimum(x), maximum(x), length(x)),
    )
end

struct BSplineKitSpline <: InterpolationMethod
    bc
    function BSplineKitSpline(bc=BSplineKit.Natural())
        return new(bc)
    end
end

function interpolate(x::AbstractVector, y::AbstractVector, method::BSplineKitSpline)
    return BSplineKit.interpolate(x, y, BSplineKit.BSplineOrder(4), method.bc)
end

end
